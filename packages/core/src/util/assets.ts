/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UUID } from './uuid.js';
import { fetchData, type DataType, type DataResponse, decodeData, readFromFile } from './data-source.js';
import { Task, type RuntimeContext } from '@molstar/core/task';
import { File_ as File } from './nodejs-shims.js';
import { iterableToArray } from '@molstar/core/data/util/array';
import { gzipDeclaredSize, isGzip, ungzip } from './zip/zip.js';

export { AssetManager, DefaultAssetManager, Asset };

type _File = File;
type Asset = Asset.Url | Asset.File;

namespace Asset {
  export type Url = {
    kind: 'url';
    id: UUID;
    url: string;
    title?: string;
    body?: string;
    headers?: Record<string, string>;
  };
  export type File = { kind: 'file'; id: UUID; name: string; file?: _File };

  export function Url(url: string, options?: { body?: string; title?: string; headers?: Record<string, string> }): Url {
    return { kind: 'url', id: UUID.create22(), url, ...options };
  }

  export function File(file: _File): File {
    return { kind: 'file', id: UUID.create22(), name: file.name, file };
  }

  export function isUrl(x?: Asset): x is Url {
    return x?.kind === 'url';
  }

  export function isFile(x?: Asset): x is File {
    return x?.kind === 'file';
  }

  export interface Wrapper<T extends DataType = DataType> {
    readonly data: DataResponse<T>;
    dispose: () => void;
  }

  export function Wrapper<T extends DataType = DataType>(data: DataResponse<T>, asset: Asset, manager: AssetManager) {
    return {
      data,
      dispose: () => {
        manager.release(asset);
      },
    };
  }

  export function getUrl(url: string | Url) {
    return typeof url === 'string' ? url : url.url;
  }

  export function getUrlAsset(
    manager: AssetManager,
    url: string | Url,
    body?: string,
    headers?: Record<string, string>,
  ) {
    if (typeof url === 'string') {
      const asset = manager.tryFindUrl(url, body, headers);
      return asset || Url(url, { body, headers });
    }
    return url;
  }
}

function urlHeadersEqual(a?: Record<string, string>, b?: Record<string, string>) {
  const aKeys = a ? Object.keys(a) : [];
  const bKeys = b ? Object.keys(b) : [];
  if (aKeys.length !== bKeys.length) return false;

  for (const key of aKeys) {
    if (a![key] !== b![key]) return false;
  }
  return true;
}

/**
 * The key of the byte cache (`DefaultAssetManager`), or `undefined` for a request that is not cached.
 *
 * Only plain GETs are cached: a request with a body (POST: GraphQL, queries) can answer differently every time, so it
 * keeps the per-asset behavior (one download per `Asset.Url`). The headers are part of the key (an `Authorization` or
 * `Accept` header changes the response); the title is not. The response type is deliberately not part of the key: the
 * cache holds the raw bytes and every type (`string`, `json`, `binary`, ...) is decoded from them (`decodeData`), so
 * a text and a binary read of one URL share one download.
 */
function urlCacheKey(asset: Asset.Url): string | undefined {
  if (asset.body) return undefined;
  const headers = asset.headers
    ? Object.keys(asset.headers)
        .sort()
        .map((k) => [k, asset.headers![k]])
    : [];
  return JSON.stringify(['GET', asset.url, headers]);
}

/** A cached download: in flight (`pending`) or done (`bytes`); `refCount` counts holders, in-flight waiters included. */
interface UrlCacheEntry extends InflatedCache {
  refCount: number;
  bytes?: Uint8Array<ArrayBuffer>;
  pending?: Promise<Uint8Array<ArrayBuffer>>;
}

/** The inflated bytes of gzipped held bytes, made on the first `decompress` read (`util/_spec/assets.md` §5). */
interface InflatedCache {
  inflated?: Promise<Uint8Array<ArrayBuffer>>;
}

/**
 * Deflate shrinks data at most about 1032 times: a gzip that declares (ISIZE) more than that of its size is corrupt or
 * crafted, and refused before its buffer is allocated (`util/_spec/assets.md` §5).
 */
const MaxDeflateRatio = 1032;

function inflateGzip(ctx: RuntimeContext, data: Uint8Array<ArrayBuffer>) {
  const size = gzipDeclaredSize(data);
  if (size > data.length * MaxDeflateRatio + 1024) {
    throw new Error(`gzip data declares ${size} bytes inflated, more than ${data.length} compressed bytes can hold`);
  }
  return ungzip(ctx, data);
}

/** Whether the task of `ctx` was aborted: its next update throws an abort. */
async function isAborted(ctx: RuntimeContext) {
  try {
    await ctx.update();
    return false;
  } catch (e) {
    return Task.isAbort(e);
  }
}

/**
 * A held asset. Its bytes are held once (`util/_spec/assets.md`): a download's bytes, or a file's read (and
 * decompressed) on its first shared read. `file` is the source file, or for a download one made from the bytes on
 * access.
 */
class HeldAsset implements AssetManager.Entry {
  refCount = 0;
  /** The decompressed bytes every read but `'zip'` decodes; `undefined` until a file's first shared read. */
  bytes: Uint8Array<ArrayBuffer> | undefined;

  get file(): File {
    return this.source ?? new File([this.bytes!], 'raw-data');
  }

  get name() {
    return this.source?.name ?? 'raw-data';
  }

  /** Lets go of bytes read from the source file; a download's bytes are all it has and stay. */
  dropBytes() {
    if (!this.source) return;
    this.bytes = undefined;
    this.inflatedCache.inflated = undefined;
  }

  /**
   * Holds no bytes in memory: a file keeps its source; a download's bytes become one (a `Blob`, which a browser may keep
   * outside the JS heap), read again from it on the next read (§6).
   */
  toSource() {
    if (!this.source) this.source = new File([this.bytes!], 'raw-data');
    this.bytes = undefined;
    this.inflatedCache.inflated = undefined;
  }

  constructor(
    readonly asset: Asset,
    /** The file the bytes are read from; `undefined` for a download, until `toSource`. */
    private source: File | undefined,
    bytes?: Uint8Array<ArrayBuffer>,
    /** The URL cache entry it holds; `undefined` once it has let go of it (`letGoOfBytes`). */
    public cacheKey?: string,
    readonly isStatic?: boolean,
    readonly tag?: string,
    /** A download's is its URL cache entry's, shared by every asset of the URL. */
    private readonly inflatedCache: InflatedCache = {},
  ) {
    this.bytes = bytes;
  }

  async read<T extends DataType>(
    ctx: RuntimeContext,
    type: T,
    shared: boolean,
    decompress: boolean,
  ): Promise<DataResponse<T>> {
    // A file is read from its source unless its bytes are held or asked for shared: holding them for every read would
    // keep in memory what the file keeps on disk. A zip keeps its entries apart, which `bytes` (one decompressed entry)
    // cannot stand for.
    const holds = this.bytes !== undefined || (shared && type === 'binary');
    if (this.source && (type === 'zip' || (!holds && !decompress))) {
      return readFromFile(this.source, type).runInContext(ctx);
    }
    // `owned`: bytes read or inflated for this read only, nobody else's: a binary read returns them without a copy
    let owned = false;
    let bytes = this.bytes;
    if (!bytes) {
      bytes = await readFromFile(this.source!, 'binary').runInContext(ctx);
      if (holds) this.bytes = bytes;
      else owned = true;
    }
    if (decompress && type !== 'zip' && isGzip(bytes)) {
      bytes = await this.inflate(ctx, bytes, holds);
      owned = !holds;
    }
    if (type === 'binary') return (shared || owned ? bytes : bytes.slice()) as DataResponse<T>;
    return decodeData(ctx, bytes, type);
  }

  /**
   * Held bytes inflate once, whoever asks first; a failure is not kept. The inflate runs in the first reader's task:
   * when that task is aborted, a reader that was not takes over, as a waiter of an aborted download does.
   */
  private async inflate(ctx: RuntimeContext, bytes: Uint8Array<ArrayBuffer>, held: boolean) {
    if (!held) return inflateGzip(ctx, bytes);
    const cache = this.inflatedCache;
    for (;;) {
      let inflated = cache.inflated;
      if (!inflated) {
        const started = inflateGzip(ctx, bytes);
        cache.inflated = inflated = started;
        started.catch(() => {
          if (cache.inflated === started) cache.inflated = undefined;
        });
      }
      try {
        return await inflated;
      } catch (e) {
        if (Task.isAbort(e) && !(await isAborted(ctx))) continue;
        throw e;
      }
    }
  }
}

/**
 * Holds assets (URL downloads and files) and their bytes, ref-counted: `resolve` takes a reference and
 * `Asset.Wrapper.dispose` (`release`) gives it back. `DefaultAssetManager` is the implementation; a plugin takes a
 * custom one from `PluginSpec.assetManager`.
 */
interface AssetManager {
  /** Every held asset with its bytes. */
  readonly assets: AssetManager.Entry[];
  tryFindUrl(url: string, body?: string, headers?: Record<string, string>): Asset.Url | undefined;
  tryFindFilename(name: string): Asset | undefined;
  /** Holds `file` as the bytes of `asset`, replacing what it held; a static asset outlives its last `release`. */
  set(asset: Asset, file: File, options?: { isStatic?: boolean; tag?: string }): void;
  get(asset: Asset): AssetManager.Entry | undefined;
  delete(asset: Asset): boolean;
  has(asset: Asset): boolean;
  /**
   * Reads the asset, taking a reference (`store` on, the default) that the returned wrapper's `dispose` gives back.
   * `shared`: a `'binary'` read may return the held bytes themselves, which the caller must not change (a manager may
   * still return a copy).
   */
  resolve<T extends DataType>(
    asset: Asset,
    type: T,
    store?: boolean | AssetManager.ResolveOptions,
  ): Task<Asset.Wrapper<T>>;
  release(asset: Asset): void;
  /** `resolve` honors `decompress` (`DefaultAssetManager` does); without it `AssetManager.resolveData` inflates. */
  readonly decompresses?: boolean;
  /**
   * The asset's bytes are not needed in memory any more (its reader holds them in another form, a response turned to
   * text): a manager may let go of them, keeping what it reads them again from. Optional, a hint.
   */
  letGoOfBytes?(asset: Asset): void;
  clearTag(tag: string): void;
  clear(): void;
  dispose(): void;
}

namespace AssetManager {
  /**
   * `manager.resolve` with `decompress` honored by any manager: one that does not inflate (`decompresses` unset, a
   * custom one) is asked for the bytes, which are inflated, uncached, and read as `type` here (`util/_spec/assets.md`
   * §5).
   */
  export function resolveData<T extends DataType>(
    manager: AssetManager,
    asset: Asset,
    type: T,
    options: ResolveOptions = {},
  ): Task<Asset.Wrapper<T>> {
    if (!options.decompress || manager.decompresses) return manager.resolve(asset, type, options);
    return Task.create(`Read ${Asset.isUrl(asset) ? asset.url : asset.name}`, async (ctx) => {
      const held = await manager.resolve(asset, 'binary', { ...options, decompress: false }).runInContext(ctx);
      try {
        let bytes = held.data;
        if (type !== 'zip' && isGzip(bytes)) bytes = await inflateGzip(ctx, bytes);
        const data = (type === 'binary' ? bytes : await decodeData(ctx, bytes, type)) as DataResponse<T>;
        return { data, dispose: held.dispose };
      } catch (e) {
        held.dispose();
        throw e;
      }
    });
  }

  export interface ResolveOptions {
    /** Hold the asset (default `true`). */
    store?: boolean;
    /** A `'binary'` read may return the held bytes, read-only, instead of a copy (default `false`). */
    shared?: boolean;
    /**
     * Gzip bytes (told by their magic bytes, not the name) are inflated before the read, once for held bytes (default
     * `false`): a server that sends `.gz` files without `Content-Encoding`. Not for `'zip'` reads.
     */
    decompress?: boolean;
  }

  export interface Entry {
    readonly asset: Asset;
    readonly file: File;
    readonly isStatic?: boolean;
    readonly tag?: string;
  }
}

/**
 * The default `AssetManager`.
 *
 * URL assets read through a URL-keyed, ref-counted byte cache (`urlCacheKey`): the first `resolve` of a URL downloads,
 * later ones (another `Asset.Url` of the same URL and headers, or the same one) are served from the held bytes, and
 * concurrent ones share the download in flight, so a URL is fetched once whoever asks first. The bytes live while any
 * asset that resolved to them is held and are dropped with the last one.
 *
 * - A failed download is not cached, the next `resolve` fetches again. The waiters of that download get its error.
 * - An aborted download (its requester was aborted) does not fail the waiters: they that are still interested retry,
 *   one of them becoming the downloader. A waiter itself notices its own abort when its task next updates.
 * - `resolve(asset, type, false)` (`store` off) and requests with a body bypass the cache.
 *
 * Bytes are held once (`util/_spec/assets.md`): a download as received, a file read and decompressed on its first
 * `shared` read; every read decodes from them, and a `shared` binary read returns them without a copy.
 */
class DefaultAssetManager implements AssetManager {
  // TODO: when serializing, check for duplicates?
  readonly decompresses = true;

  private _assets = new Map<string, HeldAsset>();
  private _urlCache = new Map<string, UrlCacheEntry>();

  get assets(): AssetManager.Entry[] {
    return iterableToArray(this._assets.values());
  }

  /** The number of URLs the byte cache holds or is downloading. */
  get urlCacheSize() {
    return this._urlCache.size;
  }

  tryFindUrl(url: string, body?: string, headers?: Record<string, string>): Asset.Url | undefined {
    const assets = this.assets.values();
    while (true) {
      const v = assets.next();
      if (v.done) return;
      const asset = v.value.asset;
      if (
        Asset.isUrl(asset) &&
        asset.url === url &&
        (asset.body || '') === (body || '') &&
        urlHeadersEqual(asset.headers, headers)
      )
        return asset;
    }
  }

  tryFindFilename(name: string): Asset | undefined {
    const it = this._assets.values();
    while (true) {
      const { done, value } = it.next();
      if (done) break;
      if (value.name === name) return value.asset;
    }
  }

  set(asset: Asset, file: File, options?: { isStatic?: boolean; tag?: string }) {
    this._remove(asset.id);
    this._assets.set(asset.id, new HeldAsset(asset, file, undefined, undefined, options?.isStatic, options?.tag));
  }

  get(asset: Asset): AssetManager.Entry | undefined {
    return this._assets.get(asset.id);
  }

  delete(asset: Asset) {
    return this._remove(asset.id);
  }

  has(asset: Asset) {
    return this._assets.has(asset.id);
  }

  resolve<T extends DataType>(
    asset: Asset,
    type: T,
    options: boolean | AssetManager.ResolveOptions = true,
  ): Task<Asset.Wrapper<T>> {
    const store = typeof options === 'boolean' ? options : options.store !== false;
    const shared = typeof options === 'boolean' ? false : !!options.shared;
    const decompress = typeof options === 'boolean' ? false : !!options.decompress;
    if (Asset.isUrl(asset)) {
      return Task.create(`Download ${asset.title || asset.url}`, async (ctx) => {
        if (this._assets.has(asset.id)) {
          const entry = this._assets.get(asset.id)!;
          entry.refCount++;
          return this._readHeld(ctx, entry, type, shared, decompress);
        }

        if (!store) {
          if (!decompress) return Asset.Wrapper(await fetchData({ ...asset, type }).runInContext(ctx), asset, this);
          const entry = new HeldAsset(
            asset,
            undefined,
            await fetchData({ ...asset, type: 'binary' }).runInContext(ctx),
          );
          return Asset.Wrapper(await entry.read(ctx, type, true, true), asset, this);
        }

        const key = urlCacheKey(asset);
        if (key === undefined) {
          const data = await fetchData({ ...asset, type: 'binary' }).runInContext(ctx);
          const entry = new HeldAsset(asset, undefined, data);
          this._assets.set(asset.id, entry);
          entry.refCount = 1;
          return this._readHeld(ctx, entry, type, shared, decompress);
        }

        // read through the byte cache: the asset takes one hold of the entry, given back when it leaves the manager
        const held = await this._acquireUrl(ctx, asset, key);
        let entry = this._assets.get(asset.id);
        if (entry) {
          // the same asset was resolved concurrently and registered first
          this._releaseUrl(key);
        } else {
          entry = new HeldAsset(asset, undefined, held.bytes!, key, undefined, undefined, held);
          this._assets.set(asset.id, entry);
        }
        entry.refCount++;
        return this._readHeld(ctx, entry, type, shared, decompress);
      });
    } else {
      return Task.create(`Read ${asset.name}`, async (ctx) => {
        if (this._assets.has(asset.id)) {
          const entry = this._assets.get(asset.id)!;
          entry.refCount++;
          return this._readHeld(ctx, entry, type, shared, decompress);
        }
        if (!(asset.file instanceof File)) {
          throw new Error(`Cannot resolve file asset '${asset.name}' (${asset.id})`);
        }
        if (!store && !decompress) {
          return Asset.Wrapper(await readFromFile(asset.file, type).runInContext(ctx), asset, this);
        }
        if (!store)
          return Asset.Wrapper(await new HeldAsset(asset, asset.file).read(ctx, type, false, true), asset, this);
        const entry = new HeldAsset(asset, asset.file);
        this._assets.set(asset.id, entry);
        entry.refCount = 1;
        return this._readHeld(ctx, entry, type, shared, decompress);
      });
    }
  }

  /** Reads an entry the caller took a hold of; a failed read gives the hold back (perf-audit-2 01/B12) */
  private async _readHeld<T extends DataType>(
    ctx: RuntimeContext,
    entry: HeldAsset,
    type: T,
    shared: boolean,
    decompress: boolean,
  ): Promise<Asset.Wrapper<T>> {
    try {
      return Asset.Wrapper(await entry.read(ctx, type, shared, decompress), entry.asset, this);
    } catch (e) {
      this.release(entry.asset);
      throw e;
    }
  }

  /** §6: the asset's entry holds no bytes in memory; not while other assets of its URL share them. */
  letGoOfBytes(asset: Asset) {
    const entry = this._assets.get(asset.id);
    if (!entry) return;
    if (entry.cacheKey !== undefined) {
      const url = this._urlCache.get(entry.cacheKey);
      if (url && url.refCount > 1) return;
      // the bytes leave the URL cache with its only hold, this entry's: the next download of the URL fetches again,
      // and a release of this entry does not reach the new one
      if (url) this._urlCache.delete(entry.cacheKey);
      entry.cacheKey = undefined;
    }
    entry.toSource();
  }

  release(asset: Asset) {
    const entry = this._assets.get(asset.id);
    if (!entry) return;
    entry.refCount--;
    if (entry.refCount > 0) return;
    if (!entry.isStatic) this._remove(asset.id);
    // a static file outlives its holders; its bytes do not, they are read again on the next resolve
    else entry.dropBytes();
  }

  clearTag(tag: string) {
    const keys = Array.from(this._assets.keys());
    for (const key of keys) {
      const entry = this._assets.get(key);
      if (entry && entry.tag === tag) {
        this._remove(key);
      }
    }
  }

  clear() {
    for (const id of Array.from(this._assets.keys())) this._remove(id);
  }

  dispose() {
    this.clear();
  }

  /** Deletes an asset entry and gives back its hold of the byte cache. */
  private _remove(id: string) {
    const entry = this._assets.get(id);
    if (!entry) return false;
    this._assets.delete(id);
    if (entry.cacheKey !== undefined) this._releaseUrl(entry.cacheKey);
    return true;
  }

  private _releaseUrl(key: string, entry = this._urlCache.get(key)) {
    if (!entry) return;
    entry.refCount--;
    if (entry.refCount <= 0 && this._urlCache.get(key) === entry) this._urlCache.delete(key);
  }

  /** A hit, the download in flight, or a new download; returns an entry with `bytes` and one hold taken for the caller. */
  private async _acquireUrl(ctx: RuntimeContext, asset: Asset.Url, key: string): Promise<UrlCacheEntry> {
    for (;;) {
      const entry = this._urlCache.get(key);
      if (!entry) return this._download(ctx, asset, key);

      entry.refCount++;
      if (entry.bytes) return entry;

      // wait for the download in flight (an abort of this task surfaces here)
      try {
        await ctx.update({ message: `Waiting for ${asset.title || asset.url}` });
        await entry.pending;
        return entry;
      } catch (e) {
        this._releaseUrl(key, entry);
        // the downloader was aborted, not this task: not a failure for us, take over
        if (Task.isAbort(e) && !(await isAborted(ctx))) continue;
        throw e;
      }
    }
  }

  private async _download(ctx: RuntimeContext, asset: Asset.Url, key: string): Promise<UrlCacheEntry> {
    const entry: UrlCacheEntry = { refCount: 1 };
    this._urlCache.set(key, entry);
    const pending = fetchData({ ...asset, type: 'binary' }).runInContext(ctx);
    entry.pending = pending;
    pending.then(
      (bytes) => {
        entry.bytes = bytes;
        entry.pending = undefined;
      },
      () => {
        // not cached: the next request fetches again
        if (this._urlCache.get(key) === entry) this._urlCache.delete(key);
      },
    );
    try {
      await pending;
      return entry;
    } catch (e) {
      this._releaseUrl(key, entry);
      throw e;
    }
  }
}
