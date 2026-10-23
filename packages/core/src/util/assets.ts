/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UUID } from './uuid.js';
import { ajaxGet, type DataType, type DataResponse, readFromFile } from './data-source.js';
import { Task, type RuntimeContext } from '@molstar/core/task';
import { File_ as File } from './nodejs-shims.js';
import { iterableToArray } from '@molstar/core/data/util/array';

export { type AssetManager, DefaultAssetManager, Asset };

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
 * cache holds the raw bytes and every type (`string`, `json`, `binary`, ...) is decoded from them (`readFromFile`),
 * so a text and a binary read of one URL share one download.
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

/** A cached download: in flight (`pending`) or done (`file`); `refCount` counts holders, in-flight waiters included. */
interface UrlCacheEntry {
  refCount: number;
  file?: File;
  pending?: Promise<File>;
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
  /** Reads the asset, taking a reference (`store` on) that the returned wrapper's `dispose` gives back. */
  resolve<T extends DataType>(asset: Asset, type: T, store?: boolean): Task<Asset.Wrapper<T>>;
  release(asset: Asset): void;
  clearTag(tag: string): void;
  clear(): void;
  dispose(): void;
}

namespace AssetManager {
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
 */
class DefaultAssetManager implements AssetManager {
  // TODO: when serializing, check for duplicates?

  private _assets = new Map<
    string,
    { asset: Asset; file: File; refCount: number; isStatic?: boolean; tag?: string; cacheKey?: string }
  >();
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
      if (value.file.name === name) return value.asset;
    }
  }

  set(asset: Asset, file: File, options?: { isStatic?: boolean; tag?: string }) {
    this._remove(asset.id);
    this._assets.set(asset.id, { asset, file, refCount: 0, tag: options?.tag, isStatic: options?.isStatic });
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

  resolve<T extends DataType>(asset: Asset, type: T, store = true): Task<Asset.Wrapper<T>> {
    if (Asset.isUrl(asset)) {
      return Task.create(`Download ${asset.title || asset.url}`, async (ctx) => {
        if (this._assets.has(asset.id)) {
          const entry = this._assets.get(asset.id)!;
          entry.refCount++;
          return this._readHeld(ctx, asset, entry.file, type);
        }

        if (!store) {
          return Asset.Wrapper(await ajaxGet({ ...asset, type }).runInContext(ctx), asset, this);
        }

        const key = urlCacheKey(asset);
        if (key === undefined) {
          const data = await ajaxGet({ ...asset, type: 'binary' }).runInContext(ctx);
          const file = new File([data], 'raw-data');
          this._assets.set(asset.id, { asset, file, refCount: 1 });
          return this._readHeld(ctx, asset, file, type);
        }

        // read through the byte cache: the asset takes one hold of the entry, given back when it leaves the manager
        const held = await this._acquireUrl(ctx, asset, key);
        const existing = this._assets.get(asset.id);
        if (existing) {
          // the same asset was resolved concurrently and registered first
          this._releaseUrl(key);
          existing.refCount++;
        } else {
          this._assets.set(asset.id, { asset, file: held.file!, refCount: 1, cacheKey: key });
        }
        const file = (existing ?? this._assets.get(asset.id)!).file;
        return this._readHeld(ctx, asset, file, type);
      });
    } else {
      return Task.create(`Read ${asset.name}`, async (ctx) => {
        if (this._assets.has(asset.id)) {
          const entry = this._assets.get(asset.id)!;
          entry.refCount++;
          return this._readHeld(ctx, asset, entry.file, type);
        }
        if (!(asset.file instanceof File)) {
          throw new Error(`Cannot resolve file asset '${asset.name}' (${asset.id})`);
        }
        if (!store) return Asset.Wrapper(await readFromFile(asset.file, type).runInContext(ctx), asset, this);
        this._assets.set(asset.id, { asset, file: asset.file, refCount: 1 });
        return this._readHeld(ctx, asset, asset.file, type);
      });
    }
  }

  /** Reads the file of an entry the caller took a hold of; a failed read gives the hold back (perf-audit-2 01/B12) */
  private async _readHeld<T extends DataType>(
    ctx: RuntimeContext,
    asset: Asset,
    file: File,
    type: T,
  ): Promise<Asset.Wrapper<T>> {
    try {
      return Asset.Wrapper(await readFromFile(file, type).runInContext(ctx), asset, this);
    } catch (e) {
      this.release(asset);
      throw e;
    }
  }

  release(asset: Asset) {
    const entry = this._assets.get(asset.id);
    if (!entry) return;
    entry.refCount--;
    if (entry.refCount <= 0 && !entry.isStatic) this._remove(asset.id);
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

  /** A hit, the download in flight, or a new download; returns an entry with `file` and one hold taken for the caller. */
  private async _acquireUrl(ctx: RuntimeContext, asset: Asset.Url, key: string): Promise<UrlCacheEntry> {
    for (;;) {
      const entry = this._urlCache.get(key);
      if (!entry) return this._download(ctx, asset, key);

      entry.refCount++;
      if (entry.file) return entry;

      // wait for the download in flight (an abort of this task surfaces here)
      try {
        await ctx.update({ message: `Waiting for ${asset.title || asset.url}` });
        await entry.pending;
        return entry;
      } catch (e) {
        this._releaseUrl(key, entry);
        // the downloader was aborted, not this task: not a failure for us, take over
        if (Task.isAbort(e) && !(await this._isAborted(ctx))) continue;
        throw e;
      }
    }
  }

  private async _download(ctx: RuntimeContext, asset: Asset.Url, key: string): Promise<UrlCacheEntry> {
    const entry: UrlCacheEntry = { refCount: 1 };
    this._urlCache.set(key, entry);
    const pending = (async () => {
      const data = await ajaxGet({ ...asset, type: 'binary' }).runInContext(ctx);
      return new File([data], 'raw-data');
    })();
    entry.pending = pending;
    pending.then(
      (file) => {
        entry.file = file;
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

  private async _isAborted(ctx: RuntimeContext) {
    try {
      await ctx.update();
      return false;
    } catch (e) {
      return Task.isAbort(e);
    }
  }
}
