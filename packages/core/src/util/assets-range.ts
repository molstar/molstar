/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '@molstar/core/task';
import type { Asset, AssetManager } from './assets.js';
import { checkMaxBytes, getFS, SizeLimitExceeded } from './data-source.js';
import { RUNNING_IN_NODEJS } from './nodejs-shims.js';

// Range reads of an asset (`_spec/assets.md` §10): byte ranges of a URL or a file, read by page and kept in a page
// cache by asset and page, shared by the readers of the asset and let go with the last of them. Exported through
// `assets.ts`; `DefaultAssetManager.openRange` is `defaultRanges`, and `openRanges` serves a manager without one.

/** A byte range: `[offset, offset + length)`. */
export interface AssetRange {
  readonly offset: number;
  readonly length: number;
}

export interface AssetRangeReader {
  /** The asset's size in bytes, known after open (`Content-Range` total, `stat`, `Blob.size`). */
  readonly size: number;
  /**
   * One result per range, in order, each a copy-free view of cached pages where one page covers it. Ranges past the
   * end are clipped; an empty result is a zero-length array. Rejects as aborted with the signal's reason.
   */
  read(ranges: readonly AssetRange[], signal?: AbortSignal): Promise<Uint8Array[]>;
  /** Gives back the reference: pages are let go with the last holder. */
  dispose(): void;
}

export interface RangeReadOptions
  extends Pick<AssetManager.ResolveOptions, 'credentials' | 'redirect' | 'signal' | 'maxBytes'> {
  /** The page size, bytes (default 64 KiB; a power of two, 4 KiB to 4 MiB). */
  pageSize?: number;
  /** Pages held at most for this asset (default 256: 16 MiB); least recently used go first. */
  maxPages?: number;
}

/**
 * A range read cannot be trusted: a `Content-Range` that is not the one asked for or has no total, a 416 on a valid
 * range, the asset's size, `ETag` or `Last-Modified` changed between pages, a body longer or shorter than its range, a
 * reader used after its manager was cleared. Matched by `name` across copies of the core, as `SizeLimitExceeded` is.
 */
export class AssetRangeError extends Error {
  constructor(
    why: string,
    /** The URL or the file name. */
    readonly source?: string,
  ) {
    super(source ? `${source}: ${why}` : why);
    this.name = 'AssetRangeError';
  }
}

const DefaultPageSize = 64 * 1024;
const DefaultMaxPages = 256;
/** Requests of one asset in flight at once. */
const Concurrency = 4;

/**
 * `manager.openRange` when the manager has it; else the whole asset is resolved (`binary`, `shared`, bounded by
 * `maxBytes`) and ranges are served as views of its bytes: correct, not lazy.
 */
export function openRanges(
  manager: AssetManager,
  asset: Asset,
  options: RangeReadOptions = {},
): Task<AssetRangeReader> {
  if (manager.openRange) return manager.openRange(asset, options);
  return wholeReader(manager, asset, options);
}

function wholeReader(manager: AssetManager, asset: Asset, options: RangeReadOptions): Task<AssetRangeReader> {
  return Task.create(`Open ${nameOf(asset)}`, async (ctx) => {
    const { credentials, redirect, signal, maxBytes } = options;
    const held = await manager
      .resolve(asset, 'binary', {
        shared: true,
        ...(credentials !== undefined ? { credentials } : {}),
        ...(redirect !== undefined ? { redirect } : {}),
        ...(signal !== undefined ? { signal } : {}),
        ...(maxBytes !== undefined ? { maxBytes } : {}),
      })
      .runInContext(ctx);
    return bytesReader(held.data, nameOf(asset), held.dispose);
  });
}

const nameOf = (asset: Asset) => (asset.kind === 'url' ? asset.url : asset.name);

function abortedBy(signal: AbortSignal) {
  const { reason } = signal;
  return Task.Aborted(typeof reason === 'string' ? reason : reason instanceof Error ? reason.message : '');
}

/** `[a, b)` of a range clipped to `size`; a range that is not one (negative, not an integer) throws. */
function clip(r: AssetRange, size: number): [number, number] {
  if (!Number.isSafeInteger(r.offset) || !Number.isSafeInteger(r.length) || r.offset < 0 || r.length < 0)
    throw new RangeError(`not a byte range: ${r.offset} + ${r.length}`);
  const a = Math.min(r.offset, size);
  return [a, Math.min(size, r.offset + r.length)];
}

/** A reader of bytes held whole: every range a view of them. */
function bytesReader(bytes: Uint8Array, name: string, release: () => void): AssetRangeReader {
  let disposed = false;
  return {
    size: bytes.length,
    async read(ranges, signal) {
      if (signal?.aborted) throw abortedBy(signal);
      if (disposed) throw new AssetRangeError('the reader was disposed', name);
      return ranges.map((r) => {
        const [a, b] = clip(r, bytes.length);
        return bytes.subarray(a, b);
      });
    },
    dispose() {
      if (disposed) return;
      disposed = true;
      release();
    },
  };
}

/**
 * Work shared by readers: started once, each reader waits on it with its own signal; when every waiter has left, the
 * work is aborted (a page nobody wants is not fetched to the end) and, failed, it is not kept.
 */
class Shared<T> {
  readonly promise: Promise<T>;
  readonly controller = new AbortController();
  private waiters = 0;
  settled = false;

  constructor(start: (signal: AbortSignal) => Promise<T>) {
    this.promise = start(this.controller.signal);
    this.promise.then(
      () => (this.settled = true),
      () => (this.settled = true),
    );
  }

  /** The result, or the abort of `signal` as soon as it comes (leaving: the last to leave aborts the work). */
  join(signal: AbortSignal | undefined): Promise<T> {
    this.waiters++;
    return new Promise<T>((resolve, reject) => {
      let left = false;
      const leave = () => {
        if (left) return false;
        left = true;
        this.waiters--;
        signal?.removeEventListener('abort', onAbort);
        return true;
      };
      const onAbort = () => {
        if (!leave()) return;
        if (this.waiters === 0 && !this.settled) this.controller.abort(signal!.reason);
        reject(abortedBy(signal!));
      };
      if (signal?.aborted) return onAbort();
      signal?.addEventListener('abort', onAbort, { once: true });
      this.promise.then(
        (v) => {
          if (leave()) resolve(v);
        },
        (e) => {
          if (leave()) reject(e);
        },
      );
    });
  }
}

/** Where a store's pages come from: `[a, b)` of the asset, exactly those bytes. */
type FetchBytes = (a: number, b: number, signal: AbortSignal) => Promise<Uint8Array>;

/**
 * The pages of one asset (one key: a URL with its options, or a file), held by its readers. `whole` when the asset's
 * bytes are all held (a server that ignores `Range`, a file the manager holds the bytes of): ranges are its views.
 */
class PageStore {
  holders = 0;
  size = 0;
  whole: Uint8Array | undefined;
  closed: Error | undefined;
  /** Least recently used first (a page read again moves to the end). */
  private readonly pages = new Map<number, Uint8Array>();
  private readonly pending = new Map<number, { run: Shared<Uint8Array>; first: number }>();
  private active = 0;
  private readonly queue: (() => void)[] = [];

  constructor(
    readonly name: string,
    readonly pageSize: number,
    readonly maxPages: number,
    readonly fetch: FetchBytes,
    private readonly onEmpty: () => void,
  ) {}

  /** Pages of `[a, b)` read before open finished (the probe's): kept as fetched pages are. */
  keep(first: number, bytes: Uint8Array) {
    const ps = this.pageSize;
    for (let i = 0, p = first; i * ps < bytes.length; i++, p++)
      this.put(p, bytes.subarray(i * ps, Math.min(bytes.length, (i + 1) * ps)));
  }

  private put(page: number, bytes: Uint8Array) {
    this.pages.delete(page);
    this.pages.set(page, bytes);
    while (this.pages.size > this.maxPages) this.pages.delete(this.pages.keys().next().value!);
  }

  private hit(page: number): Uint8Array | undefined {
    const bytes = this.pages.get(page);
    if (bytes) {
      this.pages.delete(page);
      this.pages.set(page, bytes);
    }
    return bytes;
  }

  /** A request slot (at most `Concurrency` at once), released by the returned function. */
  private async slot(signal: AbortSignal): Promise<() => void> {
    while (this.active >= Concurrency) {
      await new Promise<void>((resolve) => this.queue.push(resolve));
      if (signal.aborted) {
        // the wake-up was for a request that is gone: pass it on
        this.queue.shift()?.();
        throw abortedBy(signal);
      }
    }
    this.active++;
    return () => {
      this.active--;
      this.queue.shift()?.();
    };
  }

  /** A run of missing pages, fetched in one request; kept page by page when it arrives, forgotten when it fails. */
  private startRun(first: number, count: number): Shared<Uint8Array> {
    const ps = this.pageSize;
    const run = new Shared(async (signal) => {
      const release = await this.slot(signal);
      try {
        const a = first * ps;
        const b = Math.min(this.size, (first + count) * ps);
        const bytes = await this.fetch(a, b, signal);
        if (bytes.length !== b - a)
          throw new AssetRangeError(`${bytes.length} bytes for a range of ${b - a}`, this.name);
        return bytes;
      } finally {
        release();
      }
    });
    for (let p = first; p < first + count; p++) this.pending.set(p, { run, first });
    const forget = () => {
      for (let p = first; p < first + count; p++) if (this.pending.get(p)?.run === run) this.pending.delete(p);
    };
    run.promise.then((bytes) => {
      forget();
      if (!this.closed) this.keep(first, bytes);
    }, forget);
    return run;
  }

  async read(ranges: readonly AssetRange[], signal?: AbortSignal): Promise<Uint8Array[]> {
    if (signal?.aborted) throw abortedBy(signal);
    if (this.closed) throw this.closed;
    const spans = ranges.map((r) => clip(r, this.size));
    if (this.whole) {
      const whole = this.whole;
      return spans.map(([a, b]) => whole.subarray(a, b));
    }
    const ps = this.pageSize;
    // the pages this read needs: held ones taken now (a later eviction does not take them from this read)
    const got = new Map<number, Uint8Array>();
    const runs = new Map<Shared<Uint8Array>, number>();
    const missing: number[] = [];
    for (const [a, b] of spans) {
      if (a === b) continue;
      for (let p = Math.floor(a / ps); p <= Math.floor((b - 1) / ps); p++) {
        if (got.has(p)) continue;
        const held = this.hit(p);
        if (held) {
          got.set(p, held);
          continue;
        }
        const pending = this.pending.get(p);
        if (pending) runs.set(pending.run, pending.first);
        else if (!missing.includes(p)) missing.push(p);
      }
    }
    // adjacent missing pages are one request
    missing.sort((x, y) => x - y);
    for (let i = 0; i < missing.length; ) {
      let j = i + 1;
      while (j < missing.length && missing[j] === missing[j - 1] + 1) j++;
      runs.set(this.startRun(missing[i], j - i), missing[i]);
      i = j;
    }
    const arrived = await Promise.all(
      [...runs].map(([run, first]) => run.join(signal).then((b) => [first, b] as const)),
    );
    if (this.closed) throw this.closed;
    for (const [first, bytes] of arrived)
      for (let i = 0; i * ps < bytes.length; i++)
        got.set(first + i, bytes.subarray(i * ps, Math.min(bytes.length, (i + 1) * ps)));
    return spans.map(([a, b]) => {
      if (a === b) return new Uint8Array(0);
      const p = Math.floor(a / ps);
      if (Math.floor((b - 1) / ps) === p) return got.get(p)!.subarray(a - p * ps, b - p * ps);
      const out = new Uint8Array(b - a);
      for (let q = p; q * ps < b; q++) {
        const page = got.get(q)!;
        const from = Math.max(a, q * ps);
        const to = Math.min(b, (q + 1) * ps);
        out.set(page.subarray(from - q * ps, to - q * ps), from - a);
      }
      return out;
    });
  }

  /** Closes the store: pages let go, requests in flight stopped, later reads rejected with `why`. */
  close(why: Error) {
    if (this.closed) return;
    this.closed = why;
    this.pages.clear();
    this.whole = undefined;
    for (const { run } of this.pending.values()) if (!run.settled) run.controller.abort(why);
    this.pending.clear();
  }

  release() {
    this.holders--;
    if (this.holders > 0) return;
    this.close(new AssetRangeError('the reader was disposed', this.name));
    this.onEmpty();
  }
}

/** A reader holding one reference of a store. */
function storeReader(store: PageStore): AssetRangeReader {
  let disposed = false;
  return {
    size: store.size,
    read(ranges, signal) {
      if (disposed) return Promise.reject(new AssetRangeError('the reader was disposed', store.name));
      return store.read(ranges, signal);
    },
    dispose() {
      if (disposed) return;
      disposed = true;
      store.release();
    },
  };
}

/** The bytes of a response body, at most `limit` (`over` is thrown past it, also on a declared length over it). */
async function readBody(response: Response, limit: number, over: () => Error): Promise<Uint8Array> {
  const declared = Number(response.headers.get('Content-Length') ?? Number.NaN);
  const encoded = (response.headers.get('Content-Encoding') ?? 'identity') !== 'identity';
  if (!encoded && Number.isSafeInteger(declared) && declared > limit) {
    await response.body?.cancel().catch(() => {});
    throw over();
  }
  const chunks: Uint8Array[] = [];
  let loaded = 0;
  const reader = response.body?.getReader();
  while (reader) {
    const next = await reader.read();
    if (next.done) break;
    if (loaded + next.value.length > limit) {
      await reader.cancel().catch(() => {});
      throw over();
    }
    chunks.push(next.value);
    loaded += next.value.length;
  }
  if (chunks.length === 1) return chunks[0];
  const out = new Uint8Array(loaded);
  for (let i = 0, o = 0; i < chunks.length; o += chunks[i].length, i++) out.set(chunks[i], o);
  return out;
}

/** `bytes a-b/total` of a `Content-Range`, or undefined when it is not one with a total. */
function contentRange(header: string | null): { a: number; b: number; total: number } | undefined {
  const m = /^\s*bytes\s+(\d+)-(\d+)\/(\d+)\s*$/i.exec(header ?? '');
  if (!m) return undefined;
  const [a, b, total] = [Number(m[1]), Number(m[2]), Number(m[3])];
  return Number.isSafeInteger(total) && a <= b && b < total ? { a, b, total } : undefined;
}

type Net = Pick<RangeReadOptions, 'credentials' | 'redirect' | 'maxBytes'>;

/**
 * Opens a URL's store over HTTP: the first page is the probe, which says the size and the validator (a strong `ETag`,
 * else `Last-Modified`), or that the server ignores `Range` (a 200: the body read once, under `maxBytes`, and held
 * whole), or that it hides `Content-Range` from a cross-origin reader (the asset then read whole the same way). Later
 * pages must be what was asked of the same asset: their total and validator are compared here. No `If-Range`: RFC 9110
 * forbids it with a weak validator, and it is not CORS-safelisted (a preflight a page); `Range` alone is (§10).
 */
function httpStore(asset: Asset.Url, net: Net, pageSize: number, maxPages: number, onEmpty: () => void) {
  const url = asset.url;
  let validator: { name: 'ETag' | 'Last-Modified'; value: string } | undefined;
  const request = async (a: number, b: number | undefined, signal: AbortSignal) => {
    const headers: Record<string, string> = { ...asset.headers };
    if (b !== undefined) headers.Range = `bytes=${a}-${b - 1}`;
    return fetch(url, {
      signal,
      method: 'GET',
      headers,
      ...(net.credentials !== undefined ? { credentials: net.credentials } : {}),
      ...(net.redirect !== undefined ? { redirect: net.redirect } : {}),
    });
  };
  const fail = async (response: Response, why: string) => {
    await response.body?.cancel().catch(() => {});
    return new AssetRangeError(why, url);
  };
  const status = async (response: Response) => {
    await response.body?.cancel().catch(() => {});
    return new Error(`Download failed with status code ${response.status}`);
  };
  const store = new PageStore(
    url,
    pageSize,
    maxPages,
    async (a, b, signal) => {
      const response = await request(a, b, signal);
      if (response.status === 200)
        throw await fail(response, 'the asset changed, or its server stopped serving ranges');
      if (response.status === 416) throw await fail(response, `a range in the asset refused (416): ${a}-${b - 1}`);
      if (response.status !== 206) throw await status(response);
      const range = contentRange(response.headers.get('Content-Range'));
      if (!range || range.a !== a || range.b !== b - 1)
        throw await fail(response, `Content-Range ${response.headers.get('Content-Range')} for bytes ${a}-${b - 1}`);
      if (range.total !== store.size)
        throw await fail(response, `the size changed from ${store.size} to ${range.total}`);
      if (validator && response.headers.get(validator.name) !== validator.value)
        throw await fail(response, `the ${validator.name} changed between pages`);
      return readBody(response, b - a, () => new AssetRangeError('a body longer than its range', url));
    },
    onEmpty,
  );
  /** The probe: page 0, or the whole asset from a server that ignores `Range`. */
  const probe = async (signal: AbortSignal) => {
    let response = await request(0, pageSize, signal);
    if (response.status === 206 && response.headers.get('Content-Range') === null) {
      // a cross-origin server that does not expose `Content-Range`: its ranges cannot be checked, so the asset is read
      // whole, once, under `maxBytes`
      await response.body?.cancel().catch(() => {});
      response = await request(0, undefined, signal);
      if (response.status !== 200) throw await status(response);
    }
    if (response.status === 200) {
      const limit = net.maxBytes ?? Number.POSITIVE_INFINITY;
      const bytes = await readBody(response, limit, () => new SizeLimitExceeded(limit, url));
      store.size = bytes.length;
      store.whole = bytes;
      return;
    }
    if (response.status === 416) {
      // an empty asset has no byte 0: the only refusal of the first page that is not an error
      if (/^\s*bytes\s+\*\/0\s*$/i.test(response.headers.get('Content-Range') ?? '')) {
        await response.body?.cancel().catch(() => {});
        store.whole = new Uint8Array(0);
        return;
      }
      throw await fail(response, 'the first page refused (416)');
    }
    if (response.status !== 206) throw await status(response);
    const range = contentRange(response.headers.get('Content-Range'));
    if (!range || range.a !== 0 || range.b !== Math.min(pageSize, range.total) - 1)
      throw await fail(response, `Content-Range ${response.headers.get('Content-Range')} for bytes 0-${pageSize - 1}`);
    if (net.maxBytes !== undefined && range.total > net.maxBytes) {
      await response.body?.cancel().catch(() => {});
      throw new SizeLimitExceeded(net.maxBytes, url);
    }
    store.size = range.total;
    // a weak ETag (`W/"…"`) says the bytes may differ: no validator of bytes, so `Last-Modified` is compared instead
    const tag = response.headers.get('ETag');
    const etag = tag && !/^W\//.test(tag) ? tag : null;
    const modified = response.headers.get('Last-Modified');
    validator = etag
      ? { name: 'ETag', value: etag }
      : modified
        ? { name: 'Last-Modified', value: modified }
        : undefined;
    const bytes = await readBody(response, range.b + 1, () => new AssetRangeError('a body longer than its range', url));
    if (bytes.length !== range.b + 1)
      throw new AssetRangeError(`${bytes.length} bytes for a range of ${range.b + 1}`, url);
    store.keep(0, bytes);
  };
  return { store, probe };
}

/** A `file://` URL in Node.js: `stat` for the size (checked against `maxBytes`), then reads at an offset. */
function nodeFileStore(asset: Asset.Url, net: Net, pageSize: number, maxPages: number, onEmpty: () => void) {
  const filename = asset.url.substring('file://'.length);
  const fs = getFS();
  const call = <T>(f: (cb: (err: Error | null, v: T) => void) => void) =>
    new Promise<T>((resolve, reject) => f((err, v) => (err ? reject(err) : resolve(v))));
  const store = new PageStore(
    filename,
    pageSize,
    maxPages,
    async (a, b, signal) => {
      const fd = await call<number>((cb) => fs.open(filename, 'r', cb));
      try {
        const out = new Uint8Array(b - a);
        let at = 0;
        while (at < out.length) {
          const n = await call<number>((cb) => fs.read(fd, out, at, out.length - at, a + at, cb));
          if (n === 0) break;
          at += n;
        }
        if (signal.aborted) throw abortedBy(signal);
        return at === out.length ? out : out.subarray(0, at);
      } finally {
        await call<void>((cb) => fs.close(fd, (err) => cb(err, undefined)));
      }
    },
    onEmpty,
  );
  const probe = async (signal: AbortSignal) => {
    const { size } = await call<{ size: number }>((cb) => fs.stat(filename, cb));
    if (signal.aborted) throw abortedBy(signal);
    if (net.maxBytes !== undefined && size > net.maxBytes) throw new SizeLimitExceeded(net.maxBytes, filename);
    store.size = size;
  };
  return { store, probe };
}

/** A file or a `Blob`: its bytes as the manager holds them, else `Blob.slice` per run (never a whole read). */
function fileStore(name: string, file: Blob, pageSize: number, maxPages: number, onEmpty: () => void) {
  const store = new PageStore(
    name,
    pageSize,
    maxPages,
    async (a, b, signal) => {
      const bytes = new Uint8Array(await file.slice(a, b).arrayBuffer());
      if (signal.aborted) throw abortedBy(signal);
      return bytes;
    },
    onEmpty,
  );
  store.size = file.size;
  return store;
}

/** What `defaultRanges` reads of its manager. */
interface RangeHost {
  readonly manager: AssetManager;
  /** The bytes a held file asset has in memory, or its file. */
  held(asset: Asset.File): { bytes?: Uint8Array; file?: Blob } | undefined;
  /** The URL cache's key of a read (`urlCacheKey`), undefined for one that is not cached (a body). */
  urlKey(asset: Asset.Url, net: Net): string | undefined;
}

/**
 * @internal `DefaultAssetManager`'s range reads (`openRange`): one page store a key, its readers counted; `clear` closes
 * them all (their reads reject). Not part of the public surface: `assets.ts` re-exports what is.
 */
export function defaultRanges(host: RangeHost) {
  const stores = new Map<string, { store: PageStore; ready: Shared<void> }>();

  function open(asset: Asset, options: RangeReadOptions = {}): Task<AssetRangeReader> {
    return Task.create(`Open ${nameOf(asset)}`, async () => {
      const { credentials, redirect, signal, maxBytes } = options;
      checkMaxBytes(maxBytes);
      const pageSize = options.pageSize ?? DefaultPageSize;
      const maxPages = options.maxPages ?? DefaultMaxPages;
      if (!(pageSize >= 4096 && pageSize <= 4 * 1024 * 1024 && (pageSize & (pageSize - 1)) === 0))
        throw new Error(`pageSize must be a power of two from 4 KiB to 4 MiB, not ${pageSize}`);
      if (!(Number.isSafeInteger(maxPages) && maxPages >= 1))
        throw new Error(`maxPages must be a whole number of pages, not ${maxPages}`);
      if (signal?.aborted) throw abortedBy(signal);
      const net: Net = {
        ...(credentials !== undefined ? { credentials } : {}),
        ...(redirect !== undefined ? { redirect } : {}),
        ...(maxBytes !== undefined ? { maxBytes } : {}),
      };

      let key: string;
      let make: (onEmpty: () => void) => { store: PageStore; probe?: (signal: AbortSignal) => Promise<void> };
      if (asset.kind === 'url') {
        const urlKey = host.urlKey(asset, net);
        // a request with a body is not a resource of bytes to range over: read whole, as `resolve` reads it
        if (urlKey === undefined) return wholeReader(host.manager, asset, options).run();
        key = `${urlKey}|${pageSize}|${maxPages}`;
        make =
          RUNNING_IN_NODEJS && asset.url.startsWith('file://')
            ? (onEmpty) => nodeFileStore(asset, net, pageSize, maxPages, onEmpty)
            : (onEmpty) => httpStore(asset, net, pageSize, maxPages, onEmpty);
      } else {
        const held = host.held(asset);
        if (held?.bytes) {
          const bytes = held.bytes;
          if (maxBytes !== undefined && bytes.length > maxBytes) throw new SizeLimitExceeded(maxBytes, asset.name);
          return bytesReader(bytes, asset.name, () => {});
        }
        const file = held?.file ?? asset.file;
        if (!(file instanceof Blob)) throw new Error(`Cannot resolve file asset '${asset.name}' (${asset.id})`);
        if (maxBytes !== undefined && file.size > maxBytes) throw new SizeLimitExceeded(maxBytes, asset.name);
        key = `file|${asset.id}|${pageSize}|${maxPages}`;
        make = (onEmpty) => ({ store: fileStore(asset.name, file, pageSize, maxPages, onEmpty) });
      }

      for (;;) {
        let entry = stores.get(key);
        if (!entry) {
          const made = make(() => {
            if (stores.get(key) === created) stores.delete(key);
          });
          const probe = made.probe;
          const created = { store: made.store, ready: new Shared(async (s) => (probe ? probe(s) : undefined)) };
          // a failed or abandoned open is not kept: the next open starts again
          created.ready.promise.catch(() => {
            if (stores.get(key) === created) stores.delete(key);
          });
          stores.set(key, created);
          entry = created;
        }
        const { store, ready } = entry;
        store.holders++;
        try {
          await ready.join(signal);
          if (store.closed) throw store.closed;
          return storeReader(store);
        } catch (e) {
          // the hold goes back: the last one closes the store (an open nobody waits for any more is stopped)
          store.release();
          // the open joined was abandoned by the openers before this one: not a failure of this one, open again
          if (!signal?.aborted && ready.controller.signal.aborted) continue;
          throw e;
        }
      }
    });
  }

  function clear() {
    for (const { store } of stores.values())
      store.close(new AssetRangeError('the asset manager was cleared', store.name));
    stores.clear();
  }

  return { open, clear };
}
