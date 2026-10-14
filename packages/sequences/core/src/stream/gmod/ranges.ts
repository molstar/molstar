/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Byte-range reads: what every reader here takes in place of a URL (README, "The transport"). A host reads HTTP
// ranges (`fetchRanges`), tests read recorded bytes. Binary formats are read with a `DataView` and the file's own byte
// order, never the host's.

/** `length` bytes of the file from `offset`; fewer where the file ends. Rejects with the failure (an abort too). */
export type RangeReader = (offset: number, length: number, signal?: AbortSignal) => Promise<Uint8Array>;

/** The reference transport: HTTP `Range` requests to `url` (a server that ignores `Range` is an error, not a download). */
export function fetchRanges(url: string): RangeReader {
  return async (offset, length, signal) => {
    const r = await fetch(url, { headers: { Range: `bytes=${offset}-${offset + length - 1}` }, signal });
    if (r.status === 416) return new Uint8Array(0); // starts at or past the end of the file
    if (r.status !== 206) throw new Error(`${url}: a range request answered ${r.status} ${r.statusText}, not 206`);
    return new Uint8Array(await r.arrayBuffer());
  };
}

/**
 * Work several callers wait on (a page, a record, a file's open): it runs under a controller of its own, never a
 * caller's signal, and is aborted only once every caller waiting on it has gone (CB-01: one request's abort failed
 * every request sharing its reads).
 */
export interface Shared<T> {
  readonly promise: Promise<T>;
  readonly controller: AbortController;
  waiters: number;
  settled: boolean;
}

/** `run` started under a controller of its own. */
export function share<T>(run: (signal: AbortSignal) => Promise<T>): Shared<T> {
  const controller = new AbortController();
  const promise = run(controller.signal);
  const s: Shared<T> = { promise, controller, waiters: 0, settled: false };
  const done = () => {
    s.settled = true;
  };
  promise.then(done, done);
  return s;
}

/** Whether a shared work can still be joined: not aborted for want of waiters. */
export const live = <T>(s: Shared<T> | null | undefined): s is Shared<T> => !!s && !s.controller.signal.aborted;

/**
 * `s`'s result, waited on under `signal`: an abort rejects this caller with its reason and leaves the work to the
 * others; the last to leave aborts it.
 */
export function join<T>(s: Shared<T>, signal?: AbortSignal): Promise<T> {
  if (signal?.aborted) return Promise.reject(signal.reason);
  s.waiters++;
  return new Promise<T>((resolve, reject) => {
    let left = false;
    const leave = () => {
      if (left) return false;
      left = true;
      signal?.removeEventListener('abort', aborted);
      s.waiters--;
      return true;
    };
    const aborted = () => {
      if (!leave()) return;
      if (s.waiters === 0 && !s.settled) s.controller.abort(signal?.reason);
      reject(signal?.reason);
    };
    signal?.addEventListener('abort', aborted, { once: true });
    s.promise.then(
      (v) => leave() && resolve(v),
      (e) => leave() && reject(e),
    );
  });
}

/**
 * A reader that fetches whole `page`-aligned pages and keeps them, so the many small reads of a binary search or a
 * tree walk are one request per page, once; a page is shared by its readers (`share`). A failed page is not kept; the
 * oldest go past `keep` pages (2 MB).
 */
export function paged(read: RangeReader, page = 8192, keep = 256): RangeReader {
  const pages = new Map<number, Shared<Uint8Array>>();
  const get = (k: number, signal?: AbortSignal) => {
    // A caller already gone starts no read (LA-37: one with no waiters, which nothing would abort).
    if (signal?.aborted) return Promise.reject(signal.reason);
    let s = pages.get(k);
    if (!live(s)) {
      const mine = share((own) => read(k * page, page, own));
      s = mine;
      pages.set(k, mine);
      mine.promise.catch(() => pages.get(k) === mine && pages.delete(k));
      if (pages.size > keep) pages.delete(pages.keys().next().value as number);
    }
    return join(s, signal);
  };
  return async (offset, length, signal) => {
    const [a, b] = [Math.floor(offset / page), Math.floor((offset + length - 1) / page)];
    const got = await Promise.all(Array.from({ length: b - a + 1 }, (_, i) => get(a + i, signal)));
    const held = got.reduce((n, p) => n + p.length, 0) - (offset - a * page);
    const out = new Uint8Array(Math.max(0, Math.min(length, held)));
    let at = 0;
    for (let i = 0; i < got.length && at < out.length; i++) {
      const from = i === 0 ? offset - a * page : 0;
      const bytes = got[i].subarray(from, from + out.length - at);
      out.set(bytes, at);
      at += bytes.length;
    }
    return out;
  };
}

/**
 * `fn`'s first result, kept and shared (`share`: a manifest request's abort does not fail the data requests waiting on
 * the same open); a failure is not kept, nor an open every caller left, so the next call tries again.
 */
export function once<T>(fn: (signal: AbortSignal) => Promise<T>): (signal?: AbortSignal) => Promise<T> {
  let kept: Shared<T> | null = null;
  return (signal) => {
    if (signal?.aborted) return Promise.reject(signal.reason);
    if (!live(kept)) {
      const mine = share(fn);
      kept = mine;
      mine.promise.catch(() => {
        if (kept === mine) kept = null;
      });
    }
    return join(kept, signal);
  };
}

/**
 * `bytes` decompressed: `deflate` is zlib (a BigWig or BigBed block), `gzip` a gzip member (a BGZF block). Stopped,
 * and an error, past `max` bytes (GI-07: a 61 KB BGZF block held 60 MB of zeros, where the format promises 64 KB).
 */
export async function inflate(bytes: Uint8Array, format: 'deflate' | 'gzip', max: number): Promise<Uint8Array> {
  // A view of an ArrayBuffer (DOM lib typings reject a SharedArrayBuffer-backed one as a BlobPart).
  const stream = new Blob([bytes as Uint8Array<ArrayBuffer>]).stream().pipeThrough(new DecompressionStream(format));
  const reader = stream.getReader();
  const parts: Uint8Array[] = [];
  let n = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    n += value.length;
    if (n > max) {
      await reader.cancel();
      throw new Error(`a compressed block inflates past ${max} bytes, its format's bound`);
    }
    parts.push(value);
  }
  if (parts.length === 1) return parts[0] as Uint8Array;
  const out = new Uint8Array(n);
  for (let i = 0, at = 0; i < parts.length; at += (parts[i] as Uint8Array).length, i++)
    out.set(parts[i] as Uint8Array, at);
  return out;
}

export const view = (b: Uint8Array) => new DataView(b.buffer, b.byteOffset, b.byteLength);
/** A 64-bit file offset as a number (files stay under 2^53 bytes). */
export const u64 = (dv: DataView, at: number, le: boolean) => Number(dv.getBigUint64(at, le));
const ascii = new TextDecoder('latin1');
export const text = (b: Uint8Array) => ascii.decode(b);
/** The text of `b` up to its first NUL. */
export const cstring = (b: Uint8Array) => {
  const end = b.indexOf(0);
  return text(end < 0 ? b : b.subarray(0, end));
};
