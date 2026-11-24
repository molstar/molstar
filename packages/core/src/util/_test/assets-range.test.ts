/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { Task } from '@molstar/core/task';
import {
  Asset,
  type AssetManager,
  AssetRangeError,
  type AssetRangeReader,
  DefaultAssetManager,
  openRanges,
  type RangeReadOptions,
} from '../assets.js';
import * as DataSource from '../data-source.js';
import * as RangeModule from '../assets-range.js';
import { getFS } from '../internal/fs.js';

// Range reads (`_spec/assets.md` §10) against a fake `fetch` serving one buffer: nothing reaches the network. The
// server's toggles: an honest 206, one that ignores `Range` (200), a `Content-Range` other than asked, a total or an
// `ETag` that changes between pages, a 416.

const Page = 4096;
const opts: RangeReadOptions = { pageSize: Page };

/** Bytes `i & 255` at `i`, so a slice tells where it came from. */
function bytes(n: number) {
  return Uint8Array.from({ length: n }, (_, i) => (i * 7 + (i >> 8)) & 255);
}

interface Served {
  body: Uint8Array;
  /** Answers a range with the whole body (200), as a server that ignores `Range`. */
  ignoresRange?: boolean;
  etag?: string;
  lastModified?: string;
  /** Rewrites the `Content-Range` sent. */
  contentRange?: (a: number, b: number, total: number) => string;
  /** The response waits for this. */
  gate?: Promise<void>;
  /** Answers every range with 416. */
  refuse?: boolean;
  /** Sends no `Content-Range`, as a cross-origin server that does not expose it. */
  hideContentRange?: boolean;
}

interface Request {
  url: string;
  range?: string;
  ifRange?: string;
  init: RequestInit;
}

function serve(routes: Record<string, Served>) {
  const log: Request[] = [];
  const fake = vi.fn(async (input: RequestInfo | URL, init: RequestInit = {}) => {
    const url = String(input);
    const headers = (init.headers ?? {}) as Record<string, string>;
    log.push({ url, range: headers.Range, ifRange: headers['If-Range'], init });
    const signal = init.signal ?? undefined;
    const s = routes[url];
    if (s?.gate) {
      await new Promise<void>((resolve, reject) => {
        signal?.addEventListener('abort', () => reject(signal.reason), { once: true });
        s.gate!.then(resolve);
      });
    }
    if (signal?.aborted) throw signal.reason;
    if (!s) return new Response('', { status: 404 });
    const size = s.body.length;
    const meta: Record<string, string> = {};
    if (s.etag) meta.ETag = s.etag;
    if (s.lastModified) meta['Last-Modified'] = s.lastModified;
    const m = /^bytes=(\d+)-(\d+)$/.exec(headers.Range ?? '');
    if (!m || s.ignoresRange)
      return new Response(s.body.slice(), { status: 200, headers: { ...meta, 'Content-Length': String(size) } });
    const a = Number(m[1]);
    const b = Math.min(Number(m[2]), size - 1);
    if (s.refuse || a >= size)
      return new Response('', { status: 416, headers: { 'Content-Range': `bytes */${size}` } });
    const contentRange = s.contentRange?.(a, b, size) ?? `bytes ${a}-${b}/${size}`;
    return new Response(s.body.slice(a, b + 1), {
      status: 206,
      headers: {
        ...meta,
        ...(s.hideContentRange ? {} : { 'Content-Range': contentRange }),
        'Content-Length': String(b + 1 - a),
      },
    });
  });
  vi.stubGlobal('fetch', fake);
  return log;
}

const ranged = (log: Request[]) => log.filter((r) => r.range !== undefined).length;

async function open(m: AssetManager, asset: Asset, o: RangeReadOptions = opts) {
  return openRanges(m, asset, o).run();
}

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('range reads of a URL', () => {
  const url = 'https://x.org/a.bw';
  const data = bytes(5 * Page + 100);

  it('opening reads the first page and learns the size', async () => {
    const log = serve({ [url]: { body: data } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url));
    expect(r.size).toBe(data.length);
    expect(log.map((l) => l.range)).toEqual([`bytes=0-${Page - 1}`]);
    r.dispose();
  });

  it('ranges across pages, past the end clipped, empty ones empty', async () => {
    serve({ [url]: { body: data } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url));
    const [across, tail, past, empty] = await r.read([
      { offset: Page - 10, length: 20 },
      { offset: data.length - 50, length: 200 },
      { offset: data.length + 10, length: 5 },
      { offset: 7, length: 0 },
    ]);
    expect(across).toEqual(data.slice(Page - 10, Page + 10));
    expect(tail).toEqual(data.slice(data.length - 50));
    expect(past.length).toBe(0);
    expect(empty.length).toBe(0);
    r.dispose();
  });

  it('a range inside one page is a view of the cached page, not a copy', async () => {
    serve({ [url]: { body: data } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url));
    const [a] = await r.read([{ offset: 10, length: 5 }]);
    const [b] = await r.read([{ offset: 20, length: 5 }]);
    expect(a.buffer).toBe(b.buffer);
    r.dispose();
  });

  it('many ranges in one call: adjacent missing pages are one request, others go in parallel', async () => {
    const log = serve({ [url]: { body: data } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url));
    log.length = 0;
    // pages 1 and 2 adjacent (one request), page 4 apart (another), page 0 held already (none)
    const out = await r.read([
      { offset: Page + 5, length: 10 },
      { offset: 2 * Page + 5, length: 10 },
      { offset: 4 * Page + 1, length: 3 },
      { offset: 3, length: 3 },
    ]);
    expect(new Set(log.map((l) => l.range))).toEqual(
      new Set([`bytes=${Page}-${3 * Page - 1}`, `bytes=${4 * Page}-${5 * Page - 1}`]),
    );
    expect(log.length).toBe(2);
    expect(out[2]).toEqual(data.slice(4 * Page + 1, 4 * Page + 4));
    r.dispose();
  });

  it('a repeated read is served from the cache: no request', async () => {
    const log = serve({ [url]: { body: data } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url));
    await r.read([{ offset: 3 * Page, length: 100 }]);
    const n = log.length;
    const [again] = await r.read([{ offset: 3 * Page + 10, length: 50 }]);
    expect(log.length).toBe(n);
    expect(again).toEqual(data.slice(3 * Page + 10, 3 * Page + 60));
    r.dispose();
  });

  it('concurrent reads share a page in flight, and two readers of a URL share its pages', async () => {
    let release!: () => void;
    const log = serve({ [url]: { body: data } });
    const m = new DefaultAssetManager();
    const a = await open(m, Asset.Url(url));
    const b = await open(m, Asset.Url(url));
    expect(ranged(log)).toBe(1);
    // the server waits from now on
    serve({ [url]: { body: data, gate: new Promise<void>((resolve) => (release = resolve)) } });
    const reads = [a.read([{ offset: 2 * Page, length: 4 }]), b.read([{ offset: 2 * Page + 8, length: 4 }])];
    await new Promise((resolve) => setTimeout(resolve, 5));
    release();
    const [[x], [y]] = await Promise.all(reads);
    expect(x).toEqual(data.slice(2 * Page, 2 * Page + 4));
    expect(y).toEqual(data.slice(2 * Page + 8, 2 * Page + 12));
    expect(vi.mocked(fetch)).toHaveBeenCalledTimes(1);
    a.dispose();
    b.dispose();
  });

  it('an abort mid-read rejects as aborted; the page is not kept and the next read fetches it again', async () => {
    const m = new DefaultAssetManager();
    serve({ [url]: { body: data } });
    const r = await open(m, Asset.Url(url));
    const stop = new AbortController();
    const log = serve({ [url]: { body: data, gate: new Promise<void>(() => {}) } });
    const read = r.read([{ offset: 3 * Page, length: 4 }], stop.signal);
    await new Promise((resolve) => setTimeout(resolve, 5));
    stop.abort('panned away');
    await expect(read).rejects.toSatisfy((e) => Task.isAbort(e) && (e as Task.Aborted).reason === 'panned away');
    // the request was stopped: its only waiter left
    expect((log[0].init.signal as AbortSignal).aborted).toBe(true);
    const again = serve({ [url]: { body: data } });
    const [x] = await r.read([{ offset: 3 * Page, length: 4 }]);
    expect(x).toEqual(data.slice(3 * Page, 3 * Page + 4));
    expect(ranged(again)).toBe(1);
    r.dispose();
  });

  it('a waiter that aborts leaves a shared page to the others', async () => {
    const m = new DefaultAssetManager();
    serve({ [url]: { body: data } });
    const a = await open(m, Asset.Url(url));
    const b = await open(m, Asset.Url(url));
    let release!: () => void;
    serve({ [url]: { body: data, gate: new Promise<void>((resolve) => (release = resolve)) } });
    const stop = new AbortController();
    const left = a.read([{ offset: 3 * Page, length: 4 }], stop.signal);
    const stays = b.read([{ offset: 3 * Page, length: 4 }]);
    stop.abort();
    await expect(left).rejects.toSatisfy(Task.isAbort);
    release();
    expect((await stays)[0]).toEqual(data.slice(3 * Page, 3 * Page + 4));
    a.dispose();
    b.dispose();
  });

  it('a read that joins a request its only waiter just aborted fetches again, not failing', async () => {
    const m = new DefaultAssetManager();
    serve({ [url]: { body: data } });
    const a = await open(m, Asset.Url(url));
    const b = await open(m, Asset.Url(url));
    let release!: () => void;
    serve({ [url]: { body: data, gate: new Promise<void>((resolve) => (release = resolve)) } });
    const stop = new AbortController();
    const left = a.read([{ offset: 3 * Page, length: 4 }], stop.signal);
    await new Promise((resolve) => setTimeout(resolve, 5));
    stop.abort();
    // in the same turn, before the stopped request is forgotten
    const joins = b.read([{ offset: 3 * Page, length: 4 }]);
    release();
    await expect(left).rejects.toSatisfy(Task.isAbort);
    expect((await joins)[0]).toEqual(data.slice(3 * Page, 3 * Page + 4));
    a.dispose();
    b.dispose();
  });

  it('dispose frees the pages with the last holder; a disposed reader rejects', async () => {
    const log = serve({ [url]: { body: data } });
    const m = new DefaultAssetManager();
    const a = await open(m, Asset.Url(url));
    const b = await open(m, Asset.Url(url));
    await a.read([{ offset: 2 * Page, length: 4 }]);
    a.dispose();
    const n = log.length;
    await b.read([{ offset: 2 * Page, length: 4 }]);
    expect(log.length).toBe(n);
    b.dispose();
    await expect(a.read([{ offset: 0, length: 1 }])).rejects.toBeInstanceOf(AssetRangeError);
    const c = await open(m, Asset.Url(url));
    expect(log.length).toBe(n + 1);
    c.dispose();
  });

  it('the manager’s clear closes its readers: their reads reject', async () => {
    serve({ [url]: { body: data } });
    const m = new DefaultAssetManager();
    const r = await open(m, Asset.Url(url));
    m.clear();
    await expect(r.read([{ offset: 0, length: 1 }])).rejects.toMatchObject({ name: 'AssetRangeError' });
  });

  it('a cached page holds its own bytes, not its whole request’s: maxPages bounds memory', async () => {
    serve({ [url]: { body: data } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url), { pageSize: Page, maxPages: 1 });
    // pages 1 to 4 in one request; then page 4 alone from the cache
    await r.read([{ offset: Page, length: 4 * Page }]);
    const [x] = await r.read([{ offset: 4 * Page + 1, length: 2 }]);
    expect(x.buffer.byteLength).toBe(Page);
    r.dispose();
  });

  it('least recently used pages go first past maxPages', async () => {
    const log = serve({ [url]: { body: data } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url), { pageSize: Page, maxPages: 2 });
    await r.read([{ offset: Page, length: 1 }]);
    await r.read([{ offset: 2 * Page, length: 1 }]);
    const n = log.length;
    // page 0 (opened) is the oldest: gone; page 1 kept
    await r.read([{ offset: Page, length: 1 }]);
    expect(log.length).toBe(n);
    await r.read([{ offset: 0, length: 1 }]);
    expect(log.length).toBe(n + 1);
    r.dispose();
  });

  it('credentials, redirect and the asset’s headers reach the requests; no If-Range (not CORS-safelisted)', async () => {
    const log = serve({ [url]: { body: data, etag: '"v1"' } });
    const asset = Asset.Url(url, { headers: { Accept: 'application/octet-stream' } });
    const r = await open(new DefaultAssetManager(), asset, { ...opts, credentials: 'include', redirect: 'error' });
    await r.read([{ offset: 3 * Page, length: 1 }]);
    for (const l of log) {
      expect(l.init.credentials).toBe('include');
      expect(l.init.redirect).toBe('error');
      expect((l.init.headers as Record<string, string>).Accept).toBe('application/octet-stream');
    }
    expect(log.length).toBe(2);
    expect(log.every((l) => l.ifRange === undefined)).toBe(true);
    r.dispose();
  });

  it('reads with other credentials do not share pages', async () => {
    const log = serve({ [url]: { body: data } });
    const m = new DefaultAssetManager();
    const a = await open(m, Asset.Url(url), { ...opts, credentials: 'omit' });
    const b = await open(m, Asset.Url(url), { ...opts, credentials: 'include' });
    expect(ranged(log)).toBe(2);
    a.dispose();
    b.dispose();
  });

  it('a server that ignores Range: the body read once, under maxBytes, every range served from it', async () => {
    const log = serve({ [url]: { body: data, ignoresRange: true } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url));
    expect(r.size).toBe(data.length);
    const [x, y] = await r.read([
      { offset: 4 * Page, length: 10 },
      { offset: 1, length: 3 },
    ]);
    expect(x).toEqual(data.slice(4 * Page, 4 * Page + 10));
    expect(y).toEqual(data.slice(1, 4));
    expect(log.length).toBe(1);
    r.dispose();
    serve({ [url]: { body: data, ignoresRange: true } });
    await expect(
      open(new DefaultAssetManager(), Asset.Url(url), { ...opts, maxBytes: data.length - 1 }),
    ).rejects.toMatchObject({ name: 'SizeLimitExceeded' });
  });

  it('maxBytes bounds the asset’s size from Content-Range; a bound that is not a size refuses', async () => {
    serve({ [url]: { body: data } });
    await expect(
      open(new DefaultAssetManager(), Asset.Url(url), { ...opts, maxBytes: data.length - 1 }),
    ).rejects.toMatchObject({ name: 'SizeLimitExceeded' });
    const r = await open(new DefaultAssetManager(), Asset.Url(url), { ...opts, maxBytes: data.length });
    expect(r.size).toBe(data.length);
    r.dispose();
    await expect(open(new DefaultAssetManager(), Asset.Url(url), { ...opts, maxBytes: -1 })).rejects.toThrow();
    await expect(open(new DefaultAssetManager(), Asset.Url(url), { pageSize: 1000 })).rejects.toThrow(/pageSize/);
  });

  it.each([
    ['the start', (a: number, b: number, t: number) => `bytes ${a + 1}-${b}/${t}`],
    ['the end', (a: number, b: number, t: number) => `bytes ${a}-${b - 1}/${t}`],
    ['the total missing', (a: number, b: number) => `bytes ${a}-${b}/*`],
  ])('a Content-Range other than asked (%s) fails the read', async (_, contentRange) => {
    serve({ [url]: { body: data, contentRange } });
    await expect(open(new DefaultAssetManager(), Asset.Url(url))).rejects.toBeInstanceOf(AssetRangeError);
    // on a later page too
    const m = new DefaultAssetManager();
    serve({ [url]: { body: data } });
    const r = await open(m, Asset.Url(url));
    serve({ [url]: { body: data, contentRange } });
    await expect(r.read([{ offset: 3 * Page, length: 1 }])).rejects.toBeInstanceOf(AssetRangeError);
    r.dispose();
  });

  it('a total that changes between pages fails the read', async () => {
    const m = new DefaultAssetManager();
    serve({ [url]: { body: data } });
    const r = await open(m, Asset.Url(url));
    serve({ [url]: { body: bytes(data.length + 1) } });
    await expect(r.read([{ offset: 3 * Page, length: 1 }])).rejects.toThrow(/size changed/);
    r.dispose();
  });

  it.each([
    ['ETag', { etag: '"v1"' }, { etag: '"v2"' }],
    [
      'Last-Modified',
      { lastModified: 'Mon, 02 Nov 2026 10:00:00 GMT' },
      { lastModified: 'Fri, 20 Nov 2026 10:00:00 GMT' },
    ],
  ])('a changed %s between pages fails with AssetRangeError', async (name, before, after) => {
    const m = new DefaultAssetManager();
    serve({ [url]: { body: data, ...before } });
    const r = await open(m, Asset.Url(url));
    serve({ [url]: { body: data, ...after } });
    await expect(r.read([{ offset: 3 * Page, length: 1 }])).rejects.toThrow(new RegExp(`${name} changed`));
    // a 200 after the first page: the asset changed (or its server stopped serving ranges), the read fails too
    serve({ [url]: { body: data, ...after, ignoresRange: true } });
    await expect(r.read([{ offset: 4 * Page, length: 1 }])).rejects.toBeInstanceOf(AssetRangeError);
    r.dispose();
  });

  it('a weak ETag is no validator: Last-Modified is compared instead', async () => {
    const m = new DefaultAssetManager();
    const lastModified = 'Mon, 02 Nov 2026 10:00:00 GMT';
    serve({ [url]: { body: data, etag: 'W/"1"', lastModified } });
    const r = await open(m, Asset.Url(url));
    // a weak tag may change with the bytes the same
    serve({ [url]: { body: data, etag: 'W/"2"', lastModified } });
    expect((await r.read([{ offset: 3 * Page, length: 2 }]))[0]).toEqual(data.slice(3 * Page, 3 * Page + 2));
    serve({ [url]: { body: data, etag: 'W/"2"', lastModified: 'Fri, 20 Nov 2026 10:00:00 GMT' } });
    await expect(r.read([{ offset: 4 * Page, length: 1 }])).rejects.toThrow(/Last-Modified changed/);
    r.dispose();
  });

  it('a server that hides Content-Range (CORS): the asset read whole once, under maxBytes', async () => {
    const log = serve({ [url]: { body: data, hideContentRange: true } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url));
    expect(r.size).toBe(data.length);
    expect(log.map((l) => l.range)).toEqual([`bytes=0-${Page - 1}`, undefined]);
    expect((await r.read([{ offset: 4 * Page, length: 3 }]))[0]).toEqual(data.slice(4 * Page, 4 * Page + 3));
    expect(log.length).toBe(2);
    r.dispose();
    serve({ [url]: { body: data, hideContentRange: true } });
    await expect(
      open(new DefaultAssetManager(), Asset.Url(url), { ...opts, maxBytes: data.length - 1 }),
    ).rejects.toMatchObject({ name: 'SizeLimitExceeded' });
  });

  it('a 416 on a valid range fails; an empty asset opens empty', async () => {
    serve({ [url]: { body: data, refuse: true } });
    await expect(open(new DefaultAssetManager(), Asset.Url(url))).rejects.toBeInstanceOf(AssetRangeError);
    serve({ [url]: { body: new Uint8Array(0) } });
    const r = await open(new DefaultAssetManager(), Asset.Url(url));
    expect(r.size).toBe(0);
    expect(await r.read([{ offset: 0, length: 10 }])).toEqual([new Uint8Array(0)]);
    r.dispose();
  });

  it('a failed open is not kept: the next open fetches again', async () => {
    serve({ [url]: { body: data, refuse: true } });
    const m = new DefaultAssetManager();
    await expect(open(m, Asset.Url(url))).rejects.toBeInstanceOf(AssetRangeError);
    serve({ [url]: { body: data } });
    const r = await open(m, Asset.Url(url));
    expect(r.size).toBe(data.length);
    r.dispose();
  });

  it('a status that is not a range or the body is a failed download', async () => {
    serve({});
    await expect(open(new DefaultAssetManager(), Asset.Url(url))).rejects.toThrow(/status code 404/);
  });

  it('an aborted open rejects as aborted', async () => {
    serve({ [url]: { body: data, gate: new Promise<void>(() => {}) } });
    const stop = new AbortController();
    const opening = open(new DefaultAssetManager(), Asset.Url(url), { ...opts, signal: stop.signal });
    await new Promise((resolve) => setTimeout(resolve, 5));
    stop.abort();
    await expect(opening).rejects.toSatisfy(Task.isAbort);
  });
});

describe('range reads of a file', () => {
  const data = bytes(3 * Page + 9);

  it('a File by Blob.slice, never read whole', async () => {
    const file = new File([data], 'x.bw');
    const slice = vi.spyOn(file, 'slice');
    const whole = vi.spyOn(file, 'arrayBuffer');
    const r = await open(new DefaultAssetManager(), Asset.File(file));
    expect(r.size).toBe(data.length);
    const [x] = await r.read([{ offset: Page - 2, length: Page + 4 }]);
    expect(x).toEqual(data.slice(Page - 2, 2 * Page + 2));
    expect(whole).not.toHaveBeenCalled();
    // pages 0 to 2, one slice
    expect(slice.mock.calls).toEqual([[0, 3 * Page]]);
    await expect(
      open(new DefaultAssetManager(), Asset.File(file), { ...opts, maxBytes: data.length - 1 }),
    ).rejects.toMatchObject({ name: 'SizeLimitExceeded' });
    r.dispose();
  });

  it('a file whose bytes the manager holds is served from them', async () => {
    const m = new DefaultAssetManager();
    const file = new File([data], 'x.bw');
    const asset = Asset.File(file);
    const held = await m.resolve(asset, 'binary', { shared: true }).run();
    const slice = vi.spyOn(file, 'slice');
    const r = await open(m, asset);
    const [x] = await r.read([{ offset: 10, length: 5 }]);
    expect(x.buffer).toBe(held.data.buffer);
    expect(slice).not.toHaveBeenCalled();
    r.dispose();
    held.dispose();
  });

  it('file:// in Node.js: stat against maxBytes, then reads at an offset', async () => {
    let prev: ReturnType<typeof getFS> | undefined;
    try {
      prev = getFS();
    } catch {
      prev = undefined;
    }
    const reads: [number, number][] = [];
    const fs = {
      stat: (_: string, cb: (e: Error | null, s: { size: number }) => void) => cb(null, { size: data.length }),
      open: (_: string, __: string, cb: (e: Error | null, fd: number) => void) => cb(null, 3),
      close: (_: number, cb: (e: Error | null) => void) => cb(null),
      read: (
        _: number,
        out: Uint8Array,
        at: number,
        length: number,
        position: number,
        cb: (e: Error | null, n: number) => void,
      ) => {
        reads.push([position, length]);
        const part = data.subarray(position, position + length);
        out.set(part, at);
        cb(null, part.length);
      },
      readFile: () => {
        throw new Error('read whole');
      },
    };
    try {
      DataSource.setFSModule(fs as unknown as Parameters<typeof DataSource.setFSModule>[0]);
      const m = new DefaultAssetManager();
      await expect(
        open(m, Asset.Url('file:///data/x.bw'), { ...opts, maxBytes: data.length - 1 }),
      ).rejects.toMatchObject({ name: 'SizeLimitExceeded' });
      expect(reads).toEqual([]);
      const r = await open(m, Asset.Url('file:///data/x.bw'), { ...opts, maxBytes: data.length });
      const [x] = await r.read([{ offset: 2 * Page + 3, length: 10 }]);
      expect(x).toEqual(data.slice(2 * Page + 3, 2 * Page + 13));
      expect(reads).toEqual([[2 * Page, Page]]);
      r.dispose();
    } finally {
      DataSource.setFSModule(prev as unknown as Parameters<typeof DataSource.setFSModule>[0]);
    }
  });
});

describe('openRanges', () => {
  it('a manager without openRange is served whole: the asset resolved once, ranges views of it', async () => {
    const data = bytes(1000);
    const resolve = vi.fn((asset: Asset) =>
      Task.create('read', async () => ({ data: data as Uint8Array<ArrayBuffer>, dispose: () => release(asset) })),
    );
    const release = vi.fn();
    const custom = { resolve } as unknown as AssetManager;
    const asset = Asset.Url('https://x.org/a');
    const r: AssetRangeReader = await openRanges(custom, asset, { maxBytes: 5000, credentials: 'omit' }).run();
    expect(resolve).toHaveBeenCalledWith(asset, 'binary', { shared: true, maxBytes: 5000, credentials: 'omit' });
    expect(r.size).toBe(1000);
    const [x] = await r.read([{ offset: 990, length: 50 }]);
    expect(x).toEqual(data.slice(990));
    expect(x.buffer).toBe(data.buffer);
    r.dispose();
    r.dispose();
    expect(release).toHaveBeenCalledTimes(1);
  });

  it('a request with a body is read whole by the default manager too', async () => {
    const log = serve({ 'https://x.org/q': { body: bytes(100) } });
    const r = await open(new DefaultAssetManager(), Asset.Url('https://x.org/q', { body: '{}' }));
    expect(r.size).toBe(100);
    expect(log[0].range).toBeUndefined();
    r.dispose();
  });
});

describe('the published surface', () => {
  it('the range module publishes its public names only; the internals are in a module the package does not export', () => {
    expect(Object.keys(RangeModule).sort()).toEqual(['AssetRangeError', 'openRanges']);
    expect('getFS' in DataSource).toBe(false);
    expect('defaultRanges' in RangeModule).toBe(false);
  });
});
