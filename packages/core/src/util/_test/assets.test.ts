/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { beforeEach, describe, expect, it, vi } from 'vitest';
import { Task } from '@molstar/core/task';
import { Asset, AssetManager, DefaultAssetManager } from '../assets.js';
import { Zip } from '../zip/zip.js';
import { gzipSync } from 'node:zlib';

// the byte cache is tested against a stubbed download: nothing here reaches the network
vi.mock('../data-source.js', async (importOriginal) => ({
  ...(await importOriginal<typeof import('../data-source.js')>()),
  fetchData: vi.fn(),
}));
const { fetchData } = await import('../data-source.js');
const fetchMock = vi.mocked(fetchData) as unknown as ReturnType<typeof vi.fn>;
// the inflate is the real one unless a test makes it fail once
vi.mock('../zip/zip.js', async (importOriginal) => {
  const actual = await importOriginal<typeof import('../zip/zip.js')>();
  return { ...actual, ungzip: vi.fn(actual.ungzip) };
});
const zip = await import('../zip/zip.js');
const ungzipMock = vi.mocked(zip.ungzip);

const enc = new TextEncoder();

/** Answers every download with the bytes of `text`, after `gate` if given. */
function serve(text: string, gate?: Promise<void>) {
  fetchMock.mockImplementation((params: any) =>
    Task.create(`Download ${params.url}`, async (ctx) => {
      await ctx.update({ message: 'Downloading...' });
      if (gate) await gate;
      await ctx.update({ message: 'Done' });
      return enc.encode(text);
    }),
  );
}

describe('DefaultAssetManager URL cache', () => {
  beforeEach(() => {
    fetchMock.mockReset();
  });

  it('two resolves of the same URL fetch once', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const a = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    const b = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    expect(`${a.data}`).toBe('abc');
    expect(`${b.data}`).toBe('abc');
    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(m.urlCacheSize).toBe(1);
  });

  it('the same Asset.Url resolved twice counts two references', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/a');
    const a = await m.resolve(asset, 'string').run();
    const b = await m.resolve(asset, 'string').run();
    expect(fetchMock).toHaveBeenCalledTimes(1);
    a.dispose();
    expect(m.has(asset)).toBe(true);
    b.dispose();
    expect(m.has(asset)).toBe(false);
    expect(m.urlCacheSize).toBe(0);
  });

  it('drops the entry after the last holder disposes and fetches again', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const a = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    const b = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    a.dispose();
    expect(m.urlCacheSize).toBe(1);
    const c = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    expect(fetchMock).toHaveBeenCalledTimes(1);
    b.dispose();
    c.dispose();
    expect(m.urlCacheSize).toBe(0);

    await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    expect(fetchMock).toHaveBeenCalledTimes(2);
  });

  it('shares one download between concurrent resolves', async () => {
    let open!: () => void;
    serve('shared', new Promise<void>((r) => (open = r)));
    const m = new DefaultAssetManager();
    const p = [1, 2, 3].map(() => m.resolve(Asset.Url('https://x.org/a'), 'string').run());
    await Promise.resolve();
    open();
    const r = await Promise.all(p);
    expect(r.map((w) => `${w.data}`)).toEqual(['shared', 'shared', 'shared']);
    expect(fetchMock).toHaveBeenCalledTimes(1);

    r[0].dispose();
    r[1].dispose();
    expect(m.urlCacheSize).toBe(1);
    r[2].dispose();
    expect(m.urlCacheSize).toBe(0);
  });

  it('does not cache a failure, and every waiter of it gets the error', async () => {
    let fail!: () => void;
    const gate = new Promise<void>((_, rej) => (fail = () => rej(new Error('boom'))));
    gate.catch(() => {});
    serve('x', gate);
    const m = new DefaultAssetManager();
    const p = [0, 1].map(() => m.resolve(Asset.Url('https://x.org/a'), 'string').run());
    const settled = Promise.allSettled(p);
    await Promise.resolve();
    fail();
    const r = await settled;
    expect(r.map((x) => x.status)).toEqual(['rejected', 'rejected']);
    expect(m.urlCacheSize).toBe(0);
    expect(m.assets.length).toBe(0);

    serve('ok');
    const w = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    expect(`${w.data}`).toBe('ok');
    expect(fetchMock).toHaveBeenCalledTimes(2);
  });

  it('an aborted downloader does not fail the waiters, one takes over', async () => {
    let open!: () => void;
    serve('first', new Promise<void>((r) => (open = r)));
    const m = new DefaultAssetManager();

    let abort!: (reason?: string) => void;
    const downloader = m.resolve(Asset.Url('https://x.org/a'), 'string').run((p) => (abort = p.requestAbort));
    const downloaderResult = downloader.then(
      () => 'resolved',
      (e) => (Task.isAbort(e) ? 'aborted' : 'failed'),
    );
    const waiter = m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    await new Promise((r) => setTimeout(r, 0));
    expect(fetchMock).toHaveBeenCalledTimes(1);

    serve('second'); // what the retry downloads
    abort('test');
    open();
    expect(await downloaderResult).toBe('aborted');
    const w = await waiter;
    expect(`${w.data}`).toBe('second');
    expect(fetchMock).toHaveBeenCalledTimes(2);
    w.dispose();
    expect(m.urlCacheSize).toBe(0);
  });

  it('serves different response types from the same bytes', async () => {
    serve('{"a":1}');
    const m = new DefaultAssetManager();
    const s = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    const b = await m.resolve(Asset.Url('https://x.org/a'), 'binary').run();
    const j = await m.resolve(Asset.Url('https://x.org/a'), 'json').run();
    expect(`${s.data}`).toBe('{"a":1}');
    expect(b.data).toBeInstanceOf(Uint8Array);
    expect(Array.from(b.data)).toEqual(Array.from(enc.encode('{"a":1}')));
    expect(j.data).toEqual({ a: 1 });
    expect(fetchMock).toHaveBeenCalledTimes(1);
  });

  it('keeps URLs and headers apart', async () => {
    serve('x');
    const m = new DefaultAssetManager();
    await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    await m.resolve(Asset.Url('https://x.org/b'), 'string').run();
    await m.resolve(Asset.Url('https://x.org/a', { headers: { Authorization: 'k' } }), 'string').run();
    await m.resolve(Asset.Url('https://x.org/a', { headers: { Authorization: 'k' } }), 'string').run();
    expect(fetchMock).toHaveBeenCalledTimes(3);
    expect(m.urlCacheSize).toBe(3);
  });

  it('does not cache requests with a body, nor store: false', async () => {
    serve('x');
    const m = new DefaultAssetManager();
    const a = await m.resolve(Asset.Url('https://x.org/q', { body: '{}' }), 'string').run();
    const b = await m.resolve(Asset.Url('https://x.org/q', { body: '{}' }), 'string').run();
    expect(fetchMock).toHaveBeenCalledTimes(2);
    expect(m.urlCacheSize).toBe(0);
    a.dispose();
    b.dispose();

    await m.resolve(Asset.Url('https://x.org/a'), 'string', false).run();
    await m.resolve(Asset.Url('https://x.org/a'), 'string', false).run();
    expect(fetchMock).toHaveBeenCalledTimes(4);
    expect(m.urlCacheSize).toBe(0);
  });

  it('gives holds back when an entry is deleted or the manager is cleared', async () => {
    serve('x');
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/a');
    await m.resolve(asset, 'string').run();
    await m.resolve(Asset.Url('https://x.org/b'), 'string').run();
    expect(m.urlCacheSize).toBe(2);
    m.delete(asset);
    expect(m.urlCacheSize).toBe(1);
    m.dispose();
    expect(m.urlCacheSize).toBe(0);
  });
});

describe('DefaultAssetManager failed reads (perf-audit-2 01/B12)', () => {
  const notJson = () => Asset.File(new File(['<html>error</html>'], 'error.json'));

  it('a file that fails to read is not held', async () => {
    const m = new DefaultAssetManager();
    const asset = notJson();
    await expect(m.resolve(asset, 'json').run()).rejects.toThrow();
    expect(m.has(asset)).toBe(false);
  });

  it('a failed read of a held file gives its hold back', async () => {
    const m = new DefaultAssetManager();
    const asset = notJson();
    const text = await m.resolve(asset, 'string').run();
    await expect(m.resolve(asset, 'json').run()).rejects.toThrow();
    text.dispose();
    expect(m.has(asset)).toBe(false);
  });

  it('a URL whose bytes fail to read is not held', async () => {
    serve('<html>error</html>');
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/a.json');
    await expect(m.resolve(asset, 'json').run()).rejects.toThrow();
    expect(m.has(asset)).toBe(false);
    const text = await m.resolve(asset, 'string').run();
    await expect(m.resolve(asset, 'json').run()).rejects.toThrow();
    text.dispose();
    expect(m.has(asset)).toBe(false);
    expect(m.urlCacheSize).toBe(0);
  });
});

describe('DefaultAssetManager held bytes (util/_spec/assets.md)', () => {
  beforeEach(() => {
    fetchMock.mockReset();
  });

  it('a shared binary read returns the held bytes, other reads a copy', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/a');
    const a = await m.resolve(asset, 'binary', { shared: true }).run();
    const b = await m.resolve(asset, 'binary', { shared: true }).run();
    const c = await m.resolve(asset, 'binary').run();
    expect(b.data).toBe(a.data);
    expect(c.data).not.toBe(a.data);
    expect(Array.from(c.data)).toEqual(Array.from(a.data));
    expect(fetchMock).toHaveBeenCalledTimes(1);
    for (const w of [a, b, c]) w.dispose();
    expect(m.has(asset)).toBe(false);
    expect(m.urlCacheSize).toBe(0);
  });

  it('options store: false bypasses the cache', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    await m.resolve(Asset.Url('https://x.org/a'), 'string', { store: false }).run();
    expect(m.urlCacheSize).toBe(0);
    expect(m.assets.length).toBe(0);
  });

  it('a file read without shared is not held in memory', async () => {
    const m = new DefaultAssetManager();
    const file = new File([enc.encode('abc')], 'x.txt');
    const spy = vi.spyOn(file, 'arrayBuffer');
    const asset = Asset.File(file);
    const a = await m.resolve(asset, 'string').run();
    const b = await m.resolve(asset, 'binary').run();
    expect(`${a.data}`).toBe('abc');
    expect(Array.from(b.data)).toEqual(Array.from(enc.encode('abc')));
    // read from the file each time, as before: nothing but the file is kept
    expect(spy.mock.calls.length).toBe(2);
  });

  it('a gzipped file is read and decompressed once, then read in every form', async () => {
    const m = new DefaultAssetManager();
    const file = new File([gzipSync(enc.encode('{"a":1}'))], 'x.json.gz');
    const spy = vi.spyOn(file, 'arrayBuffer');
    const asset = Asset.File(file);
    const a = await m.resolve(asset, 'binary', { shared: true }).run();
    const s = await m.resolve(asset, 'string').run();
    const j = await m.resolve(asset, 'json').run();
    const b = await m.resolve(asset, 'binary', { shared: true }).run();
    expect(`${s.data}`).toBe('{"a":1}');
    expect(j.data).toEqual({ a: 1 });
    expect(b.data).toBe(a.data);
    expect(spy.mock.calls.length).toBeLessThanOrEqual(1);
    for (const w of [a, s, j, b]) w.dispose();
    expect(m.has(asset)).toBe(false);
  });

  it('a zip file read as zip keeps its entries', async () => {
    const m = new DefaultAssetManager();
    const zipped = await Zip({ 'a.txt': enc.encode('a'), 'b.txt': enc.encode('b') }).run();
    const asset = Asset.File(new File([zipped], 'x.zip'));
    const z = await m.resolve(asset, 'zip').run();
    expect(Object.keys(z.data).sort()).toEqual(['a.txt', 'b.txt']);
  });

  it('a static file lets go of its bytes after its last holder and reads them again', async () => {
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/static');
    m.set(asset, new File([enc.encode('abc')], 'static.txt'), { isStatic: true });
    const a = await m.resolve(asset, 'string').run();
    a.dispose();
    expect(m.has(asset)).toBe(true);
    const b = await m.resolve(asset, 'string').run();
    expect(`${b.data}`).toBe('abc');
    expect(fetchMock).not.toHaveBeenCalled();
  });

  it('an entry of a download gives a file made of its bytes', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/a');
    await m.resolve(asset, 'string').run();
    const entry = m.get(asset)!;
    expect(new Uint8Array(await entry.file.arrayBuffer())).toEqual(enc.encode('abc'));
    expect(m.tryFindFilename('raw-data')).toBe(asset);
  });
});

describe('DefaultAssetManager decompress (util/_spec/assets.md §5)', () => {
  const gz = (text: string) => new Uint8Array(gzipSync(enc.encode(text)));

  /** Answers every download with `bytes`. */
  function serveBytes(bytes: Uint8Array) {
    fetchMock.mockImplementation((params: any) => Task.create(`Download ${params.url}`, async () => bytes.slice()));
  }

  beforeEach(() => {
    fetchMock.mockReset();
  });

  it('inflates a gzipped download by its bytes, once for every reader of the URL', async () => {
    serveBytes(gz('{"a":1}'));
    const m = new DefaultAssetManager();
    const raw = await m.resolve(Asset.Url('https://x.org/a'), 'binary', { shared: true }).run();
    expect(raw.data[0]).toBe(0x1f);
    const a = await m.resolve(Asset.Url('https://x.org/a'), 'binary', { shared: true, decompress: true }).run();
    const b = await m.resolve(Asset.Url('https://x.org/a'), 'binary', { shared: true, decompress: true }).run();
    const s = await m.resolve(Asset.Url('https://x.org/a'), 'string', { decompress: true }).run();
    const j = await m.resolve(Asset.Url('https://x.org/a'), 'json', { decompress: true }).run();
    expect(new TextDecoder().decode(a.data)).toBe('{"a":1}');
    expect(b.data).toBe(a.data);
    expect(`${s.data}`).toBe('{"a":1}');
    expect(j.data).toEqual({ a: 1 });
    expect(fetchMock).toHaveBeenCalledTimes(1);
  });

  it('passes bytes that are not gzip through', async () => {
    serveBytes(enc.encode('plain'));
    const m = new DefaultAssetManager();
    const s = await m.resolve(Asset.Url('https://x.org/a'), 'string', { decompress: true }).run();
    expect(`${s.data}`).toBe('plain');
  });

  it('inflates without storing too', async () => {
    serveBytes(gz('abc'));
    const m = new DefaultAssetManager();
    const s = await m.resolve(Asset.Url('https://x.org/a'), 'string', { store: false, decompress: true }).run();
    expect(`${s.data}`).toBe('abc');
    expect(m.urlCacheSize).toBe(0);
  });

  it('inflates a gzipped file without the .gz name, and a .gz file once', async () => {
    const m = new DefaultAssetManager();
    const unnamed = await m.resolve(Asset.File(new File([gz('abc')], 'x.cif')), 'string', { decompress: true }).run();
    expect(`${unnamed.data}`).toBe('abc');
    const named = await m.resolve(Asset.File(new File([gz('abc')], 'x.cif.gz')), 'string', { decompress: true }).run();
    expect(`${named.data}`).toBe('abc');
  });

  it('resolveData inflates for a manager that does not (a custom one)', async () => {
    const inner = new DefaultAssetManager();
    const asked: unknown[] = [];
    // a custom manager: no `decompresses`, and `decompress` never reaches the one it wraps
    const custom = {
      resolve: (asset: Asset, type: any, options: any) => {
        asked.push(options);
        return inner.resolve(asset, type, { ...options, decompress: false });
      },
    } as unknown as AssetManager;
    const file = Asset.File(new File([gz('{"a":1}')], 'x.json'));
    const j = await AssetManager.resolveData(custom, file, 'json', { decompress: true }).run();
    expect(j.data).toEqual({ a: 1 });
    expect(asked).toEqual([{ decompress: false }]);
    j.dispose();
    expect(inner.has(file)).toBe(false);
  });

  it('a reader takes over the inflate when the reader that started it was aborted', async () => {
    serveBytes(gz('abc'));
    const m = new DefaultAssetManager();
    await m.resolve(Asset.Url('https://x.org/a'), 'binary', { shared: true }).run();
    ungzipMock.mockClear();
    ungzipMock.mockImplementationOnce(() => Promise.reject(Task.Aborted('the first reader')));
    const s = await m.resolve(Asset.Url('https://x.org/a'), 'string', { decompress: true }).run();
    expect(`${s.data}`).toBe('abc');
    expect(ungzipMock).toHaveBeenCalledTimes(2);
  });

  it('refuses a gzip that declares more than its data can hold, and lets go of its hold', async () => {
    const bomb = gz('abc');
    new DataView(bomb.buffer).setUint32(bomb.length - 4, 2 ** 31 + 1, true);
    serveBytes(bomb);
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/bomb');
    await expect(m.resolve(asset, 'string', { decompress: true }).run()).rejects.toThrow(/compressed bytes can hold/);
    expect(m.has(asset)).toBe(false);
    expect(m.urlCacheSize).toBe(0);
  });
});

describe('DefaultAssetManager letGoOfBytes (util/_spec/assets.md §6)', () => {
  beforeEach(() => {
    fetchMock.mockReset();
  });

  it('a file lets go of its bytes and reads them again from the file', async () => {
    const m = new DefaultAssetManager();
    const file = new File([enc.encode('abc')], 'x.txt');
    const spy = vi.spyOn(file, 'arrayBuffer');
    const asset = Asset.File(file);
    await m.resolve(asset, 'binary', { shared: true }).run();
    expect(spy).toHaveBeenCalledTimes(1);
    m.letGoOfBytes(asset);
    const again = await m.resolve(asset, 'binary', { shared: true }).run();
    expect(new TextDecoder().decode(again.data)).toBe('abc');
    expect(spy).toHaveBeenCalledTimes(2);
  });

  it('a download keeps its bytes as a blob, leaves the URL cache, and reads them again without fetching', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/a');
    await m.resolve(asset, 'binary', { shared: true }).run();
    expect(m.urlCacheSize).toBe(1);
    m.letGoOfBytes(asset);
    expect(m.urlCacheSize).toBe(0);
    expect(new Uint8Array(await m.get(asset)!.file.arrayBuffer())).toEqual(enc.encode('abc'));
    const again = await m.resolve(asset, 'string').run();
    expect(`${again.data}`).toBe('abc');
    expect(fetchMock).toHaveBeenCalledTimes(1);
    m.release(asset);
    m.release(asset);
    expect(m.has(asset)).toBe(false);
  });

  it('a download whose bytes other assets of the URL share keeps them', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const a = Asset.Url('https://x.org/a');
    const b = Asset.Url('https://x.org/a');
    await m.resolve(a, 'binary', { shared: true }).run();
    const bytes = (await m.resolve(b, 'binary', { shared: true }).run()).data;
    m.letGoOfBytes(a);
    expect(m.urlCacheSize).toBe(1);
    expect((await m.resolve(b, 'binary', { shared: true }).run()).data).toBe(bytes);
  });
});
