/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { beforeEach, describe, expect, it, vi } from 'vitest';
import { Task } from '@molstar/core/task';
import { Asset, DefaultAssetManager } from '../assets.js';

// the byte cache is tested against a stubbed download: nothing here reaches the network
vi.mock('../data-source.js', async (importOriginal) => ({
  ...(await importOriginal<typeof import('../data-source.js')>()),
  ajaxGet: vi.fn(),
}));
const { ajaxGet } = await import('../data-source.js');
const ajaxMock = vi.mocked(ajaxGet) as unknown as ReturnType<typeof vi.fn>;

const enc = new TextEncoder();

/** Answers every download with the bytes of `text`, after `gate` if given. */
function serve(text: string, gate?: Promise<void>) {
  ajaxMock.mockImplementation((params: any) =>
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
    ajaxMock.mockReset();
  });

  it('two resolves of the same URL fetch once', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const a = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    const b = await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    expect(`${a.data}`).toBe('abc');
    expect(`${b.data}`).toBe('abc');
    expect(ajaxMock).toHaveBeenCalledTimes(1);
    expect(m.urlCacheSize).toBe(1);
  });

  it('the same Asset.Url resolved twice counts two references', async () => {
    serve('abc');
    const m = new DefaultAssetManager();
    const asset = Asset.Url('https://x.org/a');
    const a = await m.resolve(asset, 'string').run();
    const b = await m.resolve(asset, 'string').run();
    expect(ajaxMock).toHaveBeenCalledTimes(1);
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
    expect(ajaxMock).toHaveBeenCalledTimes(1);
    b.dispose();
    c.dispose();
    expect(m.urlCacheSize).toBe(0);

    await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    expect(ajaxMock).toHaveBeenCalledTimes(2);
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
    expect(ajaxMock).toHaveBeenCalledTimes(1);

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
    expect(ajaxMock).toHaveBeenCalledTimes(2);
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
    expect(ajaxMock).toHaveBeenCalledTimes(1);

    serve('second'); // what the retry downloads
    abort('test');
    open();
    expect(await downloaderResult).toBe('aborted');
    const w = await waiter;
    expect(`${w.data}`).toBe('second');
    expect(ajaxMock).toHaveBeenCalledTimes(2);
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
    expect(ajaxMock).toHaveBeenCalledTimes(1);
  });

  it('keeps URLs and headers apart', async () => {
    serve('x');
    const m = new DefaultAssetManager();
    await m.resolve(Asset.Url('https://x.org/a'), 'string').run();
    await m.resolve(Asset.Url('https://x.org/b'), 'string').run();
    await m.resolve(Asset.Url('https://x.org/a', { headers: { Authorization: 'k' } }), 'string').run();
    await m.resolve(Asset.Url('https://x.org/a', { headers: { Authorization: 'k' } }), 'string').run();
    expect(ajaxMock).toHaveBeenCalledTimes(3);
    expect(m.urlCacheSize).toBe(3);
  });

  it('does not cache requests with a body, nor store: false', async () => {
    serve('x');
    const m = new DefaultAssetManager();
    const a = await m.resolve(Asset.Url('https://x.org/q', { body: '{}' }), 'string').run();
    const b = await m.resolve(Asset.Url('https://x.org/q', { body: '{}' }), 'string').run();
    expect(ajaxMock).toHaveBeenCalledTimes(2);
    expect(m.urlCacheSize).toBe(0);
    a.dispose();
    b.dispose();

    await m.resolve(Asset.Url('https://x.org/a'), 'string', false).run();
    await m.resolve(Asset.Url('https://x.org/a'), 'string', false).run();
    expect(ajaxMock).toHaveBeenCalledTimes(4);
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
