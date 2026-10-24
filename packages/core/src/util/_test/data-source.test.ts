/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import { gzipSync } from 'node:zlib';
import { describe, expect, it, vi } from 'vitest';
import { Task } from '@molstar/core/task';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { ajaxGet, ajaxGetMany, fetchData, fetchDataMany, readFromFile, setFSModule } from '../data-source.js';
import { zip } from '../zip/zip.js';

const encode = (s: string) => new TextEncoder().encode(s);

/** A response whose body streams `chunks`, one per read */
function streamed(chunks: Uint8Array[], headers?: HeadersInit, status = 200) {
  let i = 0;
  const body = new ReadableStream<Uint8Array>({
    pull(c) {
      if (i < chunks.length) c.enqueue(chunks[i++]);
      else c.close();
    },
  });
  return new Response(body, { status, headers });
}

function stubFetch(response: () => Response | Promise<Response>) {
  const fetch = vi.fn(async (_url: string, _init: RequestInit) => response());
  vi.stubGlobal('fetch', fetch);
  return fetch;
}

/** Runs the task with every progress update observed (an update rate of -1 notifies each one, with no clock) */
async function runObserved<T>(task: Task<T>) {
  const progress: Task.Progress[] = [];
  const result = await task.run((p) => progress.push(p.root.progress), -1);
  return { result, downloading: progress.filter((p) => p.message.startsWith('Downloading')) };
}

describe('ajaxGet', () => {
  it('is fetchData under its old name, as ajaxGetMany is fetchDataMany', () => {
    expect(ajaxGet).toBe(fetchData);
    expect(ajaxGetMany).toBe(fetchDataMany);
  });
});

describe('fetchData', () => {
  it('reads a string, by GET with no body', async () => {
    const fetch = stubFetch(() => streamed([encode('data_1cbs\n'), encode('_entry.id 1CBS €\n')]));
    expect(await fetchData('https://example.org/1cbs.cif').run()).toBe('data_1cbs\n_entry.id 1CBS €\n');
    // the request is let go once read
    expect(fetch.mock.calls[0][1].signal!.aborted).toBe(true);
    expect(fetch).toHaveBeenCalledWith('https://example.org/1cbs.cif', {
      signal: expect.any(AbortSignal),
      method: 'GET',
      body: undefined,
      headers: undefined,
    });
  });

  it('reads binary as bytes, without decompressing a .gz URL', async () => {
    const gz = gzipSync(encode('not decompressed'));
    stubFetch(() => streamed([gz.subarray(0, 5), gz.subarray(5)]));
    const data = await fetchData({ url: 'https://example.org/1cbs.bcif.gz', type: 'binary' }).run();
    expect(data).toBeInstanceOf(Uint8Array);
    expect(data).toEqual(new Uint8Array(gz));
  });

  it('parses JSON, by POST with the body and headers', async () => {
    const fetch = stubFetch(() => streamed([encode('{"ok":'), encode('true}')]));
    const headers = { 'Content-Type': 'application/json' };
    const result = await fetchData({ url: 'https://example.org/q', type: 'json', body: '{}', headers }).run();
    expect(result).toEqual({ ok: true });
    expect(fetch.mock.calls[0][1]).toMatchObject({ method: 'POST', body: '{}', headers });
  });

  it('gives null for a body that is not JSON, as XHR did', async () => {
    stubFetch(() => streamed([encode('<html>')]));
    expect(await fetchData({ url: 'https://example.org/q', type: 'json' }).run()).toBeNull();
  });

  it('parses XML from the decoded text', async () => {
    const parseFromString = vi.fn(() => 'document');
    vi.stubGlobal(
      'DOMParser',
      class {
        parseFromString = parseFromString;
      },
    );
    stubFetch(() => streamed([encode('<map><contour level="1.5"/>'), encode('</map>')]));
    expect(await fetchData({ url: 'https://example.org/h.xml', type: 'xml' }).run()).toBe('document');
    expect(parseFromString).toHaveBeenCalledWith('<map><contour level="1.5"/></map>', 'application/xml');
  });

  it('unzips a zip', async () => {
    const entries = { 'a.cif': encode('data_a'), 'b.cif': encode('data_b') };
    const zipped = new Uint8Array(await zip(SyncRuntimeContext, entries));
    stubFetch(() => streamed([zipped]));
    expect(await fetchData({ url: 'https://example.org/x.zip', type: 'zip' }).run()).toEqual(entries);
  });

  it('reports progress against Content-Length', async () => {
    const chunks = [encode('aaaa'), encode('bbbb'), encode('cc')];
    stubFetch(() => streamed(chunks, { 'Content-Length': '10' }));
    const { result, downloading } = await runObserved(fetchData({ url: 'https://example.org/x', type: 'binary' }));
    expect(result).toEqual(encode('aaaabbbbcc'));
    expect(downloading.map((p) => [p.current, p.max, p.isIndeterminate])).toEqual([
      [4, 10, false],
      [8, 10, false],
      [10, 10, false],
    ]);
  });

  it('reports bytes only with no Content-Length', async () => {
    stubFetch(() => streamed([new Uint8Array(1024 * 1024), new Uint8Array(1024 * 1024)]));
    const { downloading } = await runObserved(fetchData({ url: 'https://example.org/x', type: 'binary' }));
    expect(downloading.map((p) => [p.message, p.isIndeterminate])).toEqual([
      ['Downloading... 1.00 MB', true],
      ['Downloading... 2.00 MB', true],
    ]);
  });

  it('reports bytes only for a gzip-encoded body, which streams decoded bytes', async () => {
    // fetch decodes the body; Content-Length counts the encoded bytes
    stubFetch(() => streamed([encode('aaaa'), encode('bbbb')], { 'Content-Encoding': 'gzip', 'Content-Length': '5' }));
    const { result, downloading } = await runObserved(fetchData({ url: 'https://example.org/x', type: 'string' }));
    expect(result).toBe('aaaabbbb');
    expect(downloading.every((p) => p.isIndeterminate)).toBe(true);
  });

  it('keeps every byte of a body longer than its Content-Length, and of a shorter one', async () => {
    stubFetch(() => streamed([encode('aaaa'), encode('bbbb')], { 'Content-Length': '6' }));
    const longer = await runObserved(fetchData({ url: 'https://example.org/x', type: 'string' }));
    expect(longer.result).toBe('aaaabbbb');
    expect(longer.downloading.map((p) => p.isIndeterminate)).toEqual([false, true]);

    stubFetch(() => streamed([encode('aaaa')], { 'Content-Length': '6' }));
    expect(await fetchData({ url: 'https://example.org/x', type: 'string' }).run()).toBe('aaaa');
  });

  it('fails on an HTTP error with its status code', async () => {
    const fetch = stubFetch(() => streamed([encode('Not Found')], {}, 404));
    await expect(fetchData('https://example.org/missing').run()).rejects.toThrow(
      'Download failed with status code 404',
    );
    // the unread body is cancelled
    expect(fetch.mock.calls[0][1].signal!.aborted).toBe(true);
  });

  it('fails on a network error with the error Node.js gives', async () => {
    // a browser hides the cause: there it is the generic message ("Failed to download data. ..."), as with XHR
    const refused = new TypeError('fetch failed');
    stubFetch(() => Promise.reject(refused));
    await expect(fetchData('https://example.org/x').run()).rejects.toBe(refused);

    const terminated = new TypeError('terminated');
    const broken = new ReadableStream<Uint8Array>({
      pull(c) {
        c.error(terminated);
      },
    });
    stubFetch(() => new Response(broken));
    await expect(fetchData('https://example.org/x').run()).rejects.toBe(terminated);
  });

  it('stops reading when the task is aborted', async () => {
    let pulls = 0;
    const endless = new ReadableStream<Uint8Array>({
      pull(c) {
        pulls++;
        c.enqueue(new Uint8Array(16));
      },
    });
    const fetch = stubFetch(() => new Response(endless));
    const run = fetchData({ url: 'https://example.org/x', type: 'binary' }).run((p) => {
      if (p.root.progress.message.startsWith('Downloading')) p.requestAbort();
    }, -1);
    await expect(run).rejects.toSatisfy(Task.isAbort);
    // the first chunk's progress update sees the abort; the stream's queue holds one more
    expect(pulls).toBeLessThanOrEqual(2);
    expect(fetch.mock.calls[0][1].signal!.aborted).toBe(true);
  });

  it("rejects as aborted when the caller's signal aborts while waiting for the server", async () => {
    // as fetch does: the request rejects with the signal's reason, at once when it is aborted already
    const waiting = vi.fn((_url: string, init: RequestInit) => {
      const signal = init.signal!;
      return new Promise<Response>((_, reject) => {
        if (signal.aborted) reject(signal.reason);
        else signal.addEventListener('abort', () => reject(signal.reason));
      });
    });
    vi.stubGlobal('fetch', waiting);
    const controller = new AbortController();
    const run = fetchData({ url: 'https://example.org/x', signal: controller.signal }).run();
    controller.abort('closed');
    await expect(run).rejects.toSatisfy((e) => Task.isAbort(e) && e.reason === 'closed');
    expect(waiting).toHaveBeenCalledTimes(1);
  });

  it("rejects as aborted when the caller's signal aborts a stalled body", async () => {
    // as fetch does: an aborted request errors its body
    let signal: AbortSignal | undefined;
    let pulls = 0;
    const stalled = new ReadableStream<Uint8Array>({
      pull(c) {
        if (pulls++ === 0) c.enqueue(encode('aaaa'));
        else if (signal!.aborted) c.error(signal!.reason);
        else return new Promise<void>(() => signal!.addEventListener('abort', () => c.error(signal!.reason)));
      },
    });
    vi.stubGlobal(
      'fetch',
      vi.fn(async (_url: string, init: RequestInit) => {
        signal = init.signal!;
        return new Response(stalled);
      }),
    );
    const controller = new AbortController();
    const run = fetchData({ url: 'https://example.org/x', signal: controller.signal }).run((p) => {
      if (p.root.progress.message.startsWith('Downloading')) controller.abort();
    }, -1);
    await expect(run).rejects.toSatisfy((e) => Task.isAbort(e) && e.reason === 'This operation was aborted');
  });

  it("rejects as aborted, with no request, when the caller's signal is aborted already", async () => {
    const fetch = stubFetch(() => streamed([]));
    const run = fetchData({ url: 'https://example.org/x', signal: AbortSignal.abort('gone') }).run();
    await expect(run).rejects.toSatisfy((e) => Task.isAbort(e) && e.reason === 'gone');
    expect(fetch).not.toHaveBeenCalled();
  });

  it("rejects as aborted when the caller's signal aborts as an HTTP error arrives", async () => {
    const controller = new AbortController();
    stubFetch(() => {
      controller.abort('late');
      return streamed([], {}, 404);
    });
    const run = fetchData({ url: 'https://example.org/x', signal: controller.signal }).run();
    await expect(run).rejects.toSatisfy((e) => Task.isAbort(e) && e.reason === 'late');
  });

  it("lets go of the caller's signal after each download", async () => {
    // a long-lived signal reused for many downloads must not collect a listener per download
    const controller = new AbortController();
    const add = vi.spyOn(controller.signal, 'addEventListener');
    const remove = vi.spyOn(controller.signal, 'removeEventListener');
    stubFetch(() => streamed([encode('a')]));
    await fetchData({ url: 'https://example.org/x', signal: controller.signal }).run();
    stubFetch(() => streamed([], {}, 500));
    await expect(fetchData({ url: 'https://example.org/x', signal: controller.signal }).run()).rejects.toThrow();
    expect(add).toHaveBeenCalledTimes(2);
    expect(remove.mock.calls).toEqual(add.mock.calls.map(([type, listener]) => [type, listener]));
  });

  it('reads a body with a Content-Length too large to allocate in chunks', async () => {
    stubFetch(() => streamed([encode('aaaa')], { 'Content-Length': String(Number.MAX_SAFE_INTEGER) }));
    expect(await fetchData('https://example.org/x').run()).toBe('aaaa');
  });

  it('reads a file:// URL in Node.js from the file system, not by fetch', async () => {
    const fetch = stubFetch(() => streamed([]));
    // a buffer of its own, as fs.readFile gives: not a view into the pool Buffer.from uses for small strings
    const local = Buffer.alloc(5, 'local');
    setFSModule({ readFile: (_: string, cb: (e: null, d: Buffer) => void) => cb(null, local) } as any);
    try {
      expect(await fetchData('file:///tmp/x.cif').run()).toBe('local');
      expect(fetch).not.toHaveBeenCalled();
    } finally {
      setFSModule(fs);
    }
  });
});

describe('readFromFile', () => {
  it('decompresses a .gz file and a single-entry .zip file', async () => {
    const gz = new File([gzipSync(encode('data_gz'))], 'x.cif.gz');
    expect(await readFromFile(gz, 'string').run()).toBe('data_gz');
    const zipped = new File([await zip(SyncRuntimeContext, { 'x.cif': encode('data_zip') })], 'x.zip');
    expect(await readFromFile(zipped, 'string').run()).toBe('data_zip');
  });
});
