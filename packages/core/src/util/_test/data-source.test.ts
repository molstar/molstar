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
import { ajaxGet, readFromFile, setFSModule } from '../data-source.js';
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
  it('reads a string, by GET with no body', async () => {
    const fetch = stubFetch(() => streamed([encode('data_1cbs\n'), encode('_entry.id 1CBS €\n')]));
    expect(await ajaxGet('https://example.org/1cbs.cif').run()).toBe('data_1cbs\n_entry.id 1CBS €\n');
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
    const data = await ajaxGet({ url: 'https://example.org/1cbs.bcif.gz', type: 'binary' }).run();
    expect(data).toBeInstanceOf(Uint8Array);
    expect(data).toEqual(new Uint8Array(gz));
  });

  it('parses JSON, by POST with the body and headers', async () => {
    const fetch = stubFetch(() => streamed([encode('{"ok":'), encode('true}')]));
    const headers = { 'Content-Type': 'application/json' };
    const result = await ajaxGet({ url: 'https://example.org/q', type: 'json', body: '{}', headers }).run();
    expect(result).toEqual({ ok: true });
    expect(fetch.mock.calls[0][1]).toMatchObject({ method: 'POST', body: '{}', headers });
  });

  it('gives null for a body that is not JSON, as XHR did', async () => {
    stubFetch(() => streamed([encode('<html>')]));
    expect(await ajaxGet({ url: 'https://example.org/q', type: 'json' }).run()).toBeNull();
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
    expect(await ajaxGet({ url: 'https://example.org/h.xml', type: 'xml' }).run()).toBe('document');
    expect(parseFromString).toHaveBeenCalledWith('<map><contour level="1.5"/></map>', 'application/xml');
  });

  it('unzips a zip', async () => {
    const entries = { 'a.cif': encode('data_a'), 'b.cif': encode('data_b') };
    const zipped = new Uint8Array(await zip(SyncRuntimeContext, entries));
    stubFetch(() => streamed([zipped]));
    expect(await ajaxGet({ url: 'https://example.org/x.zip', type: 'zip' }).run()).toEqual(entries);
  });

  it('reports progress against Content-Length', async () => {
    const chunks = [encode('aaaa'), encode('bbbb'), encode('cc')];
    stubFetch(() => streamed(chunks, { 'Content-Length': '10' }));
    const { result, downloading } = await runObserved(ajaxGet({ url: 'https://example.org/x', type: 'binary' }));
    expect(result).toEqual(encode('aaaabbbbcc'));
    expect(downloading.map((p) => [p.current, p.max, p.isIndeterminate])).toEqual([
      [4, 10, false],
      [8, 10, false],
      [10, 10, false],
    ]);
  });

  it('reports bytes only with no Content-Length', async () => {
    stubFetch(() => streamed([new Uint8Array(1024 * 1024), new Uint8Array(1024 * 1024)]));
    const { downloading } = await runObserved(ajaxGet({ url: 'https://example.org/x', type: 'binary' }));
    expect(downloading.map((p) => [p.message, p.isIndeterminate])).toEqual([
      ['Downloading... 1.00 MB', true],
      ['Downloading... 2.00 MB', true],
    ]);
  });

  it('reports bytes only for a gzip-encoded body, which streams decoded bytes', async () => {
    // fetch decodes the body; Content-Length counts the encoded bytes
    stubFetch(() => streamed([encode('aaaa'), encode('bbbb')], { 'Content-Encoding': 'gzip', 'Content-Length': '5' }));
    const { result, downloading } = await runObserved(ajaxGet({ url: 'https://example.org/x', type: 'string' }));
    expect(result).toBe('aaaabbbb');
    expect(downloading.every((p) => p.isIndeterminate)).toBe(true);
  });

  it('keeps every byte of a body longer than its Content-Length, and of a shorter one', async () => {
    stubFetch(() => streamed([encode('aaaa'), encode('bbbb')], { 'Content-Length': '6' }));
    const longer = await runObserved(ajaxGet({ url: 'https://example.org/x', type: 'string' }));
    expect(longer.result).toBe('aaaabbbb');
    expect(longer.downloading.map((p) => p.isIndeterminate)).toEqual([false, true]);

    stubFetch(() => streamed([encode('aaaa')], { 'Content-Length': '6' }));
    expect(await ajaxGet({ url: 'https://example.org/x', type: 'string' }).run()).toBe('aaaa');
  });

  it('fails on an HTTP error with its status code', async () => {
    stubFetch(() => streamed([encode('Not Found')], {}, 404));
    await expect(ajaxGet('https://example.org/missing').run()).rejects.toThrow('Download failed with status code 404');
  });

  it('fails on a network error with the error Node.js gives', async () => {
    // a browser hides the cause: there it is the generic message ("Failed to download data. ..."), as with XHR
    const refused = new TypeError('fetch failed');
    stubFetch(() => Promise.reject(refused));
    await expect(ajaxGet('https://example.org/x').run()).rejects.toBe(refused);

    const terminated = new TypeError('terminated');
    const broken = new ReadableStream<Uint8Array>({
      pull(c) {
        c.error(terminated);
      },
    });
    stubFetch(() => new Response(broken));
    await expect(ajaxGet('https://example.org/x').run()).rejects.toBe(terminated);
  });

  it('stops reading when the task is aborted', async () => {
    let pulls = 0;
    const endless = new ReadableStream<Uint8Array>({
      pull(c) {
        pulls++;
        c.enqueue(new Uint8Array(16));
      },
    });
    stubFetch(() => new Response(endless));
    const run = ajaxGet({ url: 'https://example.org/x', type: 'binary' }).run((p) => {
      if (p.root.progress.message.startsWith('Downloading')) p.requestAbort();
    }, -1);
    await expect(run).rejects.toSatisfy(Task.isAbort);
    // the first chunk's progress update sees the abort; the stream's queue holds one more
    expect(pulls).toBeLessThanOrEqual(2);
  });

  it('reads a body with a Content-Length too large to allocate in chunks', async () => {
    stubFetch(() => streamed([encode('aaaa')], { 'Content-Length': String(Number.MAX_SAFE_INTEGER) }));
    expect(await ajaxGet('https://example.org/x').run()).toBe('aaaa');
  });

  it('reads a file:// URL in Node.js from the file system, not by fetch', async () => {
    const fetch = stubFetch(() => streamed([]));
    // a buffer of its own, as fs.readFile gives: not a view into the pool Buffer.from uses for small strings
    const local = Buffer.alloc(5, 'local');
    setFSModule({ readFile: (_: string, cb: (e: null, d: Buffer) => void) => cb(null, local) } as any);
    try {
      expect(await ajaxGet('file:///tmp/x.cif').run()).toBe('local');
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
