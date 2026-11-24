/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { DefaultAssetManager } from '@molstar/core/util/assets';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { AssetRefused, BrowserAssetPolicy, headlessAssetPolicy, policyAssets } from '../assets.js';

// The policy's range reads (_spec/README.md "Context", Ports; core's `util/_spec/assets.md` §10): the same rules as
// `readBytes` (the location, the credentials and redirect rule on every request, the size bound, the time), over core's
// manager and a fake `fetch` that honours ranges. Nothing reaches the network.

const Size = 300_000;
const body = Uint8Array.from({ length: Size }, (_, i) => (i * 13) & 255);

function serve() {
  const log: { url: string; init: RequestInit }[] = [];
  vi.stubGlobal('fetch', async (input: RequestInfo | URL, init: RequestInit = {}) => {
    log.push({ url: String(input), init });
    const m = /^bytes=(\d+)-(\d+)$/.exec(new Headers(init.headers).get('Range') ?? '');
    if (!m) return new Response(body.slice(), { status: 200 });
    const a = Number(m[1]);
    const b = Math.min(Number(m[2]), Size - 1);
    return new Response(body.slice(a, b + 1), {
      status: 206,
      headers: { 'Content-Range': `bytes ${a}-${b}/${Size}` },
    });
  });
  return log;
}

afterEach(() => vi.unstubAllGlobals());

const signal = new AbortController().signal;

describe('policyAssets.openRange', () => {
  const policy = headlessAssetPolicy({ hosts: ['data.example.org'] });

  it('reads by range under the policy: a page, not the file; its credentials and redirect rule on every request', async () => {
    const log = serve();
    const reader = policyAssets(new DefaultAssetManager(), policy);
    const r = await reader.openRange!('https://data.example.org/a.bw', signal);
    expect(r.size).toBe(Size);
    const [x] = await r.read([{ offset: 200_000, length: 10 }], signal);
    expect(x).toEqual(body.slice(200_000, 200_010));
    expect(log.length).toBe(2);
    for (const l of log) {
      expect(new Headers(l.init.headers).get('Range')).toMatch(/^bytes=/);
      expect(l.init.credentials).toBe('omit');
      expect(l.init.redirect).toBe('error');
    }
    r.dispose();
  });

  it('a location the policy refuses is refused before any request', async () => {
    const log = serve();
    const reader = policyAssets(new DefaultAssetManager(), policy);
    for (const uri of ['https://other.org/a.bw', 'http://data.example.org/a.bw', 'file:///etc/passwd'])
      await expect(reader.openRange!(uri, signal), uri).rejects.toBeInstanceOf(AssetRefused);
    expect(log).toEqual([]);
  });

  it('an asset over the policy’s size is refused; the composition’s credentials are passed', async () => {
    serve();
    const small = policyAssets(new DefaultAssetManager(), { ...policy, maxBytes: Size - 1 });
    await expect(small.openRange!('https://data.example.org/a.bw', signal)).rejects.toThrow('larger than the limit');
    const log = serve();
    const lab = policyAssets(new DefaultAssetManager(), {
      ...BrowserAssetPolicy,
      schemes: ['https:'],
      credentials: 'include',
    });
    (await lab.openRange!('https://lab.example.org/a.bw', signal)).dispose();
    expect(log[0].init.credentials).toBe('include');
  });

  it('a read that takes longer than the policy’s time is refused', async () => {
    serve();
    const reader = policyAssets(new DefaultAssetManager(), { ...policy, timeoutMs: 20 });
    const r = await reader.openRange!('https://data.example.org/a.bw', signal);
    vi.stubGlobal('fetch', (_: unknown, init: RequestInit) => {
      return new Promise((_, reject) => init.signal?.addEventListener('abort', () => reject(init.signal!.reason)));
    });
    await expect(r.read([{ offset: 100_000, length: 4 }], signal)).rejects.toThrow('took too long');
    r.dispose();
  });
});
