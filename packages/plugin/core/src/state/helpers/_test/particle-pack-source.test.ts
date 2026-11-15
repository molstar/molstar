/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { ParticleList } from '@molstar/model/model/particles/particle-list';
import { readParticlePackHeader, writeParticlePack } from '@molstar/model/formats/particles/pack';
import { openParticlePackStream } from '@molstar/model/formats/particles/pack-stream';
import { getUrlRangeSource } from '../particle-pack-source.js';

/** a small list of one target (sites) */
function createList(): ParticleList {
  const n = 50;
  return {
    count: n,
    targets: new Int32Array(n).fill(1),
    targetInfo: new Map([[1, { name: 'one' }]]),
    targetMapping: new Map([
      [1, { kind: 'sites', sites: { count: 1, positions: new Float32Array(3), radii: new Float32Array([1]) } }],
    ]),
    coordinates: Float32Array.from({ length: n * 3 }, (_, i) => i % 17),
    getParticleLabel: (index) => `${index}`,
    sourceData: { kind: 'test', name: 'test', data: {} },
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
  };
}

/**
 * a server of `bytes` answering ranges (206), a range past its end or reversed refused (416); `exposed` false: of
 * another origin not exposing `Content-Range` (its total unread), a HEAD its `Content-Length` (`head` false: none)
 */
function stubRangeServer(bytes: Uint8Array, options: { exposed?: boolean; head?: boolean } = {}) {
  const { exposed = true, head = true } = options;
  const ranges: string[] = [];
  vi.stubGlobal('fetch', async (_url: string, init?: RequestInit) => {
    if (init?.method === 'HEAD') {
      return new Response(null, { status: 200, headers: head ? { 'content-length': `${bytes.byteLength}` } : {} });
    }
    const range = (init?.headers as Record<string, string>)?.Range ?? '';
    ranges.push(range);
    const m = /^bytes=(\d+)-(\d+)$/.exec(range);
    const start = m ? Number(m[1]) : NaN;
    const end = m ? Number(m[2]) : NaN;
    if (!(start <= end) || end >= bytes.byteLength) return new Response(null, { status: 416 });
    return new Response(bytes.slice(start, end + 1), {
      status: 206,
      headers: exposed ? { 'content-range': `bytes ${start}-${end}/${bytes.byteLength}` } : {},
    });
  });
  return ranges;
}

describe('a particle pack read by HTTP ranges', () => {
  it('reads a range as asked, an empty one without a request', async () => {
    const bytes = Uint8Array.from({ length: 100 }, (_, i) => i);
    const ranges = stubRangeServer(bytes);
    const source = (await getUrlRangeSource('https://example.org/pack.mspk'))!;
    expect(source.byteLength).toBe(100);
    expect(Array.from(await source.read(10, 3))).toEqual([10, 11, 12]);
    const asked = ranges.length;
    expect((await source.read(40, 0)).byteLength).toBe(0);
    expect(ranges.length).toBe(asked);
  });

  it("takes the size by a HEAD where another origin's 'Content-Range' is not exposed; says why where neither", async () => {
    const bytes = Uint8Array.from({ length: 100 }, (_, i) => i);
    stubRangeServer(bytes, { exposed: false });
    const warned: string[] = [];
    const source = (await getUrlRangeSource('https://example.org/pack.mspk', undefined, (m) => warned.push(m)))!;
    expect(source.byteLength).toBe(100);
    expect(Array.from(await source.read(98, 2))).toEqual([98, 99]);
    expect(warned).toEqual([]);

    stubRangeServer(bytes, { exposed: false, head: false });
    expect(await getUrlRangeSource('https://example.org/pack.mspk', undefined, (m) => warned.push(m))).toBeUndefined();
    expect(warned).toHaveLength(1);
    expect(warned[0]).toMatch(/Access-Control-Expose-Headers: Content-Range/);
  });

  it('takes the size from the pack itself where no header gives it (a presigned URL), its document read once', async () => {
    const bytes = writeParticlePack(createList(), { tiles: 2 });
    const ranges = stubRangeServer(bytes, { exposed: false, head: false });
    const warned: string[] = [];
    const source = (await getUrlRangeSource('https://example.org/pack.mspk', undefined, (m) => warned.push(m)))!;
    expect(warned).toEqual([]);
    expect(source.byteLength).toBe(bytes.byteLength);
    // opened as by its size: the document's range asked once (for the size), its read answered with it
    const { jsonOffset, jsonLength } = readParticlePackHeader(bytes.subarray(0, 64));
    const stream = await openParticlePackStream(source);
    expect(stream.document.count).toBe(50);
    expect(ranges.filter((r) => r === `bytes=${jsonOffset}-${jsonOffset + jsonLength - 1}`)).toHaveLength(1);
  });

  it('says so where a server ignores ranges (a whole answer): read whole', async () => {
    vi.stubGlobal('fetch', async () => new Response(new Uint8Array(100), { status: 200 }));
    const warned: string[] = [];
    expect(await getUrlRangeSource('https://example.org/pack.mspk', undefined, (m) => warned.push(m))).toBeUndefined();
    expect(warned[0]).toMatch(/answers no range requests/);
  });
});
