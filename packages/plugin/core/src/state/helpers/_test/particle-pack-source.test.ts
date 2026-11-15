/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { getUrlRangeSource } from '../particle-pack-source.js';

/** a server of `bytes` answering ranges (206), a range past its end or reversed refused (416) */
function stubRangeServer(bytes: Uint8Array) {
  const ranges: string[] = [];
  vi.stubGlobal('fetch', async (_url: string, init?: RequestInit) => {
    const range = (init?.headers as Record<string, string>)?.Range ?? '';
    ranges.push(range);
    const m = /^bytes=(\d+)-(\d+)$/.exec(range);
    const start = m ? Number(m[1]) : NaN;
    const end = m ? Number(m[2]) : NaN;
    if (!(start <= end) || end >= bytes.byteLength) return new Response(null, { status: 416 });
    return new Response(bytes.slice(start, end + 1), {
      status: 206,
      headers: { 'content-range': `bytes ${start}-${end}/${bytes.byteLength}` },
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
});
