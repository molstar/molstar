/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { ParticleList } from '../particle-list.js';
import { type DeferredRead, DeferredPageSize, ParticlePageCache } from '../particle-page-cache.js';

// The pages read of deferred lists: a read once while in flight, kept within the cache's bytes by use, a failure not kept

const list = (label: string, count: number) => ({ label, count }) as ParticleList;
/** a page of 8 bytes a particle (a grid's two words) */
const pageBytes = DeferredPageSize * 8;

describe('ParticlePageCache', () => {
  it('reads a page once, keeps the most recently used within its bytes, a failed read not kept', async () => {
    const reads: string[] = [];
    let fail = false;
    const reader =
      (label: string): DeferredRead =>
      async (first, count) => {
        reads.push(`${label}:${first}:${count}`);
        if (fail) throw new Error('offline');
        return { positions: new Uint32Array(count * 2) };
      };
    const cache = new ParticlePageCache(2 * pageBytes);
    const a = list('a', 2 * DeferredPageSize + 10);
    const b = list('b', DeferredPageSize);
    const da = cache.deferred(a, reader('a'));
    const db = cache.deferred(b, reader('b'));
    // the last page shorter; two lists' pages apart
    expect((await da.readPage(2)).positions.length).toBe(20);
    await da.readPage(0);
    expect(cache.bytes).toBe(pageBytes + 80);
    // the first page used: the last let go when a third is read past the bytes
    expect(da.getPage(2)).toBeDefined();
    await db.readPage(0);
    expect(da.getPage(0)).toBeUndefined();
    expect(da.getPage(2)).toBeDefined();
    expect(db.getPage(0)).toBeDefined();
    expect(reads).toEqual([`a:${2 * DeferredPageSize}:10`, 'a:0:8192', 'b:0:8192']);
    expect(cache.readBytes).toBe(2 * pageBytes + 80);
    // a failure not kept: read again when asked
    fail = true;
    await expect(da.readPage(1)).rejects.toThrow('offline');
    expect(da.getPage(1)).toBeUndefined();
    fail = false;
    await da.readPage(1);
    expect(reads.filter((r) => r === `a:${DeferredPageSize}:8192`).length).toBe(2);
    await expect(da.readPage(3)).rejects.toThrow(/no page/);
  });
});
