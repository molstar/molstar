/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The reference server holds nothing (ST2): a million rows are sampled when a tile is asked, never generated and
// kept (the first version cached every row string, about 4 GB at that size). The million is the heavy lane's; the quick
// lane asks the same of 20,000 rows.

import { TILE } from '@molstar/seq-core/format';
import { describe, expect, it } from 'vitest';
import { tileClient } from '../../../src/stream/tiles/index.js';
import { proceduralCode, proceduralRow, tileServer } from '../../../src/stream/tiles/server.js';

describe('20,000 rows', () => {
  const [R, C] = [20_000, 320];
  const server = tileServer({ rowCount: R, columnCount: C });
  const client = tileClient('x', server.transport);

  it('answers the tiles of every level, sampled', async () => {
    // Level 7's samples are 128 rows apart: the 157th starts past the last row.
    const top = await client.tile(7, 1, 0, 0);
    const mid = await client.tile(4, 0, 2, 0);
    const fine = await client.tile(0, 0, 70, 1);
    expect([top.cls[0] > 0, top.cls[156 * TILE] > 0, top.cls[157 * TILE]]).toEqual([true, true, 0]);
    expect(mid.cls.every((c, i) => (i < 256 * TILE ? c > 0 : true))).toBe(true);
    expect(String.fromCharCode(fine.cls[5 * TILE + 3])).toBe(
      String.fromCharCode(proceduralCode(70 * 256 + 5, 256 + 3)),
    );
  });

  it('the same tile is the same bytes, and a row is the row it would be', async () => {
    const [a, b] = [await client.tile(5, 0, 0, 0), await client.tile(5, 0, 0, 0)];
    expect(Array.from(a.cls.subarray(0, 64))).toEqual(Array.from(b.cls.subarray(0, 64)));
    expect((await client.row(19_999)).row).toBe(proceduralRow(19_999, C));
    expect((await client.labels(19_998, 20_000)).map((l) => l.split('/')[0])).toEqual(['SEQ19998', 'SEQ19999']);
    expect(await client.stats('conservation', 0, 0)).toBeNull();
    expect(server.calls.length).toBe(8);
  });
});

describe('a million rows', { tags: ['heavy'] }, () => {
  const [R, C] = [1_000_000, 320];
  const server = tileServer({ rowCount: R, columnCount: C });
  const client = tileClient('x', server.transport);

  it('answers the tiles of every level in a blink, in memory that does not grow', async () => {
    const before = process.memoryUsage().heapUsed;
    const t0 = performance.now();
    const top = await client.tile(15, 1, 0, 0);
    const mid = await client.tile(10, 0, 2, 0);
    const fine = await client.tile(0, 0, 3906, 1);
    const ms = performance.now() - t0;
    // Samples past their rows hold no data; those inside hold a residue at a share of the cells read.
    expect([top.cls[0] > 0, top.cls[30 * TILE] > 0, top.cls[31 * TILE]]).toEqual([true, true, 0]);
    expect(mid.cls.every((c, i) => (i < 256 * TILE ? c > 0 : true))).toBe(true);
    expect(String.fromCharCode(fine.cls[5 * TILE + 3])).toBe(
      String.fromCharCode(proceduralCode(3906 * 256 + 5, 256 + 3)),
    );
    expect(ms).toBeLessThan(5000);
    expect(process.memoryUsage().heapUsed - before).toBeLessThan(100 * 2 ** 20);
  });

  it('the same tile is the same bytes, and a row is the row it would be', async () => {
    const [a, b] = [await client.tile(12, 0, 0, 0), await client.tile(12, 0, 0, 0)];
    expect(Array.from(a.cls.subarray(0, 64))).toEqual(Array.from(b.cls.subarray(0, 64)));
    expect((await client.row(999_999)).row).toBe(proceduralRow(999_999, C));
    expect((await client.labels(999_998, 1_000_000)).map((l) => l.split('/')[0])).toEqual(['SEQ999998', 'SEQ999999']);
    expect(await client.stats('conservation', 0, 0)).toBeNull();
    expect(server.calls.length).toBe(8);
  });
});
