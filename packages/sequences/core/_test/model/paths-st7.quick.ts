/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of paths-st7.heap.ts (`test:perf`): a step of reads of a sparse composed path (paths-st7-cases.ts, the
// guard's cases), with its positive control.

import { describe, expect, it } from 'vitest';
import { down, read, sums, up } from './paths-st7-cases.js';
import { quickCase } from './quick.js';

describe("ST7: a sparse path's reads, quick", () => {
  it('Paths.at, spanAt, forEachFrom and forEachTo over a 79-exon path through 200 Mb', async ({ task }) => {
    quickCase(task, 'paths-st7: at, spanAt, forEachFrom and forEachTo', read);
    expect([down.length, down.heldTo.length < 12_000, up.heldSpans?.length]).toEqual([200_000_000, true, 3 * 3685]);
    expect(Number.isFinite(sums())).toBe(true);
    expect(sums()).toBeGreaterThan(0);
  });
});
