/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for LY2's per-frame and per-move reads (seqspec:viewer#performance), on `heapScaling` (a step is one
// read): segment windows and the item count that picks items or density, the hit's `featureAt` over overlapping
// fragments, `ofPartAt` stepping and jumping from its cursor, the level choice, a column's reads of class, value and
// count pyramids, the cell key and the mode's hysteresis. Each is what a draw or a hit runs per column, item or move
// (the cases in prepare-cases.ts, shared with the quick twin, prepare.quick.ts). The positive control: the windows'
// reads allocating an object each (`leaking`), measured alike, are caught.

import { beforeAll, describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { LENGTH, prepareReads } from './prepare-cases.js';
import { itemsIn } from './prepare-fixtures.js';

const N = 50_000;
let c: Awaited<ReturnType<typeof prepareReads>>;

beforeAll(async () => {
  c = await prepareReads();
}, 120_000);

describe("LY2's reads allocate nothing", () => {
  it('segment windows and the count in view (items or density) as the view pans and zooms', () => {
    expect(itemsIn(c.genes, 0, 0, LENGTH, -1, 1e9)).toBeGreaterThan(100);
    const r = heapScaling(N, c.windows);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("the hit's `featureAt` per move, every lane, and `positionKey`", () => {
    const r = heapScaling(N, c.features);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('`ofPartAt` stepping along a frame and jumping between frames', () => {
    expect(c.of.count).toBeGreaterThan(10_000);
    const r = heapScaling(N, c.ofParts);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("a device column's reads: its plurality class, its values, its counts", () => {
    const r = heapScaling(N, c.columns);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("the windows' reads allocating an object each are caught (positive control)", () => {
    expect(controlCaught(() => heapScaling(N, leaking(c.windows)))).toBe(true);
  });
});
