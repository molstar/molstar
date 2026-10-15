/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the channels' per-move path (IN3; seqspec:interaction#channels: after a route's first use nothing
// allocates; CP2's `spans.heap`), on RB0b's stepped `heapScaling`: a hover element change per step, set through the
// hover manager, then the channels' rebuild (closed, projected into every section's hover set, its primary) and the
// overlay's reads; and a million element changes over 10,000 distinct alignment rows and their positions (a row's
// residue on its own columns read without the row's cached map, which allocates on a miss). The steps are
// channels-cases.ts's, shared with the quick twin (channels.quick.ts). The positive control: each step allocating an
// object too (`leaking`), measured alike, is caught.

import { describe, expect, it } from 'vitest';
import { alignment, genes, ROWS, sunk } from './channels-cases.js';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';

describe("the channels' hover allocates nothing per element change", () => {
  it('a nucleotide or a residue hovered, closed, projected into every section with its primary, and read', async () => {
    const { step, ch } = await genes();
    const result = heapScaling(50_000, step);
    expect(ch.rebuilds.hover).toBeGreaterThan(200_000);
    expect(sunk()).not.toBe(0);
    expect(isFlat(result), describeScaling(result, 50_000)).toBe(true);
  });

  it(`a million element changes over ${ROWS} distinct rows and positions of an alignment`, async () => {
    const { step, ch } = await alignment();
    const result = heapScaling(250_000, step);
    expect(ch.rebuilds.hover).toBeGreaterThan(1_000_000);
    expect(sunk()).not.toBe(0);
    expect(isFlat(result), describeScaling(result, 250_000)).toBe(true);
  });

  it('the same steps allocating an object each are caught (positive control)', async () => {
    const { step } = await genes();
    expect(controlCaught(() => heapScaling(50_000, leaking(step)))).toBe(true);
  });
});
