/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for LY1's hot paths (seqspec:viewer#performance), on `heapScaling` (a step is one move or one frame's
// reads): hit testing per pointer move over every region of cases 3, 6 and 9, a wrapped section with inline labels, a
// reversed and a cyclic piece, split rows and the tree; layout lookups per frame with the view moving; fits, numbers
// and cached labels; `layoutRuns` on resize; a wrapped line's lanes over 5,000 features on 35,000 residues (the cases
// in layout-cases.ts, shared with the quick twin, layout.quick.ts). The positive control: hit testing's moves
// allocating an object each (`leaking`), through the same measurement, are caught.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { fits, hitPoints, laneHeights, lookups, runsOnResize } from './layout-cases.js';

const N = 50_000;

describe("LY1's hot paths allocate nothing", () => {
  it('hit testing per pointer move: cases 3, 6, 9 and 13 (split rows), wrapped with inline labels, reversed and cyclic pieces, the tree', async () => {
    const { step, kinds, nodes } = await hitPoints();
    expect([kinds > 10, nodes > 10]).toEqual([true, true]);
    const r = heapScaling(N, step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  }, 120_000);

  it('layout lookups per frame: the layout, the axis view and its x, the visible rows and entries, with the view moving', async () => {
    const r = heapScaling(N, await lookups());
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('fits, numbers and cached labels', async () => {
    const r = heapScaling(N, await fits());
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("layoutRuns on resize (case 1's text view, and 35,000 residues)", async () => {
    const { step, runs } = await runsOnResize();
    expect(runs).toEqual([1, 35_000]);
    const r = heapScaling(2_000, step);
    expect(isFlat(r), describeScaling(r, 2_000)).toBe(true);
  });

  it("a wrapped line's lanes: 5,000 features over 35,000 residues, every line's band", async () => {
    const { step, lines, segments } = await laneHeights();
    expect([lines > 200, segments > 1]).toEqual([true, true]);
    // 600 steps: at 200 the template's control leaked under the scaled limit, blind (LC-19).
    const r = heapScaling(600, step);
    expect(isFlat(r), describeScaling(r, 600)).toBe(true);
    expect(controlCaught(() => heapScaling(600, leaking(step)))).toBe(true);
  });
});
