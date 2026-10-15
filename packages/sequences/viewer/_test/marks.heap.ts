/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the marks and the projector (IN1; seqspec:interaction#channels: after a route's first use nothing
// allocates; a hover change writes a few entries and runs into the hover MarkSet, every buffer kept), on RB0b's stepped
// `heapScaling` (a step is one per-move update): a hover element closed and projected into every section's hover set
// with its primary and the overlay's queries; loci and items converted and projected; `Coordinates` through a
// registered translator; out-of-order runs sorted at close; row rectangles; a ligand's atoms (the cases in
// marks-cases.ts, shared with the quick twin, marks.quick.ts). The positive control: the sort's step allocating an
// object too (`leaking`), measured alike, is caught.

import { beforeAll, describe, expect, it } from 'vitest';
import { MarkSet } from '../src/marks.js';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { marksCases, sorting } from './marks-cases.js';

/** Steps in the shorter run; the longer runs four times as many. */
const N = 50_000;

let c: Awaited<ReturnType<typeof marksCases>>;

beforeAll(async () => {
  c = await marksCases();
}, 60_000);

describe('marks and the projector allocate nothing per move', () => {
  it("a hover: closed, projected into each section's hover set, its primary, the overlay's queries", () => {
    const result = heapScaling(N, c.hover());
    expect(c.sink()).not.toBe(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });

  it("loci: a residue's seq loci and the domain's items, converted and projected", () => {
    const result = heapScaling(N, c.loci());
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });

  it("Coordinates through a registered translator: a chain's residue onto its unit, then each block", () => {
    const result = heapScaling(N, c.coordinates());
    expect(c.sink()).not.toBe(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });

  it('a MarkSet: 200 runs out of order sorted and merged at close, then read by has', () => {
    const result = heapScaling(N / 10, sorting(new MarkSet()));
    expect(isFlat(result), describeScaling(result, N / 10)).toBe(true);
  });

  it("row rectangles from a row's residues, and a ligand's atoms at its one position", () => {
    const result = heapScaling(N, c.rowRects());
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  });

  it('the same sort allocating an object a step is caught (positive control)', () => {
    expect(controlCaught(() => heapScaling(N / 10, leaking(sorting(new MarkSet()))))).toBe(true);
  });
});
