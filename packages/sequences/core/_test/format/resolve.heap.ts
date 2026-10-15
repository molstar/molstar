/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for resolve's per-row and per-position paths (VW4), on RB0b's stepped `heapScaling` (a step is one row;
// untimed, `heapScaling` warms by count): `cover` (what `positions` and an annotation focus mark) through a codon map,
// the filter's test of a row's value (a number column compared as numbers, no double boxed), and a row's color through
// the layer's one categorical palette. Resolve allocates per layer and per call, never per row or position. The cases
// are resolve-cases.ts's, shared with the quick twin (resolve.quick.ts).

import { beforeAll, describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../model/heap.js';
import { type ResolveCases, resolveCases } from './resolve-cases.js';

/** Steps in the shorter run; the longer runs four times as many (the cases' `ROWS`). */
const N = 20_000;
let cases: ResolveCases;

beforeAll(async () => {
  cases = await resolveCases();
}, 60_000);

describe("resolve's per-row paths allocate nothing", () => {
  it("cover: a row's codons marked through the map", () => {
    expect(cases.path.span).toBe(3);
    const r = heapScaling(N, cases.cover);
    expect(cases.mask.indexOf(0)).toBe(-1);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // The template's control (CF-05): the same step with a leak is caught.
    expect(controlCaught(() => heapScaling(N, leaking(cases.cover)))).toBe(true);
  });

  it("the filter's test of a row, a string and a number field", () => {
    const r = heapScaling(N, cases.filter);
    expect(cases.kept()).toBeGreaterThan(0);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it("a row's color through the layer's categorical palette", () => {
    const before = cases.colored();
    const r = heapScaling(N, cases.color);
    expect((cases.colored() - before) % N).toBe(0); // every row colored, in every run
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });
});
