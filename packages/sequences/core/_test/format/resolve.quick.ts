/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of resolve.heap.ts (`test:perf`): resolve's per-row paths (resolve-cases.ts, the guard's cases), each
// with its positive control.

import { beforeAll, describe, expect, it } from 'vitest';
import { quickCase } from '../model/quick.js';
import { type ResolveCases, resolveCases } from './resolve-cases.js';

let cases: ResolveCases;

beforeAll(async () => {
  cases = await resolveCases();
}, 60_000);

describe("resolve's per-row paths, quick", () => {
  it("cover: a row's codons marked through the map", async ({ task }) => {
    quickCase(task, 'resolve: cover', cases.cover);
    expect(cases.path.span).toBe(3);
    expect(cases.mask.indexOf(0)).toBe(-1);
  });

  it("the filter's test of a row, a string and a number field", async ({ task }) => {
    quickCase(task, 'resolve: filter', cases.filter);
    expect(cases.kept()).toBeGreaterThan(0);
  });

  it("a row's color through the layer's categorical palette", async ({ task }) => {
    quickCase(task, 'resolve: categorical color', cases.color);
    expect(cases.colored()).toBeGreaterThan(0);
  });
});
