/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for translate's hot paths (TR1; seqspec:libraries#translate): `forEachAt` at pointer rate and `translate`
// in batches, through cached routes (one map, composed with a codon step either way, copies, complex members), on
// RB0b's stepped guard: runs of `n` and `4n` steps, the longer growing at most 64 KB beyond the shorter, and no GC in
// either. Untimed: `heapScaling` warms by count. The cases are translate-cases.ts's, shared with the quick twin
// (translate.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../model/heap.js';
import { resetVisits, TRANSLATE_CASES, translatedRuns, visited } from './translate-cases.js';

describe('translate heap guards', () => {
  for (const { name, n, step } of TRANSLATE_CASES) {
    it(name, () => {
      resetVisits();
      const r = heapScaling(n, step);
      expect(visited(), 'the steps visit images').toBeGreaterThan(n);
      expect(isFlat(r), describeScaling(r, n)).toBe(true);
      // The template's control (CB-05): the same step with a leak is caught.
      expect(controlCaught(() => heapScaling(n, leaking(step)))).toBe(true);
    });
  }

  it('translated what it should', () => {
    expect(translatedRuns()).toBe(3);
  });
});
