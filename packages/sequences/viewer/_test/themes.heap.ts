/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for LY3a's per-element path (seqspec:viewer#performance; seqspec:viewer#extension-points): `theme.value`
// (and `alpha`) per position for every built-in and every named rule table compiled into `rules`, on `heapScaling` (a
// step one element; it warms by count, n steps before each measured run). The positive control: `letter`'s step
// allocating an object too (`leaking`), measured alike, is caught. The cases are themes-cases.ts's, shared with the
// quick twin (themes.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';
import { CASES, reading } from './themes-cases.js';

const N = 50_000;

describe('theme values per element allocate nothing', () => {
  for (const { name, theme } of CASES) {
    it(name, () => {
      expect(theme).not.toBeNull();
      const r = heapScaling(N, reading(theme));
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
    });
  }

  it('the same elements allocating an object each are caught (positive control)', () => {
    expect(controlCaught(() => heapScaling(N, leaking(reading(CASES[0].theme))))).toBe(true);
  });
});
