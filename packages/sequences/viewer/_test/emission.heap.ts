/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for emission (seqspec:viewer#performance): a step a set, a bump or an emission, over runs of n and 4n
// (`heapScaling`, RB0b), each path warmed up first (the measurement's own warm-up of n steps). The steps are
// emission-cases.ts's, shared with the quick twin (emission.quick.ts). The positive control: the state sets allocating
// an object a step (`leaking`), measured alike, are caught.

import { describe, expect, it } from 'vitest';
import { bumping, emitting, N, queued, setting, sunk } from './emission-cases.js';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from './heap.js';

describe('emission allocates nothing', () => {
  it('state sets, with rxjs subscribers', () => {
    const r = heapScaling(N, setting());
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('state sets queued during delivery', () => {
    const r = heapScaling(N, queued());
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('live version bumps', () => {
    const r = heapScaling(N, bumping());
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('events with a reused record, by priority', () => {
    const r = heapScaling(N, emitting());
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(sunk()).toBeGreaterThanOrEqual(0);
  });

  it('the same sets allocating an object each are caught (positive control)', () => {
    expect(controlCaught(() => heapScaling(N, leaking(setting())))).toBe(true);
  });
});
