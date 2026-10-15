/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guard for LS1's per-hover path in `mount`: a section's `highlights` (`AxisMeta.highlights`) kept under a host's
// highlight (`keepLit`) on every change, after the host's loci are projected, on RB0b's stepped `heapScaling`; the
// step is the case of kept-cases.ts, shared with the quick twin (kept.quick.ts). The positive control is the
// template's: the same step with a leak, at the guard's own length (NG-11: the strings it kept on a list could cancel
// in the shorter run).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../../core/_test/model/heap.js';
import { kept, N } from './kept-cases.js';

describe("a section's highlights kept under a host's", () => {
  it('allocate nothing per change', async () => {
    const c = await kept();
    const r = heapScaling(N, c.step);
    expect(c.lit()).toBeGreaterThan(40 * N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
  });

  it('the same step with a leak is caught (the positive control)', async () => {
    const c = await kept();
    expect(controlCaught(() => heapScaling(N, leaking(c.step)))).toBe(true);
  });
});
