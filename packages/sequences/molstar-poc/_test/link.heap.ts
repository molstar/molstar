/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The link's 3D-hover path (seqspec:libraries#the-mol-adapter, Rates): a Mol* loci to runs in the member systems and on
// to Seq* units, into reused `Extents`, allocates nothing, on RB0b's stepped guard; the step is link-cases.ts's, shared
// with the quick twin (link.quick.ts). The positive control: the same step allocating an object each (`leaking`),
// measured alike, is caught. The Seq*-to-3D hover makes a Mol* loci per hover change by design (the link's `react`),
// so it has no guard.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../../core/_test/model/heap.js';
import { molHover, N } from './link-cases.js';

describe('link heap guards', () => {
  it('a Mol* hover to member runs and Seq* runs', async () => {
    const { step, sink } = await molHover();
    const r = heapScaling(N, step);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(sink()).not.toBe(0.5);
  });

  it('the same step allocating an object each is caught (positive control)', async () => {
    const { step } = await molHover();
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });
});
