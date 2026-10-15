/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the `follow` behavior (LK1; seqspec:viewer#performance), on `heapScaling`: a leader with two
// followers (a protein through the codon step and the gene reversed), driven at input rate by a drag and by a Ctrl
// wheel zoom; each step moves all three. Fractional coordinates, as an engine's; every default behavior subscribed. The
// steps are follow-cases.ts's, shared with the quick twin (follow.quick.ts); each is warmed and settled by count
// (`frameScaling`), then `N` and `4N` steps are measured. The positive control: the wheel's steps allocating an object
// each (`leaking`), warmed and measured alike, are caught.

import { describe, expect, it } from 'vitest';
import { dragging, N, wheeling } from './follow-cases.js';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';

describe('follow at input rate', () => {
  it('a drag on the leader moves both followers each step: flat', async () => {
    const { step, protein, rev, finish } = await dragging();
    const [p0, r0] = [protein.live.version.value, rev.live.version.value];
    const r = frameScaling(step, N);
    finish();
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(protein.live.version.value - p0).toBeGreaterThan(N);
    expect(rev.live.version.value - r0).toBeGreaterThan(N);
    expect(protein.live.current.origin).toBe('behavior');
  }, 300_000);

  it('a Ctrl wheel zoom on the leader, in and out, zooms both followers each step: flat', async () => {
    const { step, dna, protein, rev } = await wheeling();
    const [p0, r0] = [protein.live.version.value, rev.live.version.value];
    const r = frameScaling(step, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(protein.live.version.value - p0).toBeGreaterThan(N);
    expect(rev.live.version.value - r0).toBeGreaterThan(N);
    expect(rev.live.current.cell).toBe(dna.live.current.cell);
  }, 300_000);

  it('the same wheel steps allocating an object each are caught (positive control)', async () => {
    const { step } = await wheeling();
    expect(controlCaught(() => frameScaling(leaking(step), N))).toBe(true);
  }, 300_000);
});
