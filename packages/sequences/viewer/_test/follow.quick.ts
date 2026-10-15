/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of follow.heap.ts (`test:perf`): a drag and a Ctrl wheel zoom on the leader, each moving both
// followers a step, the guard's steps (follow-cases.ts), each with its positive control.

import { describe, expect, it } from 'vitest';
import { dragging, wheeling } from './follow-cases.js';
import { quickCase } from './heap.js';

describe('follow at input rate, quick', () => {
  it('a drag on the leader moves both followers each step', async ({ task }) => {
    const { step, protein, rev, finish } = await dragging();
    const [p0, r0] = [protein.live.version.value, rev.live.version.value];
    // The first case of the file warms the whole input pipeline: its tier-up takes about 10,000 steps.
    quickCase(task, 'follow: a drag on the leader', step, {
      warmSteps: 2_000,
      maxWarmSteps: 40_000,
    });
    finish();
    expect(protein.live.version.value - p0).toBeGreaterThan(1_000);
    expect(rev.live.version.value - r0).toBeGreaterThan(1_000);
  });

  it('a Ctrl wheel zoom on the leader zooms both followers each step', async ({ task }) => {
    const { step, dna, protein, rev } = await wheeling();
    const [p0, r0] = [protein.live.version.value, rev.live.version.value];
    quickCase(task, 'follow: a Ctrl wheel zoom on the leader', step);
    expect(protein.live.version.value - p0).toBeGreaterThan(1_000);
    expect(rev.live.version.value - r0).toBeGreaterThan(1_000);
    expect(rev.live.current.cell).toBe(dna.live.current.cell);
  });
});
