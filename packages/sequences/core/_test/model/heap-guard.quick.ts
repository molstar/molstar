/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of heap-guard.heap.ts (`test:perf`): the self-test of the quick measurement (`quickCase`,
// `measureBytesPerStep`) on the same synthetic steps (heap-guard-cases.ts). A flat step reads about nothing and its
// positive control (`leaking`) is caught at `leak`'s 88 bytes; the two views a step read 208; an allocation at the first
// index only is behind the warm-up; garbage that collects many times in a window is accounted whole.

import { describe, expect, it } from 'vitest';
import { flatStep, garbage, onceStep, viewsStep } from './heap-guard-cases.js';
import { LEAK_BYTES, quickCase } from './quick.js';

describe("the quick measurement, on the guards' synthetic steps", () => {
  it('reads a step that allocates nothing as nothing, and catches its leak', async ({ task }) => {
    const { r, c } = quickCase(task, 'heap-guard: a step that allocates nothing', flatStep());
    expect(r.bytesPerStep).toBeLessThan(1);
    expect(c?.bytesPerStep ?? 0).toBeGreaterThan(0.75 * LEAK_BYTES);
    expect(c?.bytesPerStep ?? 0).toBeLessThan(1.25 * LEAK_BYTES);
  });

  it("reads BK1's case, two subarray views a step, at 208 bytes", async ({ task }) => {
    const { r, c } = quickCase(task, 'heap-guard: two subarray views a step', viewsStep);
    expect(r.bytesPerStep).toBeGreaterThan(190);
    expect(r.bytesPerStep).toBeLessThan(230);
    expect((c?.bytesPerStep ?? 0) - r.bytesPerStep).toBeGreaterThan(0.75 * LEAK_BYTES);
  });

  it('reads an allocation at the first index, behind the warm-up, as nothing', async ({ task }) => {
    const { r } = quickCase(task, 'heap-guard: 512 KB once', onceStep);
    expect(r.bytesPerStep).toBeLessThan(1);
  });

  it('accounts garbage that collects many times in a window', async ({ task }) => {
    // 40 KB a step: a window of 250 steps is 10 MB over the quick config's 1 MB young generation, collected many
    // times; the exact accounting still reads 40,000 bytes a step.
    const g = garbage(1000);
    const { r } = quickCase(task, 'heap-guard: 40 KB of garbage a step', g, {
      warmSteps: 100,
      steps: 250,
      windows: 3,
    });
    expect(r.collections).toBeGreaterThan(5);
    expect(Math.abs(r.bytesPerStep - 40_000)).toBeLessThan(2_000);
  });
});
