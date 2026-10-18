/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of pairs-kinds.heap.ts (`test:perf`): a pan, a hover on a glyph of a cell of several kinds and the
// overlay, a frame (the guard's case, pairs-kinds-cases.ts), with its control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { kindsFrame } from './pairs-kinds-cases.js';

describe('frames over cells of several kinds, quick', () => {
  it('pans, a hover on each glyph and its overlay', async ({ task }) => {
    const { frame, multi, inSection } = await kindsFrame();
    expect(multi).toBe(true);
    expect(inSection).toBe(true);
    quickCase(task, 'pairs-kinds: pan, hover and overlay', frame, {
      warmSteps: 2_000,
      steps: 1_000,
      windows: 3,
      maxWarmSteps: 60_000,
    });
  });
});
