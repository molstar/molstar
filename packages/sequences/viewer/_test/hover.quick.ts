/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of hover.heap.ts (`test:perf`): the hover sweeps over ranges, points, variants and a matrix, an
// overlay frame a move, each on a viewer loaded afresh (hover-sweep.ts, the guard's cases), each with its control.

import { describe, expect, it, vi } from 'vitest';
import { quickCase } from './heap.js';

describe('hover emphasis frames, quick', () => {
  it.for(['t.ranges', 't.points', 't.variants', 't.matrix'] as const)(
    'a hover sweep over %s',
    async (track, { task }) => {
      vi.resetModules();
      const { hoverSweep } = await import('./hover-sweep.js');
      const { frame, hovers, drawn } = await hoverSweep(track);
      quickCase(task, `hover: ${track}`, frame);
      expect(hovers()).toBeGreaterThan(1_000);
      expect(drawn()).toBeGreaterThan(1_000);
    },
  );
});
