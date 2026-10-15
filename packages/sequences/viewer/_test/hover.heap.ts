/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The hover emphasis allocates nothing per hover frame after warm-up (HV1; seqspec:interaction#hover,
// seqspec:viewer#performance): the pointer sweeping over ranges, points, variants and a matrix, alternating between
// their items so the hovered element changes on most moves, an overlay frame per move on the recording backend (the
// channels' hover rebuilt, the item's features found through the per-state index, its outline in the cached hover ink).
// Not vacuous: the hovered element changed on most moves, and most frames drew an emphasis. The positive control: the
// same frames allocating an object each (`leaking`), through the same warm-up, settle and measurement, are caught.
//
// Each case loads the viewer afresh (`vi.resetModules`), as overlay.heap.ts does, and is warmed and settled by count
// (`frameScaling`): the composer's passes are calls across its modules, a longer climb to TurboFan. The cases are
// hover-sweep.ts's, shared with the quick twin (hover.quick.ts).

import { describe, expect, it, vi } from 'vitest';
import { controlCaught, describeScaling, FRAME_WARM, frameScaling, isFlat, leaking } from './heap.js';

const N = 1_000;

/** The guard's measurement: warmed and settled by count, then `N` frames and `4N` measured. */
const measure = (frame: (i: number) => void) => frameScaling(frame, N);

describe('hover emphasis frames allocate nothing after warm-up', () => {
  it.each(['t.ranges', 't.points', 't.variants', 't.matrix'] as const)(
    'a hover sweep over %s',
    async (track) => {
      vi.resetModules();
      const { hoverSweep } = await import('./hover-sweep.js');
      const { frame, hovers, drawn } = await hoverSweep(track);
      const r = measure(frame);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      // Not vacuous: the element changed on most moves, and most frames outlined it.
      expect(hovers()).toBeGreaterThan(FRAME_WARM);
      expect(drawn()).toBeGreaterThan(FRAME_WARM);
    },
    300_000,
  );

  it('the same frames allocating an object each are caught (positive control)', async () => {
    vi.resetModules();
    const { hoverSweep } = await import('./hover-sweep.js');
    const { frame } = await hoverSweep('t.ranges');
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
  }, 300_000);
});
