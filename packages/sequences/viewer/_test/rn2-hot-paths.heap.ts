/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN2's hot paths allocate nothing per frame after warm-up (RN1's pattern): the shared density (points stacked and
// majority, raw and from the count pyramid; ranges' coverage; a block counted through its items beside one through its
// pyramid) through the test entries `ranges`, `points` and `variants` will become; `pairs`' draw and `tiny` through the
// kit at each of its zooms, and its hit over a grid of points; the matrix's tiles at the kit's zooms. The positive
// control goes through the kit's `heap` hook on raw points: at each frame the same warm-up and measurement of the frame
// allocating an object each (`leaking`) must catch it. Each case loads the viewer afresh (`vi.resetModules`;
// draws.heap.ts says why: Maglev and the module runner). The cases are rn2-hot-paths-cases.ts's, shared with the quick
// twin (rn2-hot-paths.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { CASES, caseName, kitRun, litPairsFrame, pairsHit } from './rn2-hot-paths-cases.js';

const N = FRAME_N;
/**
 * The guards' measurement: warmed by count (`warmFrame`: `FRAME_WARM` frames, then windows until the frame stops
 * allocating, at most `FRAME_WARM_CAP`), then `N` frames and `4N`, measured twice when the first is not flat: the
 * first run after the kit switches contexts can grow once by about 3 MB (measured: the short run only, the long run
 * and every later run flat; none in one context), which no frame repeats. A frame that allocates fails both.
 */
const measure = (frame: (i: number) => void) => {
  warmFrame(frame);
  const r = heapScaling(N, frame);
  return isFlat(r) ? r : heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  if (!isFlat(r)) throw new Error(describeScaling(r, N));
};
/** The positive control as the kit's `heap` hook: a frame's check fails unless its leaking frame is caught. */
const control = (frame: () => void) => {
  if (!controlCaught(() => measure(leaking(frame)))) throw new Error('the leaking frame was not caught');
};

describe("RN2's hot paths allocate nothing after warm-up", () => {
  for (const c of CASES)
    it(caseName(c), async () => {
      const heap = (await kitRun(c, guard)).filter((x) => x.check === 'heap');
      expect(heap.length).toBeGreaterThan(0);
      expect(heap.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
    }, 240_000);

  it('points, raw, an object allocated a frame (the positive control)', async () => {
    const heap = (await kitRun({ entry: 'points', doc: 'points, raw' }, control)).filter((x) => x.check === 'heap');
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
  }, 240_000);

  it('pairs: the hit over a grid of fractional points beside one long pair', async () => {
    guard(await pairsHit());
  }, 240_000);

  it("pairs: the overlay's lit pairs through the channels' runs and the end index", async () => {
    guard(await litPairsFrame());
  }, 240_000);
});
