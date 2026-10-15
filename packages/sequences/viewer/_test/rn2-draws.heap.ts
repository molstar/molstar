/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// RN2's feature, pair-overlay and logo draws through the conformance kit's heap guard (seqspec:viewer#extension-points;
// seqspec:viewer#the-render-backend): a frame of `ranges`, `points`, `variants`, `secondary` and `logo` allocates
// nothing after warm-up, at each of the kit's zooms and below the row's minimum, linear and wrapped, on documents that
// reach their paths (every feature layer; boundaries pinned over a track's letters; thousands of items through their
// density, raw and from the count pyramid; transcripts with hairlines; 8,000 alignment columns for the logo's thin path,
// raw and from its pyramid); the `pairs` overlay in both of its modes over many lit ends. The positive control goes
// through the kit's `heap` hook on ranges over the features document: at each frame the same warm-up and measurement of
// the frame allocating an object each (`leaking`) must catch it. Each case loads the viewer afresh (`vi.resetModules`;
// draws.heap.ts says why: Maglev and the module runner). The cases are rn2-draws-cases.ts's, shared with the quick twin
// (rn2-draws.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { CASES, kitRun, PAIRS_MODES, pairsOverlay } from './rn2-draws-cases.js';

const N = FRAME_N;
/**
 * The guards' measurement: warmed by count (`warmFrame`: `FRAME_WARM` frames, then windows until the frame stops
 * allocating, at most `FRAME_WARM_CAP`; a tiny draw climbs to TurboFan over tens of thousands), then `N` frames and
 * `4N` measured, again when not flat.
 */
const measure = (step: (i: number) => void) => {
  let i = 0;
  const frame = () => step(i++);
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

describe("RN2's feature and logo frames allocate nothing after warm-up", () => {
  for (const c of CASES)
    it(`${c.entry} on ${c.doc}`, async () => {
      const heap = (await kitRun(c, guard)).filter((x) => x.check === 'heap');
      expect(heap.length).toBeGreaterThan(0);
      expect(heap.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
    }, 240_000);

  it('ranges on features, an object allocated a frame (the positive control)', async () => {
    const heap = (await kitRun({ entry: 'ranges', doc: 'features' }, control)).filter((x) => x.check === 'heap');
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((x) => !x.ok).map((x) => x.detail)).toEqual([]);
  }, 240_000);
});

describe('the pairs overlay allocates nothing after warm-up', () => {
  it.each(PAIRS_MODES)(
    'arcs %s, many lit ends',
    async (mode) => {
      const { frame, paths } = await pairsOverlay(mode);
      expect(paths()).toBeGreaterThan(0);
      guard(frame);
    },
    240_000,
  );
});
