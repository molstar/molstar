/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LM1's letter maps through the conformance kit's heap guard (seqspec:viewer#extension-points;
// seqspec:viewer#the-render-backend): the map lookup is per drawn letter and per cell, so `letters` (its own, `of`
// another polymer's, and 20,000 positions into classes) and `rows` (their tiles and letters) with a complement map, and
// `emphasis` redrawing its track's mapped letters, forward and as a reverse complement, linear and wrapped, allocate
// nothing after warm-up. Each frame the kit draws is warmed by count until it stops allocating (`warmFrame`), then
// measured over `FRAME_N` steps and `4 * FRAME_N`, again when not flat. The positive control goes through the kit's
// `heap` hook on letters over every type: at each frame the same warm-up and measurement of the frame allocating an
// object each (`leaking`) must catch it. Each case loads the viewer afresh (`vi.resetModules`; draws.heap.ts says why);
// the cases are letter-maps-cases.ts's, shared with the quick twin (letter-maps.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, FRAME_N, HEAP_DEEP, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { CASES, REVERSED, run } from './letter-maps-cases.js';

const N = FRAME_N;
/** The guards' measurement: warmed by count, then `N` frames and `4N` measured, again when not flat. */
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

describe('LM1: frames drawn through a letter map allocate nothing after warm-up', () => {
  for (const reversed of [false, true])
    for (const [name, doc] of CASES)
      if (!reversed || HEAP_DEEP || REVERSED.has(`${name} on ${doc}`))
        it(`${name} on ${doc}, ${reversed ? 'as a reverse complement' : 'complemented'}`, async () => {
          const heap = (await run(name, doc, reversed, guard)).filter((c) => c.check === 'heap');
          expect(heap.length).toBeGreaterThan(0);
          expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
        }, 240_000);

  it('letters on all types, an object allocated a frame (the positive control)', async () => {
    const heap = (await run('letters', 'all types', false, control)).filter((c) => c.check === 'heap');
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  }, 240_000);
});
