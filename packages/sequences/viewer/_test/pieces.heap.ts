/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// AX1's reversed draws through the conformance kit's heap guard (seqspec:viewer#extension-points;
// seqspec:viewer#the-render-backend): every built-in renderer's frame on reversed and ranged pieces (every section
// `direction: -1`, every block of a sequence a `range`) allocates nothing after warm-up, at each of the kit's zooms and
// below the row's minimum, linear and wrapped, on documents that reach their paths: every type on two chains, 40,000
// positions zoomed out (classes, pyramids), dense means, thousands of features through their density, an alignment's
// rows and logo, a matrix and a miniature. The run conversions are the axis's (`runToShown`, `shownToRuns`),
// allocation-free in either direction. Then the overlay's hover frames (overlay.heap.ts's) on reversed pieces. Each
// frame the kit draws is warmed by count until it stops allocating (`warmFrame`), then measured over `FRAME_N` steps
// and `4 * FRAME_N`, again when not flat. The positive control goes through the kit's `heap` hook on reversed letters
// over every type: at each frame the same warm-up and measurement of the frame allocating an object each (`leaking`)
// must catch it. Each case loads the viewer afresh (`vi.resetModules`; draws.heap.ts says why: Maglev and the module
// runner); the cases are pieces-cases.ts's, shared with the quick twin (pieces.quick.ts).

import { describe, expect, it } from 'vitest';
import {
  cheapScaling,
  controlCaught,
  describeScaling,
  FRAME_N,
  FRAME_WARM,
  HEAP_DEEP,
  heapScaling,
  isFlat,
  leaking,
  warmFrame,
} from './heap.js';
import { CASES, DEEP_CASES, reversedHover, run } from './pieces-cases.js';

const N = FRAME_N;

/** The kit's guards' measurement: warmed by count, then `N` frames and `4N` measured, again when not flat. */
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

describe("AX1: reversed and ranged pieces' frames allocate nothing after warm-up", () => {
  for (const [name, doc] of HEAP_DEEP ? [...CASES, ...DEEP_CASES] : CASES)
    it(`${name} on ${doc}, reversed`, async () => {
      const heap = (await run(name, doc, guard)).filter((c) => c.check === 'heap');
      expect(heap.length).toBeGreaterThan(0);
      expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
    }, 240_000);

  it('letters on all types, reversed, an object allocated a frame (the positive control)', async () => {
    const heap = (await run('letters', 'all types', control)).filter((c) => c.check === 'heap');
    expect(heap.length).toBeGreaterThan(0);
    expect(heap.filter((c) => !c.ok).map((c) => c.detail)).toEqual([]);
  }, 240_000);
});

describe("AX1: the overlay's hover frames on reversed pieces allocate nothing after warm-up", () => {
  it.each(['scale', 'rows'] as const)(
    'hover frames over %s, reversed, with a selection and a highlight',
    async (kind) => {
      const { frame, hovers } = await reversedHover(kind);
      const r = cheapScaling(frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(hovers()).toBeGreaterThan(FRAME_WARM);
    },
    300_000,
  );
});
