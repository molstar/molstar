/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A data hover's row band (the `hoverRow` behavior; seqspec:interaction#hover, seqspec:viewer#performance) allocates
// nothing per hover frame after warm-up: the pointer sweeping down the hover document's tracks (letters, ranges, points,
// values, a matrix's bands through its renderer's `rowHover`) at DPR 2, the element changing on most moves, an overlay
// frame per move on the recording backend. Counted over the measured frames alone, not vacuous: the element changed on
// most moves, and most frames drew the band (its header part at x 0 in the hover tint). The positive control: the same
// frames allocating an object each (`leaking`), through the same warm-up, settle and measurement, are caught. The
// viewer is loaded afresh per case (`vi.resetModules`), warmed and settled by count (`warmFrame`); the cases are
// row-hover-cases.ts's, shared with the quick twin (row-hover.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking, warmFrame } from './heap.js';
import { rowHover } from './row-hover-cases.js';

const N = 1_000;

/** The guard's measurement: warmed and settled by count (`warmFrame`), `warmed` called, then `N` frames and `4N`. */
function measure(frame: (i: number) => void, warmed = () => {}) {
  warmFrame(frame);
  warmed();
  return heapScaling(N, frame);
}

describe("a data hover's row band allocates nothing after warm-up", () => {
  it('a hover sweep down the tracks, the band drawn', async () => {
    const { frame, hovers, drawn } = await rowHover();
    let hovers0 = 0;
    let drawn0 = 0;
    let measured = 0;
    const r = measure(
      (i) => {
        frame(i);
        measured++;
      },
      () => {
        hovers0 = hovers();
        drawn0 = drawn();
        measured = 0;
      },
    );
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // Not vacuous, over the measured frames: the element changed on most moves, and most drew the band.
    expect(measured).toBeGreaterThanOrEqual(5 * N);
    expect(hovers() - hovers0).toBeGreaterThan(0.8 * measured);
    expect(drawn() - drawn0).toBeGreaterThan(0.8 * measured);
  }, 300_000);

  it('the same frames allocating an object each are caught (positive control)', async () => {
    const { frame } = await rowHover();
    expect(controlCaught(() => measure(leaking(frame)))).toBe(true);
  }, 300_000);
});
