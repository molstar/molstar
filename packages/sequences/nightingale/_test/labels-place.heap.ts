/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guard for the label layer (NG-12), on RB0b's stepped `heapScaling`: `placeLabel` at a panned view each step
// (every zoom frame; the style writes, strings by nature, aside), its width, cut and height fractional (LA-13). A file
// of its own (LA-28). The positive control is the template's: the same step with a leak.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../../core/_test/model/heap.js';
import { type LabelBox, layLabels, type NightingaleAxis, placeLabel } from '../src/labels.js';
import { COUNT, HEIGHT, labels, layout, measure, N, view } from './labels-heap-fixture.js';

describe('the label layer allocates nothing per frame', () => {
  it('placeLabel at a panned view each step', () => {
    // Laid once, for the labels' measured widths and cuts.
    layLabels(labels, view, layout, measure);
    const [panned, out]: [NightingaleAxis, LabelBox] = [{ ...view }, { mode: 'none', left: 0, width: 0 }];
    // Zoomed in and out every other step: a label's 4-position fragment is about 200 px wide over 60 positions (its
    // text fits: "full") and 20 px over 600 (its cut or nothing). Every mode is placed.
    const modes = { full: 0, cut: 0, none: 0 };
    const step = (i: number) => {
      panned.start = 1 + (i % 500);
      panned.end = panned.start + ((i & 1) === 0 ? 60 : 600);
      const l = labels[i % COUNT];
      modes[placeLabel(out, panned, l.fragments, l.width, l.cut, HEIGHT).mode]++;
    };
    const r = heapScaling(N, step);
    expect([modes.full > 0, modes.cut > 0, modes.none > 0]).toEqual([true, true, true]);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });
});
