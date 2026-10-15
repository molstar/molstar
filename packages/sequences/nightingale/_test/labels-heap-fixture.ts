/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The label layer's cases (NG-12), for its heap guards and their quick twins (labels-lay.heap.ts, labels-place.heap.ts,
// labels-lay.quick.ts, labels-place.quick.ts): a few hundred labels, widths and heights fractional as production's
// (`measureText().width`), so a double crossing a call that is not inlined boxes and the guards see it (LA-13), and the
// two steps, `layLabels` at an unchanged view and `placeLabel` at a panned one.

import { axisView, type LabelBox, type Laid, layLabels, type NightingaleAxis, placeLabel } from '../src/labels.js';

export const N = 50_000;
export const COUNT = 300;
/** A label's height, fractional. */
export const HEIGHT = 14.5;
const style = () => ({ display: '', left: '', width: '', top: '', height: '', lineHeight: '' });
export const labels: Laid[] = Array.from({ length: COUNT }, (_, k) => ({
  feature: { k },
  text: `feature ${k}`,
  fragments: [
    { start: 10 * k + 1, end: 10 * k + 4 },
    { start: 10 * k + 6, end: 10 * k + 9 },
  ],
  style: style(),
  width: -1,
  cut: -1,
  box: { mode: 'none', left: 0, width: 0 },
  top: 0,
  height: 0,
  visible: false,
}));
export const layout = { getFeatureYPos: () => 4.25, getFeatureHeight: () => HEIGHT };
export const measure = (text: string) => 6.37 * text.length;
export const view: NightingaleAxis = axisView({ start: 0, end: 0, left: 0, right: 0 }, 1, -1, 10 * COUNT, 3000, 10, 10);

/** `layLabels` at an unchanged view: laid, compared, nothing written. `writes` counts what it wrote (none). */
export function laying(): { step: () => void; writes: () => number } {
  layLabels(labels, view, layout, measure);
  let writes = 0;
  return {
    step: () => {
      writes += layLabels(labels, view, layout, measure);
    },
    writes: () => writes,
  };
}

/**
 * `placeLabel` at a panned view each step, zoomed in and out every other step: a label's 4-position fragment is about
 * 200 px wide over 60 positions (its text fits: "full") and 20 px over 600 (its cut or nothing). `modes` counts every
 * mode placed.
 */
export function placing(): { step: (i: number) => void; modes: { full: number; cut: number; none: number } } {
  // Laid once, for the labels' measured widths and cuts.
  layLabels(labels, view, layout, measure);
  const [panned, out]: [NightingaleAxis, LabelBox] = [{ ...view }, { mode: 'none', left: 0, width: 0 }];
  const modes = { full: 0, cut: 0, none: 0 };
  const step = (i: number) => {
    panned.start = 1 + (i % 500);
    panned.end = panned.start + ((i & 1) === 0 ? 60 : 600);
    const l = labels[i % COUNT];
    modes[placeLabel(out, panned, l.fragments, l.width, l.cut, HEIGHT).mode]++;
  };
  return { step, modes };
}
