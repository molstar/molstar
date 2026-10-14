/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The label layer's heap guards' fixture (NG-12): a few hundred labels, widths and heights fractional as production's
// (`measureText().width`), so a double crossing a call that is not inlined boxes and the guards see it (LA-13).

import { axisView, type Laid, type NightingaleAxis } from '../src/labels.js';

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
