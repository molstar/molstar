/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { referenceExtent } from './values-reference.js';
import type { ValuesLay } from './values.js';

// The values layer's open domain (_spec/README.md "Render", Layer types, `values`): over the values in view and its
// reference lines.

/** The values in view, min to max; `[0, 1]` when there are none, widened around a single value. */
export function openDomain(lay: ValuesLay): [number, number] {
  // the lines are part of what an open end is computed from (a threshold's absence from a plot is the finding)
  const seed = referenceExtent(lay);
  let lo = seed[0];
  let hi = seed[1];
  let pos = seed[2];
  for (let j = 0; j < lay.n; j++) {
    const v = lay.v[j];
    if (Number.isNaN(v)) continue;
    if (lay.lo[j] < lo) lo = lay.lo[j];
    if (v > hi) hi = v;
    const p = lay.lo[j] > 0 ? lay.lo[j] : v;
    if (p > 0 && p < pos) pos = p;
  }
  // a log scale's open low end: the smallest positive value in view (a tenth of its top when there is one alone)
  if (lay.yLog) return hi > 0 ? [pos < hi ? pos : hi / 10, hi] : [1, 10];
  if (lo > hi) return [0, 1];
  if (lo === hi) return lo >= 0 ? [0, hi || 1] : [lo, 0];
  return [Math.min(lo, 0), hi];
}
