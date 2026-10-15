/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of the font metrics measurer's guards (metrics.heap.ts) and their quick twin (metrics.quick.ts): `metrics`
// over two font strings, `measure` over strings of ASCII and of other code points. The heap guard and the quick twin
// share them.

import { DefaultFontMetrics, type EmBox, metricsMeasurer } from '../../src/svg/index.js';

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 100_000;
const m = metricsMeasurer(DefaultFontMetrics);
const fonts = ['400 12px ui-sans-serif, system-ui, sans-serif', '600 11px ui-monospace, Menlo, monospace'];
const texts = ['chr12:1,234,567-1,299,999', 'Zinc ligands (UniProt) é', 'ABCDEFGHIKLMNPQRSTVWY'];

/** `metrics` for a font string already known; `sum` reads the ascents added. */
export function metrics(): { step: (i: number) => void; sum: () => number } {
  const sum = new Float64Array(1);
  return {
    step: (i) => {
      sum[0] += (m.metrics(fonts[i & 1]) as EmBox).ascent;
    },
    sum: () => sum[0],
  };
}

/** `measure`; `sum` reads the widths added (a captured `let` double would be boxed on every step, so a typed array). */
export function measuring(): { step: (i: number) => void; sum: () => number } {
  const sum = new Float64Array(1);
  return {
    step: (i) => {
      sum[0] += m.measure(fonts[i & 1], texts[i % 3]);
    },
    sum: () => sum[0],
  };
}
