/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// EC-116 (.seq/sources/seqstar/edge-cases.md): a row of letters shows text for all its columns or none; as density
// grows it goes to cells, then dots, then bars. The steps are `densityStep`'s; the letters step never starts below the
// width the smallest letter needs, so there is no step with some cells lettered and some not.

import type { Node } from '@molstar/core/repr/tree';
import type { FormatRuntime } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { type DensityStep, DefaultRenderConfig, densityStep } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { BIG, BigRuntime, bigSpec, p53TracksSpec } from './fixtures.js';
import { blocks } from './helpers.js';

const DataWidth = 600;
const { inset, labelWidth, labelGap } = DefaultRenderConfig;
const Size = { width: DataWidth + inset + labelWidth + labelGap + inset, height: 120 };
const Order: DensityStep[] = ['letters', 'cells', 'dots', 'bars', 'density'];

/** The letters track's counts over `n` whole positions in view of a spec (a window at its middle). */
async function lettersOf(spec: Node, n: number, start: number, runtimes: readonly FormatRuntime[] = []) {
  const r = await renderSvg(spec, { size: Size, viewport: { start, end: start + n }, runtimes });
  return { step: densityStep(DataWidth / n, 1, DefaultRenderConfig), counts: blocks(r.lists)[0] };
}

describe('EC-116: letters, then cells, dots, bars as density grows', () => {
  it('the steps come in order, and the letters step starts where the smallest letter fits its cell', () => {
    let last = -1;
    let firstLetters = Infinity;
    for (let px = 40; px >= 1e-5; px *= 0.97) {
      const step = densityStep(px, 1, DefaultRenderConfig);
      const k = Order.indexOf(step);
      expect(k).toBeGreaterThanOrEqual(last);
      last = k;
      if (step === 'letters') firstLetters = Math.min(firstLetters, px);
    }
    expect(last).toBe(Order.length - 1);
    const { letterSizes, monoAdvance } = DefaultRenderConfig;
    expect(firstLetters).toBeGreaterThanOrEqual(letterSizes[0] * monoAdvance + 2);
  });

  it('text for all the columns in view or none: never a lettered few, from the widest to the cells', async () => {
    let lettered = 0;
    let plain = 0;
    // 20 to 120 positions in 600 px: 30 px a position down to 5
    for (const n of [20, 30, 40, 50, 60, 70, 80, 90, 100, 120]) {
      const { step, counts } = await lettersOf(p53TracksSpec(), n, 100);
      expect(counts.rects).toBe(n);
      if (step === 'letters') {
        expect(counts.letters).toBe(n);
        lettered++;
      } else {
        expect(counts.letters).toBe(0);
        plain++;
      }
    }
    // the sweep crosses the edge: both sides were seen
    expect(lettered).toBeGreaterThan(0);
    expect(plain).toBeGreaterThan(0);
  });

  it('then dots a position, then one bar a run: dots keep a mark a position, bars one for the letters', async () => {
    const dots = await lettersOf(p53TracksSpec(), 390, 1);
    expect(dots.step).toBe('dots');
    expect(dots.counts).toMatchObject({ letters: 0, texts: 0, rects: 390 });
    const mid = Math.floor(BIG / 2);
    const bars = await lettersOf(bigSpec(), 3000, mid, [BigRuntime]);
    expect(bars.step).toBe('bars');
    expect(bars.counts).toMatchObject({ letters: 0, rects: 1 });
  });
});
