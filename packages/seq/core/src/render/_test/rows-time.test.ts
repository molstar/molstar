/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows' targets (the design's "Performance and size"): over a 10,000 × 300 alignment, a first frame (the join, the
// measures, the order, the theme's classes, the reduction, the layout) under 300 ms; a pan frame under 8 ms at the fit
// (image tiles), at letters and between two level pairs; a theme change's rebuild under 50 ms; folding every clade
// under 16 ms; and a frame of 10,000,000 lazy rows at letters within 1.2× of the held one. Each bound is the heavy
// lane's (`itWithin`, best of three); the quick lane runs the same frames unbounded over 1,000 rows.

import { BuiltInLayerTypes, BuiltInLayouts } from '../built-ins.js';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { DarkTokens, type FrameInput, LightTokens, Renderer, type RowsLay } from '../index.js';
import { settle } from './helpers.js';
import { Columns, rowsSpec as spec, Section, Synthetic, Track } from './rows-fixtures.js';

// the bounds are the design's, at 10,000 rows; the quick lane runs the same bodies over 1,000: the 10,000-row joins and
// reductions are slow by nature, and the heavy lane (and `all`) keeps them
const Rows = ['heavy', 'all'].includes((import.meta as { env?: { MODE?: string } }).env?.MODE ?? '') ? 10_000 : 1_000;

const Size = { width: 1440, height: 800 };
const at = (bandTop: number, rowScale: number, start = 0, end = Columns, folds?: readonly string[]): FrameInput => ({
  size: Size,
  dpr: 2,
  settled: true,
  viewports: { [Section]: { start, end, bands: { [Track]: { bandTop, rowScale, ...(folds ? { folds } : {}) } } } },
});

function rowsLay(f: {
  layout: { sections: readonly { lays: readonly (readonly { lay: unknown; type: { name: string } }[])[] }[] };
}) {
  for (const s of f.layout.sections)
    for (const b of s.lays) for (const l of b ?? []) if (l.type.name === 'rows') return l.lay as RowsLay;
  throw new Error('no rows');
}

/** The median of `n` frames, each its own input (a pan: every frame lays out again). */
function median(renderer: Renderer, inputs: readonly FrameInput[], n = 30): number {
  const t: number[] = [];
  for (let k = 0; k < n; k++) {
    const t0 = performance.now();
    const f = renderer.frame(inputs[k % inputs.length]);
    t.push(performance.now() - t0);
    if (!f.rebuilt) throw new Error('a pan frame lays out again');
  }
  t.sort((a, b) => a - b);
  return t[t.length >> 1];
}

describe('the rows’ time over 10,000 × 300', () => {
  itWithin(
    'a first frame: the join, the measures, the order, the classes, the reduction, the layout',
    async (within) => {
      const s = await settle(spec(Rows, true), [Synthetic]);
      await s.frame(at(0, 0.06));
      const first = () => {
        const r = new Renderer({
          types: BuiltInLayerTypes,
          layouts: BuiltInLayouts,
          model: s.model,
          view: s.view,
          data: s.store,
          tokens: LightTokens,
          text: s.renderer.options.text,
        });
        const t0 = performance.now();
        r.frame(at(0, 0.06));
        return performance.now() - t0;
      };
      const ms = within(300, first);
      console.info(`rows: first frame ${ms.toFixed(1)} ms`);
      expect(ms).toBeGreaterThan(0);
    },
  );

  itWithin('a pan frame at the fit (image tiles), at letters and between two level pairs', async (within) => {
    const s = await settle(spec(Rows, true), [Synthetic]);
    const fit = [0, 1, 2, 3].map((k) => at(k, 0.06));
    const letters = [0, 20, 40, 60].map((k) => at(4000 + k / 20, 20, 100, 160));
    const pairs = [0, 4, 8, 12].map((k) => at(k * 20, 0.5, 0, Columns));
    for (const inputs of [fit, letters, pairs]) for (const i of inputs) await s.frame(i);
    const f = await s.frame(fit[0]);
    expect(rowsLay(f).images.length).toBeGreaterThan(0);
    const l = rowsLay(await s.frame(letters[0]));
    expect(l.letters).toBe(true);
    for (const [name, inputs] of [
      ['fit', fit],
      ['letters', letters],
      ['level pairs', pairs],
    ] as const) {
      median(s.renderer, inputs);
      const ms = within(8 * 30, () => median(s.renderer, inputs));
      console.info(`rows: pan frame at ${name}: median ${ms.toFixed(3)} ms`);
    }
  });

  itWithin('a theme change rebuilds the classes and the reduction under 50 ms', async (within) => {
    const s = await settle(spec(Rows, false), [Synthetic]);
    await s.frame(at(0, 0.06));
    let dark = false;
    const change = () => {
      dark = !dark;
      s.renderer.tokens = dark ? DarkTokens : LightTokens;
      const t0 = performance.now();
      s.renderer.frame(at(0, 0.06));
      return performance.now() - t0;
    };
    change();
    const ms = within(50, change);
    console.info(`rows: theme change ${ms.toFixed(1)} ms`);
  });

  itWithin('folding every clade under 16 ms', async (within) => {
    const s = await settle(spec(Rows, true), [Synthetic]);
    const f = await s.frame(at(0, 20));
    const l = rowsLay(f);
    const names = l.clades!.name.filter((_, c) => c > 0 && l.clades!.hi[c] - l.clades!.lo[c] > 1);
    let k = 0;
    const fold = () => {
      // a new set each time (the folds compare by identity): every clade, then none
      const folds = k++ % 2 ? [] : [...names];
      const t0 = performance.now();
      s.renderer.frame(at(0, 20, 0, 60, folds));
      return performance.now() - t0;
    };
    fold();
    fold();
    const ms = within(16, fold);
    console.info(`rows: fold every clade ${ms.toFixed(2)} ms`);
    expect(rowsLay(s.renderer.frame(at(0, 20, 0, 60, [...names]))).display.count).toBeLessThan(10);
  });

  it('a frame of 10,000,000 lazy rows at letters within 1.2× of 10,000 held', { tags: ['heavy'] }, async () => {
    const held = await settle(spec(10_000, false), [Synthetic]);
    const lazy = await settle(spec(10_000_000, false), [Synthetic]);
    const hi = [0, 20, 40, 60].map((k) => at(4000 + k / 20, 20, 100, 160));
    const li = [0, 20, 40, 60].map((k) => at(5_000_000 + k / 20, 20, 100, 160));
    for (const i of hi) await held.frame(i);
    for (const i of li) await lazy.frame(i);
    let best = Infinity;
    for (let i = 0; i < 3; i++) best = Math.min(best, median(lazy.renderer, li) / median(held.renderer, hi));
    console.info(`rows: lazy 10M against held 10k at letters, ratio ${best.toFixed(2)}`);
    expect(best).toBeLessThan(1.2);
  });
});
