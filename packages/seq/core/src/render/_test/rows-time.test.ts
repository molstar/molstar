/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows' targets (the design's "Performance and size"): over a 10,000 × 300 alignment, a first frame (the join, the
// measures, the order, the theme's classes, the reduction, the layout) under 300 ms; a pan frame under 8 ms at the fit
// (image tiles), at letters and between two level pairs; a theme change's rebuild under 50 ms; folding every clade
// under 16 ms; and a frame of 10,000,000 lazy rows at letters within 1.2× of the held one. Each bound is the heavy
// lane's (`itWithin`, best of three); the quick lane runs a fit frame and a letters frame unbounded.

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { cellsSource, type FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { DarkTokens, type FrameInput, LightTokens, Renderer, type RowsLay } from '../index.js';
import { settle } from './helpers.js';

const Residues = 'ACDEFGHIKLMNPQRSTVWY-';
const Columns = 300;
const Track = '/views/v/s/rows';
const Section = '/views/v/s';

/** A synthetic alignment source (`synthetic`): rows made by index, each a biased draw a column. */
const Synthetic: FormatRuntime = {
  name: 'synthetic',
  formats: ['synthetic'],
  async open(source) {
    const n = (source.value as { rows: number }).rows;
    return cellsSource(n, Columns, (start, end) => {
      const codes = new Uint8Array((end - start) * Columns);
      for (let i = start; i < end; i++) {
        const rng = new PCG((i * 2654435761) >>> 0);
        for (let j = 0; j < Columns; j++)
          codes[(i - start) * Columns + j] = Residues.charCodeAt(
            rng.float() < 0.6 ? (j * 7) % Residues.length : rng.int() % Residues.length,
          );
      }
      const label = Array.from({ length: end - start }, (_, k) => `R${start + k}`);
      return { codes, label, first: new Int32Array(end - start).fill(1) };
    });
  },
};

/** A balanced tree over the rows, its leaves `R<i>`. */
function balanced(n: number): string {
  const out: string[] = [];
  const grow = (a: number, b: number) => {
    if (b - a === 1) return void out.push(`R${a}:1`);
    const m = (a + b) >> 1;
    out.push('(');
    grow(a, m);
    out.push(',');
    grow(m, b);
    out.push('):1');
  };
  grow(0, n);
  return `${out.join('')};`;
}

function spec(rows: number, tree: boolean): Node {
  const b = createSeqSpec();
  const src = b.data.source('src', { format: 'synthetic', value: { rows } });
  const aln = b.data.alignment('aln', { alphabet: 'protein', source: src });
  const t = tree ? b.data.tree('tree', { of: aln, newick: balanced(rows) }) : undefined;
  b.views.view('v', (v) =>
    v
      .section('s', { of: aln })
      .track('rows', { label: '', height: 600 })
      .layer('r', { type: 'rows', data: aln, ...(t ? { tree: t, order: 'tree' } : {}) } as never),
  );
  return b.build();
}

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
      const s = await settle(spec(10_000, true), [Synthetic]);
      await s.frame(at(0, 0.06));
      const first = () => {
        const r = new Renderer({
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
    const s = await settle(spec(10_000, true), [Synthetic]);
    const fit = [0, 1, 2, 3].map((k) => at(k, 0.06));
    const letters = [0, 20, 40, 60].map((k) => at(20 * 4000 + k, 20, 100, 160));
    const pairs = [0, 4, 8, 12].map((k) => at(k * 10, 0.5, 0, Columns));
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
    const s = await settle(spec(10_000, false), [Synthetic]);
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
    const s = await settle(spec(10_000, true), [Synthetic]);
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
    const hi = [0, 20, 40, 60].map((k) => at(20 * 4000 + k, 20, 100, 160));
    const li = [0, 20, 40, 60].map((k) => at(20 * 5_000_000 + k, 20, 100, 160));
    for (const i of hi) await held.frame(i);
    for (const i of li) await lazy.frame(i);
    let best = Infinity;
    for (let i = 0; i < 3; i++) best = Math.min(best, median(lazy.renderer, li) / median(held.renderer, hi));
    console.info(`rows: lazy 10M against held 10k at letters, ratio ${best.toFixed(2)}`);
    expect(best).toBeLessThan(1.2);
  });
});
