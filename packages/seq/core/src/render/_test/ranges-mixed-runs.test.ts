/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A merged run of items of mixed categories (_spec/README.md "Render", Layer types, Ranges): it takes its majority
// category's fill (ties by the legend's order, never the neutral colour), and the hover over it lists its categories.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, type RangesLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay, readRects } from './small-helpers.js';

const SIZE = { width: 400, height: 40 } as const;
const A = 0xaa0000;
const B = 0x0000bb;
const ORDER = ['alpha', 'beta'] as const;

/** Items one position wide, adjacent, over 4000 positions in 400 px (a tenth of a pixel each: one merged bar). */
function spec(kinds: string[], values: readonly string[] = ORDER, widths?: number[]): Node {
  const b = createSeqSpec({ title: 'mixed runs' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGT'.repeat(1000) });
  // `widths`: item i spans `widths[i]` positions, adjacent to the next
  const start: number[] = [];
  const end: number[] = [];
  let at = 100;
  kinds.forEach((_, i) => {
    start.push(at);
    at += widths?.[i] ?? 1;
    end.push(at);
  });
  const rows = {
    start,
    end,
    id: kinds.map((_, i) => `i${i}`),
    kind: kinds,
  };
  const items = b.data.ranges('items', { target: s, rows: rows as never });
  b.views.view('main', (v) => {
    v.section('s', { of: s })
      .track('items')
      .layer('items', {
        type: 'ranges',
        data: items,
        color: theme.field('kind', {
          palette: theme.categorical({
            colors: Object.fromEntries(values.map((v) => [v, v === 'alpha' ? '#aa0000' : '#0000bb'])),
          }),
        }),
      } as never);
  });
  return b.build();
}

async function drawn(kinds: string[], values?: readonly string[], widths?: number[]) {
  const f = await (await settle(spec(kinds, values, widths))).frame({ size: SIZE, settled: true });
  const lay = layerLay<RangesLay>(f, 0);
  const rects = readRects(f.lists).filter((r) => r.y >= lay.box.y && r.y < lay.box.y + lay.box.height);
  return { f, lay, rects };
}

const rgb = (fill: number) => fill >>> 8;

describe('ranges: a merged run of mixed categories', () => {
  it('takes the majority category’s fill, never the neutral colour', async () => {
    const { lay, rects } = await drawn(['beta', 'beta', 'alpha', 'beta', 'alpha']);
    expect(lay.step).toBe('bars');
    expect(lay.m).toBe(1);
    expect(rects.map((r) => rgb(r.fill))).toEqual([B]);
  });

  it('weighs each item by the width it covers: a long item outweighs several short ones', async () => {
    // a 100-position beta (10 px) among four 1-position alphas (a pixel each at least): beta covers most of it
    const { lay, rects } = await drawn(['alpha', 'alpha', 'beta', 'alpha', 'alpha'], undefined, [1, 1, 100, 1, 1]);
    expect(lay.m).toBe(1);
    expect(rects.map((r) => rgb(r.fill))).toEqual([B]);
  });

  it('equal total widths from different sums are a tie: the legend’s first category', async () => {
    // beta 3 px; alpha 1 px and 2 px: the sums differ in the last float digit, not as widths
    const { rects } = await drawn(['beta', 'alpha', 'alpha'], undefined, [30, 10, 20]);
    expect(rects.map((r) => rgb(r.fill))).toEqual([A]);
  });

  it('a tie takes the legend’s first category', async () => {
    const tie = await drawn(['beta', 'alpha']);
    expect(tie.rects.map((r) => rgb(r.fill))).toEqual([A]);
    const flipped = await drawn(['alpha', 'beta'], ['beta', 'alpha']);
    expect(flipped.rects.map((r) => rgb(r.fill))).toEqual([B]);
  });

  it('a run of one category is unchanged', async () => {
    const { lay, rects } = await drawn(['beta', 'beta', 'beta']);
    expect(lay.m).toBe(1);
    expect(rects.map((r) => rgb(r.fill))).toEqual([B]);
  });

  it('the hover over a mixed run lists its categories in the legend’s order; over one category, none', async () => {
    const mixed = await drawn(['beta', 'alpha', 'beta']);
    const x = (mixed.lay.mx0[0] + mixed.lay.mx1[0]) / 2;
    const y = mixed.lay.box.y + mixed.lay.box.height / 2;
    expect(hitView(mixed.f.layout, x, y)?.categories).toBe('alpha, beta');
    const one = await drawn(['beta', 'beta']);
    const h = hitView(one.f.layout, (one.lay.mx0[0] + one.lay.mx1[0]) / 2, one.lay.box.y + one.lay.box.height / 2);
    expect(h?.element).toBeDefined();
    expect(h?.categories).toBeUndefined();
  });
});
