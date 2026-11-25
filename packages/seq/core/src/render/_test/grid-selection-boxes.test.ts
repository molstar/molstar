/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { location, rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import { type DrawLists, type GridPlacement, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';
import { letters } from './pairs-fixtures.js';

// A selection's boxes on a grid over blocks (_spec/README.md "Render", The grid layer): each piece of the selection (a
// row band, a rect) lies inside its cell, stroke included, at a block boundary and at the grid's outer edge; a row
// band's fill and a column band's meet in the cell once.

const Size = { width: 720, height: 900 };
const Ns = [10, 20, 15];

/** Blocks of lengths `Ns` on systems `s0…`, a grid a pair of them, the grid following them on both axes. */
function spec(): Node {
  const b = createSeqSpec();
  const seqs = Ns.map((n, k) => b.data.sequence(`s${k}`, { alphabet: 'rna', letters: letters(n, k + 1) }));
  const grids: string[] = [];
  Ns.forEach((nx, x) =>
    Ns.forEach((ny, y) => {
      const v = new Float32Array(nx * ny);
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) v[i * ny + j] = i * 100 + j;
      grids.push(b.data.grid(`g${x}${y}`, { x: seqs[x], y: seqs[y], values: v }));
    }),
  );
  b.views.view('v', (v) => {
    const s = v.section('chains', {});
    seqs.forEach((q, k) => s.block(`b${k}`, { of: q, label: `chain ${k}` }));
    s.track('letters').layer('letters', { type: 'letters' });
    const g = v.section('plot', { layout: 'grid', follow: s.address });
    g.track('plot', {}, (t) => grids.forEach((d, k) => t.layer(`g${k}`, { type: 'grid', data: d })));
  });
  return b.build();
}

interface Box {
  x0: number;
  y0: number;
  x1: number;
  y1: number;
  fill: boolean;
}

/** The selection pass's rects in the plot, as the area they paint (a stroke centred on its edge reaches half out). */
function footprints(lists: DrawLists, top: number): Box[] {
  const out: Box[] = [];
  for (const l of recordDrawLists(lists, ['selection']).split('\n')) {
    if (!l.startsWith('rect ')) continue;
    const t = l.split(' ');
    const [x, y, w, h] = t.slice(1, 5).map(Number);
    if (y < top - 1) continue;
    const fill = t[5] !== 'none';
    const si = t.indexOf('stroke');
    const sw = si >= 0 ? Number(t[si + 2]) : 0;
    out.push({ x0: x - sw / 2, y0: y - sw / 2, x1: x + w + sw / 2, y1: y + h + sw / 2, fill });
  }
  return out;
}

/** The pairs of blocks' cells as `[x0, y0, x1, y1]`. */
function cells(plot: GridPlacement): number[][] {
  const xr = plot.xRow!;
  const yr = plot.yRow!;
  const out: number[][] = [];
  for (let i = 0; i < xr.n; i++)
    for (let j = 0; j < yr.n; j++) out.push([xr.x[i], yr.x[j], xr.x[i] + xr.w[i], yr.x[j] + yr.w[j]]);
  return out;
}

function expectInside(boxes: Box[], plot: GridPlacement, label: string) {
  const cs = cells(plot);
  expect(boxes.length, label).toBeGreaterThan(0);
  for (const b of boxes) {
    const inside = cs.some(
      (c) => b.x0 >= c[0] - 1e-6 && b.x1 <= c[2] + 1e-6 && b.y0 >= c[1] - 1e-6 && b.y1 <= c[3] + 1e-6,
    );
    expect(inside, `${label}: ${JSON.stringify(b)}`).toBe(true);
  }
}

/** No two fills overlap by an area (a cell's tint is drawn once). */
function expectTintedOnce(boxes: Box[], label: string) {
  const fills = boxes.filter((b) => b.fill);
  for (let i = 0; i < fills.length; i++)
    for (let j = i + 1; j < fills.length; j++) {
      const w = Math.min(fills[i].x1, fills[j].x1) - Math.max(fills[i].x0, fills[j].x0);
      const h = Math.min(fills[i].y1, fills[j].y1) - Math.max(fills[i].y0, fills[j].y0);
      expect(w > 1e-6 && h > 1e-6, `${label}: ${JSON.stringify([fills[i], fills[j]])}`).toBe(false);
    }
}

describe('a selection’s boxes on a grid over blocks', () => {
  const selections: [string, ReturnType<typeof rangesLocation>][] = [
    ['the first block (the grid’s left and top edge)', rangesLocation('/data/s0', [0, 10])],
    ['the middle block (both boundaries)', rangesLocation('/data/s1', [0, 20])],
    ['the last block (the grid’s right and bottom edge)', rangesLocation('/data/s2', [0, 15])],
    ['a range inside a block', rangesLocation('/data/s1', [3, 12])],
    [
      'a rect across the middle and the last block',
      location({ rects: [{ a: { system: '/data/s1', ranges: [5, 20] }, b: { system: '/data/s2', ranges: [0, 9] } }] }),
    ],
  ];

  for (const dpr of [1, 2])
    for (const [name, selection] of selections)
      it(`${name}: every piece inside its cell, tinted once (dpr ${dpr})`, async () => {
        const st = await settle(spec());
        const f = await st.frame({ size: Size, settled: true, dpr });
        const plot = f.layout.sections[1].placement as GridPlacement;
        const g = st.renderer.frame({ size: Size, dpr, emphasis: { selection } });
        const boxes = footprints(g.lists, plot.box.y);
        expectInside(boxes, plot, name);
        expectTintedOnce(boxes, name);
      });

  it('the export draws the same boxes', async () => {
    const r = await renderSvg(spec(), { size: Size, selection: '/data/s1/(range :from 4 :to 13)' });
    const st = await settle(spec());
    const f = await st.frame({ size: Size, settled: true });
    const plot = f.layout.sections[1].placement as GridPlacement;
    const boxes = footprints(r.lists, plot.box.y);
    expectInside(boxes, plot, 'export');
    expectTintedOnce(boxes, 'export');
  });
});
