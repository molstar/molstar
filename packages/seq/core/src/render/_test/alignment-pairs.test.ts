/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A family's consensus pairs and ranges on an alignment's columns, drawn over its rows (_spec/README.md "Render",
// The pairs layer, Layer types): each arc's feet over the columns its pair names, at the fit view and across an image
// tile's edge; a hover on a cell lights the arcs ending in its column, a hover on an arc reaches a chain's residues
// through its map onto a seed row and a residue's comes back to the arc; a selected pair drawn; the chain's ranges
// shown on the columns through that map.

import type { Node } from '@molstar/core/repr/tree';
import { cellsSource, type FormatRuntime } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { keySet, location, rangesLocation } from '@molstar/seqx-core/query';
import { projectOntoView } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { DrawLists, drawView, hitView, type PairsLay, type SectionLayout } from '../index.js';
import { ImageTile } from '../images.js';
import { settle } from './helpers.js';

const Rows = 40;
const Columns = 600;
const Size = { width: 720, height: 420 };

/** An alignment source of `Rows` rows over `Columns` columns: row i gaps column j where (i + j) % 7 is 0. */
const Family: FormatRuntime = {
  name: 'family',
  formats: ['family'],
  async open() {
    return cellsSource(Rows, Columns, (start, end) => {
      const codes = new Uint8Array((end - start) * Columns);
      for (let i = start; i < end; i++)
        for (let j = 0; j < Columns; j++)
          codes[(i - start) * Columns + j] = ((i + j) % 7 === 0 ? '-' : 'ACGU'[(i * 3 + j) % 4]).charCodeAt(0);
      const label = Array.from({ length: end - start }, (_, k) => `R${start + k}`);
      return { codes, label, first: new Int32Array(end - start).fill(1) };
    });
  },
};

/**
 * Nested stems around the middle, pairs whose ends lie either side of the first image tile's edge (256), and one end on
 * a column row 0 gaps (7: the chain has no residue there).
 */
const Pairs: [number, number][] = [
  [7, 120],
  [10, 590],
  [11, 589],
  [100, 300],
  [101, 299],
  [250, 262],
  [255, 256],
  [400, 500],
];

/** Row 0's residues onto the columns (its gapped columns 0, 7, 14, …): residue k at column k + 1 + floor(k / 6). */
function rowZero() {
  const from: number[] = [];
  const to: number[] = [];
  const length: number[] = [];
  for (let c = 1, k = 0; c < Columns; c += 7, k += 6) {
    from.push(k);
    to.push(c);
    length.push(Math.min(6, Columns - c));
  }
  return { from, to, length };
}

/** With `onChain`, the pairs drawn on the chain's section too: not their system, so said, never projected. */
function spec(o: { onChain?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'consensus pairs over a family' });
  const file = b.data.source('file', { format: 'family', uri: 'family:seed' });
  const seed = b.data.alignment('seed', { alphabet: 'rna', source: file });
  const pairs = b.data.pairs('ss', {
    a: seed,
    rows: {
      a_start: Pairs.map((p) => p[0]),
      a_end: Pairs.map((p) => p[0] + 1),
      b_start: Pairs.map((p) => p[1]),
      b_end: Pairs.map((p) => p[1] + 1),
      id: Pairs.map((p) => `${p[0]}-${p[1]}`),
      type: Pairs.map(() => 'canonical'),
    },
  });
  // the chain: row 0's residues, its stems on its own positions, mapped onto the columns through row 0
  const chain = b.data.sequence('chain', { alphabet: 'rna', letters: 'A'.repeat(514) });
  const stems = b.data.ranges('stems', { target: chain, rows: { start: [0, 12], end: [6, 18], id: ['s1', 's2'] } });
  b.data.map('onto', { from: chain, to: seed, method: 'segments', segments: rowZero() });
  b.views.view('v', (v) => {
    const s = v.section('family', { of: seed });
    s.track('arcs', { label: 'SS_cons' }).layer('arcs', { type: 'pairs', data: pairs });
    s.track('rows', { label: '', height: 200 }).layer('rows', { type: 'rows', data: seed });
    s.track('stems', { label: 'Chain stems' }).layer('stems', { type: 'ranges', data: stems });
    const c = v.section('chain', { of: chain });
    c.track('letters').layer('letters', { type: 'letters' });
    if (o.onChain) c.track('arcs').layer('arcs', { type: 'pairs', data: pairs });
  });
  return b.build();
}

// one settled family a file: the cases read it, none changes it
let shared: ReturnType<typeof settle> | undefined;
const settled = () => (shared ??= settle(spec(), [Family]));

const arcs = (s: SectionLayout) => s.sectionLays[0].lay as PairsLay;

/** The middle of track `t`'s band. */
const bandMiddle = (s: SectionLayout, t: number) => s.placement.strips.bandY[t] + s.placement.strips.bandH[t] / 2;

/** Each drawn arc's two feet, by pair: `[xa, xb]` with the pair's columns. */
function feet(lay: PairsLay): Map<string, [number, number]> {
  const out = new Map<string, [number, number]>();
  for (let i = 0; i < lay.n; i++) {
    const k = lay.pair[i];
    const d = lay.prep!.data!;
    out.set(`${d.aStart[k]}-${d.bStart[k]}`, [lay.curve[8 * i], lay.curve[8 * i + 6]]);
  }
  return out;
}

describe('pairs on an alignment’s columns, over its rows', () => {
  it('at the fit view and across an image tile’s edge: each foot over the column its pair names', async () => {
    const st = await settled();
    for (const viewport of [undefined, { start: 240, end: 272 }, { start: ImageTile - 3, end: ImageTile + 3 }]) {
      const f = await st.frame({ size: Size, ...(viewport ? { viewport } : {}) });
      expect(f.layout.issues).toEqual([]);
      const sec = f.layout.sections[0];
      const lay = arcs(sec);
      const y = bandMiddle(sec, 1);
      const a0 = sec.placement.strips.axes[0];
      const inView = (c: number) => c >= a0.u0 && c < a0.u1;
      const drawn = feet(lay);
      for (const [i, j] of Pairs) {
        const x = drawn.get(`${i}-${j}`);
        if (!inView(i) && !inView(j)) continue;
        expect(x, `${i}-${j}`).toBeDefined();
        // a foot in view is over its column of the rows (the rows' hit there names it)
        if (inView(i)) expect(hitView(f.layout, x![0], y)?.position).toBe(i);
        if (inView(j)) expect(hitView(f.layout, x![1], y)?.position).toBe(j);
      }
    }
  });

  it('a hover on a cell lights the arcs ending in its column; a hover on an arc and on a residue, carried through the map', async () => {
    const st = await settled();
    const f = await st.frame({ size: Size });
    const sec = f.layout.sections[0];
    const lay = arcs(sec);
    const x = feet(lay).get('100-300')![0];
    const cell = hitView(f.layout, x, bandMiddle(sec, 1))!;
    expect([cell.system, cell.position]).toEqual(['/data/seed', 100]);
    const lists = new DrawLists();
    drawView(f.layout, lists, { hover: cell }, ['hover']);
    expect(lists.passes.hover.paths.n).toBe(1);
    // the arc's apex: the pair, its extent its two columns
    const i = [...lay.pair.slice(0, lay.n)].findIndex((k) => lay.prep!.data!.aStart[k] === 400);
    const c = lay.curve;
    const o = 8 * i;
    // the cubic at t = 1/2
    const mid = (k: number) => (c[o + k] + 3 * c[o + 2 + k] + 3 * c[o + 4 + k] + c[o + 6 + k]) / 8;
    const apex = hitView(f.layout, mid(0), mid(1));
    expect(apex?.element).toEqual({ owner: '/data/ss', key: '400-500' });
    expect(apex?.extent?.parts.map((p) => [p.system, Array.from(p.ranges)])).toEqual([
      ['/data/seed', [400, 401, 500, 501]],
    ]);
    // carried onto the chain's section through the map (what a viewer projects), its residues lit there; a hover on
    // a residue carried back lights the arc ending at its column
    const onChain = projectOntoView(st.model, st.view, apex!.extent!);
    expect(onChain.parts.map((p) => [p.system, Array.from(p.ranges)])).toEqual([['/data/chain', [342, 343, 428, 429]]]);
    const projected = new DrawLists();
    drawView(f.layout, projected, { projected: { hover: onChain, hoverSection: sec.section.address } }, ['hover']);
    // the chain's section is the one lit: its two residues' cells, nothing over the family's
    const r = projected.passes.hover.rects;
    const chainBox = f.layout.sections[1].placement.box;
    const ax = f.layout.sections[1].placement.strips.axes[0];
    expect(r.n).toBe(2);
    for (let k = 0; k < r.n; k++) {
      expect(r.y[k]).toBeGreaterThanOrEqual(chainBox.y - 2);
      expect(r.y[k] + r.h[k]).toBeLessThanOrEqual(chainBox.y + chainBox.height + 2);
    }
    const centres = Array.from({ length: r.n }, (_, k) => r.x[k] + r.w[k] / 2).sort((a, b) => a - b);
    expect(centres[0]).toBeCloseTo(ax.xOfUnit(342.5), 0);
    expect(centres[1]).toBeCloseTo(ax.xOfUnit(428.5), 0);
    expect(projected.passes.hover.paths.n).toBe(0);
    // an end on a column row 0 gaps reaches no residue: the pair 7-120 lights its other end alone
    const gap = projectOntoView(
      st.model,
      st.view,
      location({ parts: [{ system: '/data/seed', ranges: [7, 8, 120, 121] }] }),
    );
    expect(gap.parts.map((p) => [p.system, Array.from(p.ranges)])).toEqual([['/data/chain', [102, 103]]]);
    const back = projectOntoView(st.model, st.view, rangesLocation('/data/chain', [85, 86]));
    expect(back.parts.map((p) => [p.system, Array.from(p.ranges)])).toEqual([['/data/seed', [100, 101]]]);
    const fromChain = new DrawLists();
    drawView(f.layout, fromChain, { projected: { hover: back, hoverSection: '/views/v/chain' } }, ['hover']);
    expect(fromChain.passes.hover.paths.n).toBe(1);
  });

  it('drawn on the chain’s section, the pairs on the columns are said (`seq.render.pairs.systems`), never projected', async () => {
    const st = await settle(spec({ onChain: true }), [Family]);
    const f = await st.frame({ size: Size });
    expect(f.layout.issues.map((i) => [i.code, i.address])).toEqual([
      ['seq.render.pairs.systems', '/views/v/chain/arcs/arcs'],
    ]);
    expect(arcs(f.layout.sections[1]).n).toBe(0);
    expect(arcs(f.layout.sections[0]).n).toBe(Pairs.length);
  });

  it('a selected pair is drawn in the accent; a selection on the columns draws over the rows', async () => {
    const st = await settled();
    const f = await st.frame({ size: Size });
    const lists = new DrawLists();
    const elements = keySet([{ owner: '/data/ss', key: '250-262' }]);
    const selection = location({ parts: [{ system: '/data/seed', ranges: [250, 251, 262, 263] }] });
    drawView(f.layout, lists, { selection, elements }, ['selection']);
    expect(lists.passes.selection.paths.n).toBe(1);
    expect(lists.passes.selection.rects.n).toBeGreaterThan(0);
  });

  it('a chain’s ranges drawn on the columns its seed row puts its residues at', async () => {
    const st = await settled();
    const f = await st.frame({ size: Size, viewport: { start: 0, end: 40 } });
    const sec = f.layout.sections[0];
    const y = bandMiddle(sec, 2);
    const a0 = sec.placement.strips.axes[0];
    // residues [0, 6) at columns [1, 7); residues [12, 18) at [15, 21): column 0 and 14 are row 0's gaps
    const key = (c: number) => hitView(f.layout, a0.xOfUnit(c + 0.5), y)?.element?.key;
    expect([key(0), key(1), key(6), key(7), key(14), key(15), key(20), key(21)]).toEqual([
      undefined,
      's1',
      's1',
      undefined,
      undefined,
      's2',
      's2',
      undefined,
    ]);
  });
});
