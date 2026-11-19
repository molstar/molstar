/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { location } from '@molstar/seqx-core/query';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import { type GridLay, GridPlacement, gridSummary, recordDrawLists } from '../index.js';
import { settle } from './helpers.js';
import { gridSpec, letters } from './pairs-fixtures.js';

// A grid over a section of blocks (_spec/README.md "Resolved view", A grid over blocks; "Render", The grid layout and
// The grid layer): its resolution from the leader's blocks and its refusals, a strip a pair of blocks with its own lays,
// the chain boundaries' hairlines, and `gridSummary` (each cell once, the triangle kept, missing values left out).

const Size = { width: 720, height: 900 };

/** Blocks of lengths `ns` on systems `s0…`, a grid a pair of them (x-major, i·100 + j), the grid following them. */
function blocksSpec(ns: readonly number[], o: { arrange?: 'stack'; ownBlocks?: boolean; layout?: 'rings' } = {}): Node {
  const b = createSeqSpec();
  const seqs = ns.map((n, k) => b.data.sequence(`s${k}`, { alphabet: 'rna', letters: letters(n, k + 1) }));
  const grids: string[] = [];
  ns.forEach((nx, x) =>
    ns.forEach((ny, y) => {
      const v = new Float32Array(nx * ny);
      for (let i = 0; i < nx; i++) for (let j = 0; j < ny; j++) v[i * ny + j] = i * 100 + j;
      grids.push(b.data.grid(`g${x}${y}`, { x: seqs[x], y: seqs[y], values: v }));
    }),
  );
  b.views.view('v', (v) => {
    const s = v.section('chains', {
      ...(o.arrange ? { arrange: o.arrange } : {}),
      ...(o.layout ? { layout: o.layout } : {}),
    });
    seqs.forEach((q, k) => s.block(`b${k}`, { of: q, label: `chain ${k}` }));
    s.track('letters').layer('letters', { type: 'letters' });
    const g = v.section('plot', { layout: 'grid', follow: s.address });
    if (o.ownBlocks) g.block('own', { of: seqs[0] });
    g.track('plot', {}, (t) => grids.forEach((d, k) => t.layer(`g${k}`, { type: 'grid', data: d })));
  });
  return b.build();
}

const issuesOf = (spec: Node) => resolveViews(resolveSpec(spec))[0].issues.map((i) => `${i.code}: ${i.message}`);

describe('a grid over blocks: resolution', () => {
  it('takes the leader’s blocks and axis, its follow the identity on their units, its y the same blocks', () => {
    const view = resolveViews(resolveSpec(blocksSpec([10, 20])))[0];
    expect(view.issues).toEqual([]);
    const [lead, plot] = view.sections;
    expect(plot.blocks!.map((b) => b.axis.system)).toEqual(['/data/s0', '/data/s1']);
    expect(plot.axis).toBe(lead.axis);
    expect(plot.follow?.leader).toBe(lead.address);
    expect([plot.follow?.point(7), plot.follow?.scale]).toEqual([7, 1]);
    expect(plot.y?.axis).toBe(plot.axis);
    expect(plot.y?.followsX).toBe(true);
  });

  it('refuses a stack or a ring of blocks to follow, and a grid with blocks of its own', () => {
    for (const o of [{ arrange: 'stack' as const }, { layout: 'rings' as const }])
      expect(issuesOf(blocksSpec([10, 20], o))).toContain(
        'seq.view.section.system: A grid follows a section of blocks in a row only',
      );
    expect(issuesOf(blocksSpec([10, 20], { ownBlocks: true }))).toContain(
      'seq.view.section.system: A grid section has no blocks of its own: it follows a section of blocks',
    );
  });

  it('a follower that is no grid still cannot follow a section of blocks on several systems', async () => {
    const b = createSeqSpec();
    const s0 = b.data.sequence('s0', { alphabet: 'rna', letters: letters(10) });
    const s1 = b.data.sequence('s1', { alphabet: 'rna', letters: letters(20, 2) });
    b.views.view('v', (v) => {
      const s = v.section('chains', {});
      s.block('a', { of: s0 });
      s.block('b', { of: s1 });
      s.track('letters').layer('letters', { type: 'letters' });
      v.section('one', { of: s0, follow: s.address }).track('letters').layer('letters', { type: 'letters' });
    });
    const st = await settle(b.build());
    const f = await st.frame({ size: Size, settled: true });
    expect(f.layout.issues.map((i) => i.code)).toContain('seq.render.follow.blocks');
  });
});

describe('a grid over blocks: layout', () => {
  it('a strip a pair of blocks, each with lays of its own (kept by its y block); a grid on its own pair only', async () => {
    const st = await settle(blocksSpec([10, 20]));
    const f = await st.frame({ size: Size, settled: true });
    expect(f.layout.issues).toEqual([]);
    const g = f.layout.sections[1];
    expect(g.placement).toBeInstanceOf(GridPlacement);
    const strips = g.placement.strips;
    expect(strips.n).toBe(4);
    expect([0, 1, 2, 3].map((i) => `${strips.block[i]}${strips.yBlock[i]}@${strips.yFirst[i]}`)).toEqual([
      '00@0',
      '01@10',
      '10@0',
      '11@10',
    ]);
    // the two strips of x block 0 hold lays of their own (their y blocks differ); a strip lays the grids on its x block's
    // system (a layer placed per block), one drawn, the other on its pair of systems only, off here
    expect(g.lays[0][0]).not.toBe(g.lays[1][0]);
    for (let i = 0; i < 4; i++) {
      const lays = g.lays[i].map((l) => l.lay as GridLay);
      expect(lays.length).toBe(2);
      expect(lays.filter((l) => l.grid).length).toBe(1);
      expect(lays.filter((l) => l.off).length).toBe(1);
    }
  });

  it('three blocks: a hairline down the plot and one across it at each of the two boundaries', async () => {
    const st = await settle(blocksSpec([10, 20, 15]));
    const f = await st.frame({ size: Size, settled: true });
    const plot = f.layout.sections[1].placement;
    const hair = recordDrawLists(f.lists, ['labels'])
      .split('\n')
      .filter((l) => l.startsWith('rect '))
      .map((l) => l.split(' ').slice(1, 5).map(Number))
      .filter(([, y, w, h]) => y >= plot.box.y && (w === 1 || h === 1));
    expect(hair.filter(([, , w, h]) => w === 1 && h === plot.contentHeight).length).toBe(2);
    expect(hair.filter(([, , w, h]) => h === 1 && w > 100).length).toBe(2);
  });
});

describe('gridSummary', () => {
  it('each cell once under overlapping rects (an added drag), missing values left out', async () => {
    const values = new Float32Array(16).map((_, k) => k);
    values[5] = Number.NaN;
    const st = await settle(gridSpec({ n: 4, values }));
    const v = { model: st.model, view: st.view };
    const s = '/data/s';
    const one = location({ rects: [{ a: { system: s, ranges: [0, 2] }, b: { system: s, ranges: [0, 2] } }] });
    // cells (0,0) 0, (0,1) 1, (1,0) 4, (1,1) NaN
    expect(gridSummary(v, one)).toEqual({ mean: 5 / 3, min: 0, max: 4, count: 3 });
    const twice = location({
      rects: [
        { a: { system: s, ranges: [0, 2] }, b: { system: s, ranges: [0, 2] } },
        { a: { system: s, ranges: [1, 3] }, b: { system: s, ranges: [0, 2] } },
      ],
    });
    // the union: rows 0–2 × columns 0–1 less the missing (1,1): 0, 1, 4, 8, 9
    expect(gridSummary(v, twice)).toEqual({ mean: 22 / 5, min: 0, max: 9, count: 5 });
  });

  it('within the layer’s triangle; an empty rect, or one of only missing cells, counts nothing', async () => {
    const values = new Float32Array(16).map((_, k) => k);
    const st = await settle(gridSpec({ n: 4, values, triangle: 'upper' }));
    const v = { model: st.model, view: st.view };
    const s = '/data/s';
    const all = location({ rects: [{ a: { system: s, ranges: [0, 4] }, b: { system: s, ranges: [0, 4] } }] });
    // upper: y before x (j < i): (1,0) 4, (2,0) 8, (2,1) 9, (3,0) 12, (3,1) 13, (3,2) 14
    expect(gridSummary(v, all)).toEqual({ mean: 60 / 6, min: 4, max: 14, count: 6 });
    const none = gridSummary(v, location({}));
    expect(none.count).toBe(0);
    expect([none.mean, none.min, none.max].every(Number.isNaN)).toBe(true);
    const missing = new Float32Array(16).fill(Number.NaN);
    const st2 = await settle(gridSpec({ n: 4, values: missing }));
    expect(gridSummary({ model: st2.model, view: st2.view }, all).count).toBe(0);
  });
});
