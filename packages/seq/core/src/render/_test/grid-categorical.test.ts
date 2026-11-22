/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Value } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it, vi } from 'vitest';
import { colorOf } from '../grid.js';
import { type GridLay, GridPlacement, hitView, recordDrawLists } from '../index.js';
import { MemberHooks } from '../members-hooks.js';
import { Pyramid } from '../pyramid.js';
import { kinship } from './grid-categorical-fixtures.js';
import { settle } from './helpers.js';

// Categorical grid axes (EC-35; _spec/README.md "Resolved view", Grids: categorical axes; "Render", The grid layer):
// a grid section's `of`/`of_y` a rows system, its members in `member_order`/`member_order_y`; the axis transform holds
// the permutation, so cells, the reduction in display order and hits read members through it; the labels' chrome.

const Size = { width: 640, height: 900 };

/** The permutation a field order gives (ascending, ties in data order). */
const argsort = (v: readonly number[]) =>
  v
    .map((x, i) => [x, i] as const)
    .sort((a, b) => a[0] - b[0] || a[1] - b[1])
    .map((p) => p[1]);

describe('categorical grid axes (EC-35): resolution', () => {
  it('`of` and `of_y` a rows system: each axis its members’ order; no param left unread', () => {
    const m = resolveSpec(kinship(4, 1, { field: 'score', direction: 'descending' }, 'tree').spec);
    expect(m.issues).toEqual([]);
    const view = resolveViews(m)[0];
    expect(view.issues).toEqual([]);
    const s = view.sections[0];
    expect(s.members).toEqual({
      owner: '/data/aln',
      kind: 'rows',
      order: { field: 'score', direction: 'descending' },
    });
    expect(s.y?.members?.order).toBe('tree');
    expect(s.tracks[0].layers[0].placed).toBe(true);
  });

  it('a matrix’s states as both axes (samples beside a genotype matrix), labelled by its `states`', async () => {
    const b = createSeqSpec();
    const chr = b.data.sequence('chr', { alphabet: 'dna', letters: 'ACGTACGT' });
    const samples = ['NA2', 'NA0', 'NA1'];
    b.data.matrix('gt', { target: chr, states: samples, codes: new Uint8Array(24) } as never);
    const states = '/data/gt/(states)';
    const g = b.data.grid('kin', { x: states, y: states, values: new Float32Array(9).map((_, i) => i) } as never);
    b.views.view('v', (v) => {
      v.section('plot', { layout: 'grid', of: states, of_y: states, member_order: { field: 'label' } } as never)
        .track('k')
        .layer('g', { type: 'grid', data: g } as never);
    });
    const st = await settle(b.build());
    const f = await st.frame({ size: Size, settled: true });
    expect(f.layout.issues).toEqual([]);
    const a = f.layout.sections[0].placement.strips.axes[0];
    expect(a.members!.labels).toEqual(samples);
    expect([...a.members!.order]).toEqual([1, 2, 0]);
  });

  it('a categorical y over positions on x (genes × samples) keeps a window of its own, said nowhere', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
    const labels = ['a', 'b'];
    const aln = b.data.alignment('aln', { rows: { label: labels, letters: labels.map(() => 'AC') } } as never);
    const rows = `${aln}/(rows)`;
    const g = b.data.grid('e', { x: s, y: rows, values: new Float32Array(16) } as never);
    b.views.view('v', (v) => {
      v.section('plot', { layout: 'grid', of: s, of_y: rows } as never)
        .track('k')
        .layer('g', { type: 'grid', data: g } as never);
    });
    const view = resolveViews(resolveSpec(b.build()))[0];
    expect(view.issues).toEqual([]);
    expect([view.sections[0].members, view.sections[0].y?.followsX, view.sections[0].y?.members?.kind]).toEqual([
      undefined,
      false,
      'rows',
    ]);
  });
});

describe('categorical grid axes (EC-35): cells and hits through the permutation', () => {
  it('a property over random orders: the axis holds the order and its inverse, each cell its two members’ value', async () => {
    for (let seed = 1; seed <= 6; seed++) {
      const n = 6 + (seed % 5);
      const { spec, scores } = kinship(n, seed, { field: 'score' }, { field: 'score', direction: 'descending' });
      const st = await settle(spec);
      const f = await st.frame({ size: Size, settled: true });
      expect(f.layout.issues).toEqual([]);
      const s = f.layout.sections[0];
      const a = s.placement.strips.axes[0];
      const ay = s.placement.strips.ys[0]!;
      const up = argsort(scores);
      const down = argsort(scores.map((x) => -x));
      expect([...a.members!.order], `seed ${seed}`).toEqual(up);
      expect([...ay.members!.order], `seed ${seed}`).toEqual(down);
      for (let d = 0; d < n; d++) expect(a.members!.rank[a.members!.order[d]]).toBe(d);
      const lay = s.lays[0][0].lay as GridLay;
      expect(lay.n).toBe(n * n);
      for (let k = 0; k < lay.n; k++) {
        const x = lay.rx[4 * k] + lay.rx[4 * k + 2] / 2;
        const y = lay.rx[4 * k + 1] + lay.rx[4 * k + 3] / 2;
        // the display offsets under the cell, the members there, and the hit naming them
        const dx = Math.floor(a.unitAt(x));
        const dy = Math.floor(ay.unitAt(y));
        const h = hitView(f.layout, x, y)!;
        expect([h.position, h.yPosition]).toEqual([up[dx], down[dy]]);
        expect(lay.rc[k]).toBe(colorOf(lay, up[dx] * n + down[dy]));
      }
    }
  });

  it('a tree order: its leaves by label, rows the tree lacks after them in data order', async () => {
    const { spec } = kinship(5, 3, 'tree', 'data', { tree: '((S3,S1),(S4,S0));' });
    const st = await settle(spec);
    const f = await st.frame({ size: Size, settled: true });
    const s = f.layout.sections[0];
    expect([...s.placement.strips.axes[0].members!.order]).toEqual([3, 1, 4, 0, 2]);
    expect([...s.placement.strips.ys[0]!.members!.order]).toEqual([0, 1, 2, 3, 4]);
  });

  it('a reorder rebuilds the reduction in display order once; frames after it build none', async () => {
    const built = vi.spyOn(Pyramid.prototype, 'setValues');
    const n = 900;
    for (const order of [{ field: 'score' }, { field: 'score', direction: 'descending' }] as Value[]) {
      const st = await settle(kinship(n, 7, order).spec);
      built.mockClear();
      const lay = async (start: number, end: number) => {
        const f = await st.frame({ size: Size, settled: true, viewports: { '/views/v/plot': { start, end } } });
        return f.layout.sections[0].lays[0][0].lay as GridLay;
      };
      // image tiles (cells under `rectMinPx`): the reduction in display order, made once
      expect((await lay(0, n)).images.length).toBeGreaterThan(0);
      const once = built.mock.calls.length;
      expect(once).toBe(1);
      await lay(100, 700);
      await lay(0, n);
      expect(built.mock.calls.length).toBe(once);
    }
    built.mockRestore();
  });

  it('a warm pan looks no order up: the members are kept by the resolved axis', async () => {
    const st = await settle(kinship(300, 2).spec);
    const at = (k: number) => ({
      size: Size,
      settled: true,
      viewports: { '/views/v/plot': { start: 10 + (k % 2), end: 200 + (k % 2), y: { start: 0, end: 300 } } },
    });
    await st.frame(at(0));
    await st.frame(at(1));
    const made = vi.spyOn(MemberHooks, 'membersOf');
    for (let k = 0; k < 4; k++) st.renderer.frame(at(k));
    expect(made).not.toHaveBeenCalled();
    made.mockRestore();
  });
});

describe('categorical grid axes (EC-35): chrome', () => {
  it('members’ labels down the label column and at 45° along the top, thinned past `rows.labelMinPx`; no positions', async () => {
    const st = await settle(kinship(200, 2).spec);
    const f = await st.frame({ size: Size, settled: true });
    const s = f.layout.sections[0];
    const p = s.placement as GridPlacement;
    const a = p.strips.axes[0];
    // the box grew by the x labels' room above the plot
    expect(p.members!.room).toBeGreaterThan(0);
    expect(p.plotTop).toBeGreaterThanOrEqual(p.box.y + p.members!.room);
    const rec = recordDrawLists(f.lists, ['labels']);
    const ys = rec.split('\n').filter((l) => /^text "S\d+"/.test(l));
    // a label a step of at least 8 px: about the plot's height over 8
    expect(ys.length).toBeGreaterThan(5);
    expect(ys.length).toBeLessThanOrEqual(Math.ceil(p.plotHeight / 8) + 1);
    expect(p.yFirst).toBe('');
    // the x labels are placed letters, slanted
    expect(rec).toMatch(/glyphs/);
    expect(a.scale).toBeLessThan(8);
  });

  it('a click on a label is its member: a `row` hit on the rows system', async () => {
    const st = await settle(kinship(8, 4).spec);
    const f = await st.frame({ size: Size, settled: true });
    const s = f.layout.sections[0];
    const p = s.placement as GridPlacement;
    const a = p.strips.axes[0];
    const ay = p.strips.ys[0]!;
    const d = 3;
    const y = (ay.xOfUnit(d) + ay.xOfUnit(d + 1)) / 2;
    const left = hitView(f.layout, a.x0 - 20, y)!;
    expect([left.target, left.system, left.position]).toEqual(['row', '/data/aln/(rows)', ay.members!.order[d]]);
    const cx = (a.xOfUnit(d) + a.xOfUnit(d + 1)) / 2;
    // along the slanted label, a few px up its run
    const top = hitView(f.layout, cx + 5, p.plotTop - 4 - 5)!;
    expect([top.target, top.position]).toEqual(['row', a.members!.order[d]]);
  });
});

describe('categorical grid axes (EC-35): a reorder over 2,000 members', () => {
  itWithin('a frame after a reorder (the reduction in display order) under its bound', async (within) => {
    const n = 2000;
    const st = await settle(kinship(n, 9, { field: 'score' }).spec);
    await within(1500, () => st.renderer.frame({ size: Size, settled: true }));
  });
});
