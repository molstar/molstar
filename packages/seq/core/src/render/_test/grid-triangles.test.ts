/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node, Value } from '@molstar/core/repr/tree';
import { SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it, vi } from 'vitest';
import type { GridLay, PairsLay, SectionLayout } from '../index.js';
import { Pyramid } from '../pyramid.js';
import { ChrLength, Chromosomes, sideBlocks, sidedSource } from './bridge-sides-fixtures.js';
import { twoSectionsSpec } from './grid-two-sections-fixtures.js';
import { settle } from './helpers.js';
import { letters } from './pairs-fixtures.js';

// Triangles and reductions of grids across sections and on categorical axes (_spec/README.md "Render", The grid layer:
// triangles; categorical axes): a triangle halves one system against itself along the diagonal as drawn (a categorical
// axis's display order), never two sections or two systems; over one section of blocks a block pair is whole on its
// side; a layer's reductions are its own (two layers on one grid build one each); the dot plot's bins honour `mirror`
// and the triangle, their opacity rising with the count.

const Size = { width: 720, height: 900 };

const gridLays = (s: SectionLayout, i: number) =>
  s.lays[i].map((l) => l.lay as GridLay).filter((l) => 'grid' in l && l.grid !== undefined);

/** The cells a lay drew as rects, `i,j` by their display offsets on the strip. */
function drawnCells(s: SectionLayout, i: number, lay: GridLay): Set<string> {
  const a = s.placement.strips.axes[i];
  const ay = s.placement.strips.ys[i]!;
  const out = new Set<string>();
  for (let k = 0; k < lay.n; k++) {
    const dx = Math.floor(a.unitAt(lay.rx[4 * k] + lay.rx[4 * k + 2] / 2));
    const dy = Math.floor(ay.unitAt(lay.rx[4 * k + 1] + lay.rx[4 * k + 3] / 2));
    out.add(`${dx},${dy}`);
  }
  return out;
}

/** Blocks of lengths `ns` on one section, the grid over every pair following them on both axes, with `triangle`. */
function oneSection(ns: readonly number[], triangle: 'upper' | 'lower', layers: readonly string[] = ['g']): Node {
  const b = createSeqSpec();
  const seqs = ns.map((n, k) => b.data.sequence(`s${k}`, { alphabet: 'rna', letters: letters(n, k + 1) }));
  const grids: string[] = [];
  ns.forEach((nx, x) =>
    ns.forEach((ny, y) => {
      const v = new Float32Array(nx * ny).fill(x * 10 + y + 1);
      grids.push(b.data.grid(`g${x}${y}`, { x: seqs[x], y: seqs[y], values: v }));
    }),
  );
  b.views.view('v', (v) => {
    const s = v.section('chains', {});
    seqs.forEach((q, k) => s.block(`b${k}`, { of: q, label: `s${k}` }));
    s.track('letters').layer('letters', { type: 'letters' });
    v.section('plot', { layout: 'grid', follow: s.address }).track('plot', {}, (t) =>
      grids.forEach((d, k) => layers.forEach((l) => t.layer(`${l}${k}`, { type: 'grid', data: d, triangle }))),
    );
  });
  return b.build();
}

describe('a triangle across sections and blocks', () => {
  it('over two sections of blocks never halves a block pair', async () => {
    const xs = [6, 8];
    const ys = [5, 7, 4];
    const st = await settle(twoSectionsSpec({ xs, ys, triangle: 'upper' }));
    const f = await st.frame({ size: Size, settled: true });
    const g = f.layout.sections[2];
    const s = g.placement.strips;
    for (let i = 0; i < s.n; i++) expect(gridLays(g, i)[0].n, `strip ${i}`).toBe(xs[s.block[i]] * ys[s.yBlock[i]]);
  });

  it('over one section of blocks: a pair of two blocks whole on its side, none on the other; a block against itself halved, the upper the lower mirrored', async () => {
    const ns = [5, 7];
    const sets: Record<string, Set<string>[]> = {};
    for (const t of ['upper', 'lower'] as const) {
      const st = await settle(oneSection(ns, t));
      const f = await st.frame({ size: Size, settled: true });
      const g = f.layout.sections[1];
      const s = g.placement.strips;
      sets[t] = [];
      for (let i = 0; i < s.n; i++) {
        const bx = s.block[i];
        const by = s.yBlock[i];
        const [lay] = gridLays(g, i);
        const n = lay?.n ?? 0;
        if (bx === by) {
          // the strict half: n(n − 1) / 2 cells
          expect(n, `${t} ${bx}${by}`).toBe((ns[bx] * (ns[bx] - 1)) / 2);
          sets[t][bx] = drawnCells(g, i, lay);
        } else {
          const whole = t === 'upper' ? by < bx : by > bx;
          expect(n, `${t} ${bx}${by}`).toBe(whole ? ns[bx] * ns[by] : 0);
        }
      }
    }
    for (let b = 0; b < ns.length; b++) {
      const mirrored = new Set([...sets.upper[b]].map((c) => c.split(',').reverse().join(',')));
      expect(mirrored).toEqual(sets.lower[b]);
    }
  });

  it('on two systems with no blocks (chr1 × chr4) never halves the plot', async () => {
    const b = createSeqSpec();
    const c1 = b.data.sequence('c1', { alphabet: 'dna', letters: letters(12, 1) });
    const c4 = b.data.sequence('c4', { alphabet: 'dna', letters: letters(9, 4) });
    const g = b.data.grid('g', { x: c1, y: c4, values: new Float32Array(108).fill(1) });
    b.views.view('v', (v) => {
      v.section('plot', { layout: 'grid', of: c1, of_y: c4, follow_y: 'none' })
        .track('t')
        .layer('g', { type: 'grid', data: g, triangle: 'upper' });
    });
    const st = await settle(b.build());
    const f = await st.frame({ size: Size, settled: true });
    expect((f.layout.sections[0].lays[0][0].lay as GridLay).n).toBe(12 * 9);
  });
});

/** A kinship-like grid over `n` rows in a random field order, with `triangle`, in one or two sections (two orders). */
function kinship(n: number, o: { triangle?: 'upper'; orders?: readonly Value[]; layers?: number } = {}): Node {
  const r = new PCG(5);
  const scores = Array.from({ length: n }, () => r.int() % 1000);
  const b = createSeqSpec();
  const labels = Array.from({ length: n }, (_, i) => `S${i}`);
  const aln = b.data.alignment('aln', {
    rows: { label: labels, letters: labels.map(() => 'AC'), score: scores },
  } as never);
  const rows = `${aln}/(rows)`;
  const values = new Float32Array(n * n).map((_, k) => k);
  const g = b.data.grid('kin', { x: rows, y: rows, values } as never);
  b.views.view('v', (v) =>
    (o.orders ?? [{ field: 'score' }]).forEach((order, s) =>
      v
        .section(`plot${s}`, {
          layout: 'grid',
          of: rows,
          of_y: rows,
          member_order: order,
          member_order_y: order,
        } as never)
        .track('k', {}, (t) => {
          for (let l = 0; l < (o.layers ?? 1); l++)
            t.layer(`g${l}`, {
              type: 'grid',
              data: g,
              ...(o.triangle ? { triangle: l === 0 ? o.triangle : 'lower' } : {}),
            } as never);
        }),
    ),
  );
  return b.build();
}

describe('a triangle on categorical axes', () => {
  it('keeps the half above the drawn diagonal: by display order, not by member index', async () => {
    const st = await settle(kinship(9, { triangle: 'upper' }));
    const f = await st.frame({ size: Size, settled: true });
    const s = f.layout.sections[0];
    const cells = drawnCells(s, 0, s.lays[0][0].lay as GridLay);
    expect(cells.size).toBe((9 * 8) / 2);
    // upper: y before x in display order
    for (const c of cells) {
      const [dx, dy] = c.split(',').map(Number);
      expect(dy).toBeLessThan(dx);
    }
  });
});

describe('a grid layer’s reductions', () => {
  it('two layers on one grid (upper and lower), frames zoomed out: one reduction built each, never again', async () => {
    const built = vi.spyOn(Pyramid.prototype, 'setValues');
    const st = await settle(kinship(900, { triangle: 'upper', layers: 2 }));
    built.mockClear();
    for (let k = 0; k < 4; k++) {
      const f = await st.frame({
        size: Size,
        settled: true,
        viewports: { '/views/v/plot0': { start: 10 * k, end: 900 - 10 * k } },
      });
      const lays = f.layout.sections[0].lays[0].map((l) => l.lay as GridLay);
      expect(lays.every((l) => l.images.length > 0)).toBe(true);
    }
    expect(built.mock.calls.length).toBe(2);
    built.mockRestore();
  });

  it('one grid in a data-order section and a field-order section: one reduction each, never again', async () => {
    const built = vi.spyOn(Pyramid.prototype, 'setValues');
    const st = await settle(kinship(900, { orders: ['data', { field: 'score' }] }));
    built.mockClear();
    for (let k = 0; k < 3; k++) await st.frame({ size: { width: 600, height: 2400 }, settled: true });
    expect(built.mock.calls.length).toBe(2);
    built.mockRestore();
  });
});

describe('follows onto a categorical axis, and onto no blocks', () => {
  it('a `follow` or `follow_y` from positions onto members is `seq.view.follow.unplaced`', () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'ACGTACGT' });
    const labels = ['a', 'b', 'c'];
    const aln = b.data.alignment('aln', { rows: { label: labels, letters: labels.map(() => 'AC') } } as never);
    const rows = `${aln}/(rows)`;
    b.views.view('v', (v) => {
      v.section('pos', { of: s }).track('l').layer('l', { type: 'letters' });
      v.section('x', { layout: 'grid', of: rows, of_y: rows, follow: '/views/v/pos' } as never);
      v.section('y', { layout: 'grid', of: s, of_y: rows, follow_y: '/views/v/pos' } as never);
    });
    const codes = resolveViews(resolveSpec(b.build()))[0].issues.map((i) => `${i.code} ${i.address}`);
    expect(codes).toContain('seq.view.follow.unplaced /views/v/x');
    expect(codes).toContain('seq.view.follow.unplaced /views/v/y');
  });

  it('a grid over blocks whose `follow_y` names a section with no blocks says so, y the blocks across', () => {
    const b = createSeqSpec();
    const s0 = b.data.sequence('s0', { alphabet: 'dna', letters: 'ACGTACGT' });
    b.views.view('v', (v) => {
      const a = v.section('a', {});
      a.block('b0', { of: s0 });
      a.track('l').layer('l', { type: 'letters' });
      v.section('plain', { of: s0 }).track('l').layer('l', { type: 'letters' });
      v.section('plot', { layout: 'grid', follow: a.address, follow_y: '/views/v/plain' } as never);
    });
    const view = resolveViews(resolveSpec(b.build()))[0];
    expect(view.issues.map((i) => `${i.code}: ${i.message}`)).toContain(
      'seq.view.follow.unplaced: The grid’s y follows no section of blocks: it shows the blocks across',
    );
    expect(view.sections[2].y?.followsX).toBe(true);
  });
});

describe('a dot plot’s bins', () => {
  const blocks = sideBlocks(60, 4);
  /** P against itself over one section of blocks (x and y its blocks), a sided self source, bins (`maxArcs` 4). */
  const spec = (layer: Record<string, unknown>) => {
    const b = createSeqSpec();
    const seqs = Chromosomes.map((name) => ({ name, length: ChrLength }));
    const P = b.data.assembly('P', assembly('P', seqs, { order: 'source' } as never));
    const src = b.data.source('src', { format: 'seq', value: seqToValue(sidedSource(blocks)) as never });
    const aln = b.data.pairs('aln', { a: P, b: P, source: src } as never);
    b.views.view('v', (v) => {
      genome(v, 'top', P);
      v.section('plot', { layout: 'grid', follow: '/views/v/top' } as never)
        .track('t')
        .layer('dots', { type: 'pairs', data: aln, ...layer } as never);
    });
    return b.build();
  };
  const binsOf = async (layer: Record<string, unknown>) => {
    const st = await settle(spec(layer), [SeqRuntime], { pairs: { maxArcs: 4 } });
    const f = await st.frame({ size: Size, settled: true });
    const g = f.layout.sections.find((s) => s.section.address.endsWith('/plot'))!;
    return { g, lays: g.lays.map((ls) => ls[0]?.lay as PairsLay | undefined) };
  };

  it('honour `mirror` (each off-diagonal bin twice) and the triangle (above the drawn diagonal alone)', async () => {
    const count = (lays: (PairsLay | undefined)[]) => lays.reduce((n, l) => n + (l?.nd ?? 0), 0);
    const plain = await binsOf({});
    const mirrored = await binsOf({ mirror: true });
    expect(count(plain.lays)).toBeGreaterThan(0);
    expect(count(mirrored.lays)).toBeGreaterThan(count(plain.lays));
    const upper = await binsOf({ mirror: true, triangle: 'upper' });
    const s = upper.g.placement.strips;
    let kept = 0;
    upper.lays.forEach((l, i) => {
      for (let k = 0; l && k < l.nd; k++) {
        kept++;
        const cx = l.dots[4 * k] + l.dots[4 * k + 2] / 2;
        const cy = l.dots[4 * k + 1] + l.dots[4 * k + 3] / 2;
        if (s.block[i] === s.yBlock[i]) expect(cy - s.top[i]).toBeLessThan(cx - s.axes[i].x0 + 1e-6);
        else expect(s.yBlock[i]).toBeLessThan(s.block[i]);
      }
    });
    expect(kept).toBeGreaterThan(0);
    expect(kept).toBeLessThan(count(mirrored.lays));
  });

  it('their opacity rises with the count', async () => {
    const { lays } = await binsOf({});
    let checked = 0;
    for (const l of lays) {
      if (!l || l.nd < 2) continue;
      const w = l.prep!.weight!;
      const by = [...Array(l.nd).keys()].sort((p, q) => w[l.dotPair[p]] - w[l.dotPair[q]]);
      for (let k = 1; k < by.length; k++) {
        const a0 = l.dotColor[by[k - 1]] & 0xff;
        const a1 = l.dotColor[by[k]] & 0xff;
        expect(a1).toBeGreaterThanOrEqual(a0);
        checked++;
      }
    }
    expect(checked).toBeGreaterThan(0);
  });
});
