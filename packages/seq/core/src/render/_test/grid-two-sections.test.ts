/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { resolveViews } from '@molstar/seqx-core/view';
import { describe, expect, it } from 'vitest';
import {
  type GridLay,
  GridPlacement,
  hitView,
  type PairsLay,
  recordDrawLists,
  type SectionLayout,
  type Viewport,
} from '../index.js';
import { BottomSection, ChrLength, Chromosomes, sideBlocks, sidedSource, TopSection } from './bridge-sides-fixtures.js';
import { twoSectionsSpec } from './grid-two-sections-fixtures.js';
import { settle } from './helpers.js';

// A grid over two sections of blocks (EC-209; _spec/README.md "Resolved view", A grid over blocks; "Render", The grid
// layout): `follow_y` naming a second section of blocks gives the plot its y blocks; a strip a block pair, its cells
// the whole grid's clipped to both blocks; a reversed y block flips its side (EC-145); a dot plot of two blocks' pairs.

const Size = { width: 720, height: 900 };

const gridLays = (s: SectionLayout, i: number) =>
  s.lays[i].map((l) => l.lay as GridLay).filter((l) => 'grid' in l && l.grid !== undefined);

describe('a grid over two sections of blocks (EC-209): resolution', () => {
  it('takes x from the leader’s blocks and y from `follow_y`’s, each its identity follow', () => {
    const view = resolveViews(resolveSpec(twoSectionsSpec({ xs: [6, 8], ys: [5, 7, 4] })))[0];
    expect(view.issues).toEqual([]);
    const [a, b, plot] = view.sections;
    expect(plot.blocks!.map((k) => k.label)).toEqual(['s0', 's1']);
    expect(plot.y?.blocks?.map((k) => k.label)).toEqual(['t0', 't1', 't2']);
    expect(plot.y?.axis).toBe(b.axis);
    expect([plot.follow?.leader, plot.y?.follow?.leader]).toEqual([a.address, b.address]);
    expect([plot.y?.follow?.point(9), plot.y?.followsX]).toEqual([9, false]);
  });

  it('refuses a second section of blocks in a stack', () => {
    const view = resolveViews(resolveSpec(twoSectionsSpec({ xs: [6, 8], ys: [5, 7], stack: true })))[0];
    expect(view.issues.map((i) => `${i.code}: ${i.message}`)).toContain(
      'seq.view.section.system: A grid follows a section of blocks in a row only',
    );
  });
});

describe('a grid over two sections of blocks (EC-209): layout', () => {
  it('a strip a block pair (x by y), one grid drawn on each, the y blocks’ labels down the label column', async () => {
    const st = await settle(twoSectionsSpec({ xs: [6, 8], ys: [5, 7, 4] }));
    const f = await st.frame({ size: Size, settled: true });
    expect(f.layout.issues).toEqual([]);
    const g = f.layout.sections[2];
    const strips = g.placement.strips;
    expect(strips.n).toBe(6);
    expect([...Array(6).keys()].map((i) => `${strips.block[i]}${strips.yBlock[i]}`)).toEqual([
      '00',
      '01',
      '02',
      '10',
      '11',
      '12',
    ]);
    for (let i = 0; i < 6; i++) {
      const lays = gridLays(g, i);
      expect(lays.length).toBe(1);
      expect([lays[0].system, lays[0].ySystem]).toEqual([`/data/s${strips.block[i]}`, `/data/t${strips.yBlock[i]}`]);
    }
    const labels = recordDrawLists(f.lists, ['labels']);
    for (const t of ['t0', 't1', 't2']) expect(labels).toContain(t);
  });

  it('a block pair’s cells equal the whole grid’s clipped to both blocks (random windows on each side)', async () => {
    const xs = [6, 8];
    const ys = [5, 7, 4];
    const st = await settle(twoSectionsSpec({ xs, ys }));
    const r = new PCG(5);
    const windows: Record<string, Viewport>[] = [{}];
    for (let k = 0; k < 8; k++) {
      const w = (total: number) => {
        const a = (r.int() % (total * 10)) / 10;
        const b = Math.min(total, a + 2 + (r.int() % (total * 10)) / 10);
        return { start: a, end: b };
      };
      windows.push({ '/views/v/a': w(14), '/views/v/b': w(16) });
    }
    for (const viewports of windows) {
      const f = await st.frame({ size: Size, settled: true, viewports });
      const g = f.layout.sections[2];
      const s = g.placement.strips;
      for (let i = 0; i < s.n; i++) {
        const [lay] = gridLays(g, i);
        if (!lay) continue;
        const ax = s.axes[i];
        const ay = s.ys[i]!;
        const nx = xs[s.block[i]];
        const ny = ys[s.yBlock[i]];
        const box = [ax.x0, s.top[i], ax.x0 + ax.width, s.top[i] + s.height[i]];
        const expected: string[] = [];
        for (let x = 0; x < nx; x++)
          for (let y = 0; y < ny; y++) {
            if (!ax.cell(x)) continue;
            const x0 = Math.max(box[0], Math.min(ax.sx0, ax.sx1));
            const x1 = Math.min(box[2], Math.max(ax.sx0, ax.sx1));
            if (!ay.cell(y)) continue;
            const y0 = Math.max(box[1], Math.min(ay.sx0, ay.sx1));
            const y1 = Math.min(box[3], Math.max(ay.sx0, ay.sx1));
            if (x1 > x0 && y1 > y0) expected.push([x0, y0, x1 - x0, y1 - y0].map((v) => v.toFixed(3)).join(' '));
          }
        const drawn: string[] = [];
        for (let k = 0; k < lay.n; k++)
          drawn.push([...lay.rx.subarray(4 * k, 4 * k + 4)].map((v) => v.toFixed(3)).join(' '));
        expect(drawn.sort(), JSON.stringify(viewports)).toEqual(expected.sort());
      }
    }
  });

  it('a hit names the x block, the y block and the cell on both systems', async () => {
    const st = await settle(twoSectionsSpec({ xs: [6, 8], ys: [5, 7, 4] }));
    const f = await st.frame({ size: Size, settled: true });
    const s = f.layout.sections[2].placement.strips;
    const i = 4; // s1 × t1
    s.axes[i].cell(3);
    const x = (s.axes[i].sx0 + s.axes[i].sx1) / 2;
    s.ys[i]!.cell(2);
    const y = (s.ys[i]!.sx0 + s.ys[i]!.sx1) / 2;
    const h = hitView(f.layout, x, y)!;
    expect([h.block, h.yBlock, h.system, h.position, h.ySystem, h.yPosition]).toEqual([
      1,
      1,
      '/data/s1',
      3,
      '/data/t1',
      2,
    ]);
  });

  it('a reversed y block flips its side of the plot (EC-145)', async () => {
    const st = await settle(twoSectionsSpec({ xs: [6], ys: [5, 7], reversed: [0] }));
    const f = await st.frame({ size: Size, settled: true });
    const g = f.layout.sections[2];
    const p = g.placement as GridPlacement;
    expect(p).toBeInstanceOf(GridPlacement);
    const s = p.strips;
    const yOf = (i: number, j: number) => {
      s.ys[i]!.cell(j);
      return (s.ys[i]!.sx0 + s.ys[i]!.sx1) / 2;
    };
    // t0 reversed: its first position at its block's bottom; t1 forward: at its top
    expect(yOf(0, 0)).toBeGreaterThan(yOf(0, 4));
    expect(yOf(1, 0)).toBeLessThan(yOf(1, 6));
  });
});

describe('a dot plot of two blocks’ pairs (EC-209)', () => {
  it('runs drawn on their block pair alone, a run against the strand an anti-diagonal; nothing said elsewhere', async () => {
    const st = await settle(
      twoSectionsSpec({
        xs: [20, 30],
        ys: [25, 15],
        noGrids: true,
        pairs: {
          x: 1,
          y: 0,
          rows: [
            [2, 10, 4, 12, 1],
            [15, 25, 5, 15, -1],
            [27, 28, 20, 21, 1],
          ],
        },
      }),
    );
    const f = await st.frame({ size: Size, settled: true });
    expect(f.layout.issues).toEqual([]);
    const g = f.layout.sections[2];
    const s = g.placement.strips;
    for (let i = 0; i < s.n; i++) {
      // placed on its x block (s1) alone; its strip on t1 draws nothing
      const lay = g.lays[i][0]?.lay as PairsLay | undefined;
      const mine = s.block[i] === 1 && s.yBlock[i] === 0;
      expect([lay?.ns ?? 0, lay?.nd ?? 0], `strip ${i}`).toEqual(mine ? [2, 1] : [0, 0]);
      if (!mine || !lay) continue;
      // with the strand: x and y grow together; against it: y falls as x grows (an inversion, an anti-diagonal)
      const seg = (k: number) => [...lay.segs.subarray(5 * k, 5 * k + 4)];
      const [ax0, ay0, ax1, ay1] = seg(0);
      expect([ax1 > ax0, ay1 > ay0]).toEqual([true, true]);
      const [bx0, by0, bx1, by1] = seg(1);
      expect([bx1 > bx0, by1 < by0]).toEqual([true, true]);
    }
  });
});

describe('a dot plot of two genomes over a sided source (EC-209)', () => {
  const blocks = sideBlocks(30, 3);
  const spec = () => {
    const b = createSeqSpec({ title: 'dot plot' });
    const seqs = Chromosomes.map((name) => ({ name, length: ChrLength }));
    const P = b.data.assembly('P', assembly('P', seqs, { order: 'source' } as never));
    const G = b.data.assembly('G', assembly('G', seqs, { order: 'source' } as never));
    const src = b.data.source('src', { format: 'seq', value: seqToValue(sidedSource(blocks)) as never });
    const aln = b.data.pairs('aln', { a: P, b: G, source: src } as never);
    b.views.view('v', (v) => {
      genome(v, 'top', P);
      genome(v, 'bottom', G);
      v.section('plot', { layout: 'grid', follow: TopSection, follow_y: BottomSection } as never)
        .track('t')
        .layer('dots', { type: 'pairs', data: aln } as never);
    });
    return b.build();
  };
  const plotOf = (layout: { sections: readonly SectionLayout[] }) =>
    layout.sections.find((s) => s.section.address.endsWith('/plot'))!;

  it('items: each block drawn once, on the strip of its two chromosomes', async () => {
    const st = await settle(spec(), [SeqRuntime], { pairs: { maxArcs: 4096 } });
    const f = await st.frame({ size: Size, settled: true });
    expect(f.layout.issues).toEqual([]);
    const g = plotOf(f.layout);
    const s = g.placement.strips;
    let drawn = 0;
    for (let i = 0; i < s.n; i++) {
      const lay = g.lays[i][0].lay as PairsLay;
      const a = Chromosomes[s.block[i]];
      const b = Chromosomes[s.yBlock[i]];
      const mine = blocks.filter((x) => x.a === a && x.b === b);
      expect(lay.nd + lay.ns, `${a} × ${b}`).toBe(mine.length);
      drawn += lay.nd + lay.ns;
    }
    expect(drawn).toBe(blocks.length);
  });

  it('bins past the gate: a rect a bin pair, its opacity by count within `ribbons.binOpacity`', async () => {
    const st = await settle(spec(), [SeqRuntime], { pairs: { maxArcs: 4 } });
    const f = await st.frame({ size: Size, settled: true });
    const g = plotOf(f.layout);
    const s = g.placement.strips;
    let bins = 0;
    for (let i = 0; i < s.n; i++) {
      const lay = g.lays[i][0].lay as PairsLay;
      expect(lay.ns).toBe(0);
      for (let k = 0; k < lay.nd; k++) {
        const alpha = (lay.dotColor[k] & 0xff) / 255;
        expect(alpha).toBeGreaterThanOrEqual(0.15 - 0.01);
        expect(alpha).toBeLessThanOrEqual(0.8 + 0.01);
      }
      bins += lay.nd;
    }
    expect(bins).toBeGreaterThan(0);
  });
});
