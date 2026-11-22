/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, type LogoLay, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import type { RowsLay } from '../rows.js';
import { settle } from './helpers.js';
import { Columns, rowsSpec, Section, Synthetic, Track } from './rows-fixtures.js';

const Labels = ['HBA/1-6', 'HBB/2-7', 'MYG/1-6', 'GLB/3-8', 'CYG/1-6'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLSDG', 'MKAVLG', 'MEKVPG'];
const Newick = '((HBA:1,HBB:1)hb:1,(MYG:1,(GLB:1,CYG:1):1):1)root;';

function spec(params: Record<string, unknown> = {}, track: Record<string, unknown> = {}) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: Labels, letters: Letters },
    reference: 'HBA/1-6',
  });
  const tree = b.data.tree('tree', { of: aln, newick: Newick });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }, (s) =>
      s.track('rows', { label: '', ...track }, (t) =>
        t.layer('r', { type: 'rows', data: aln, tree, ...params } as never),
      ),
    ),
  );
  return b;
}

const Size = { width: 480, height: 300 };

function rowsLay(f: {
  layout: { sections: readonly { lays: readonly (readonly { lay: unknown; type: { name: string } }[])[] }[] };
}) {
  for (const s of f.layout.sections)
    for (const l of s.lays) for (const x of l ?? []) if (x.type.name === 'rows') return x.lay as RowsLay;
  throw new Error('no rows layer');
}

describe('the rows layer', () => {
  it('draws every row in the band in data order, each cell hit as its row and column (letters)', async () => {
    const b = spec();
    expect(b.check()).toEqual([]);
    const s = await settle(b.build());
    const f = await s.frame({ size: Size, dpr: 1, settled: true });
    const lay = rowsLay(f as never);
    expect(lay.held?.rows).toBe(5);
    expect(lay.display.count).toBe(5);
    expect(lay.letters).toBe(true);
    // every row in the window has its band of cells
    const ys = new Set(Array.from(lay.y.subarray(0, lay.n)));
    expect(ys.size).toBe(5);
    expect(lay.labelText).toEqual(Labels);
    // every cell at letters is hit as its row and column
    for (let i = 0; i < 5; i++)
      for (let j = 0; j < 6; j++) {
        lay.axis.cell(j);
        const x = (lay.axis.sx0 + lay.axis.sx1) / 2;
        const y = lay.rowY(i) + lay.rowScale / 2;
        const h = hitView(f.layout, x, y);
        expect(h).toMatchObject({ target: 'cell', position: j, yPosition: i, ySystem: '/data/aln/(rows)' });
      }
    // a header's label is the row
    const h = hitView(f.layout, lay.labelX + 2, lay.rowY(2) + 3);
    expect(h).toMatchObject({ target: 'row', position: 2, system: '/data/aln/(rows)' });
  });

  it('in tree order folds clades from `collapsed`, and the band’s folds win; a node hits as `tree-node`', async () => {
    const s = await settle(spec({ order: 'tree', collapsed: ['hb'] }).build());
    const f = await s.frame({ size: Size, dpr: 1, settled: true });
    const lay = rowsLay(f as never);
    expect(lay.display.count).toBe(4);
    expect(lay.labelText[0]).toBe('HBA/1-6 +1');
    expect(lay.wedges.length).toBe(5);
    expect(hitView(f.layout, lay.labelX + 2, lay.rowY(0) + 3)).toMatchObject({
      target: 'summary-row',
      element: { owner: '/data/aln', key: 'hb' },
    });
    // a summary row's cells: the clade's consensus, opaque by its share
    expect(lay.n).toBeGreaterThan(0);
    // the band's folds (session state) win over `collapsed`
    const g = await s.frame({
      size: Size,
      dpr: 1,
      settled: true,
      viewports: { '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { bandTop: 0, folds: [] } } } },
    });
    const l2 = rowsLay(g as never);
    expect(l2.display.count).toBe(5);
    // a node within reach of its branch
    const [node, nx, ny] = l2.nodes;
    expect(hitView(g.layout, nx, ny)).toMatchObject({ target: 'tree-node', element: { key: l2.clades!.name[node] } });
  });

  it('a tree node’s target reaches 8 px from its branch (the junction’s), no further', async () => {
    const s = await settle(spec({ order: 'tree' }).build());
    const f = await s.frame({ size: Size, dpr: 1, settled: true });
    const lay = rowsLay(f as never);
    // the root: the node whose junction is leftmost, its branch the leftmost vertical
    let root = 0;
    for (let k = 0; k < lay.nodes.length; k += 3) if (lay.nodes[k + 1] < lay.nodes[root + 1]) root = k;
    const [node, rx, ry] = [lay.nodes[root], lay.nodes[root + 1], lay.nodes[root + 2]];
    expect(hitView(f.layout, rx - 7.5, ry)).toMatchObject({
      target: 'tree-node',
      element: { key: lay.clades!.name[node] },
    });
    expect(hitView(f.layout, rx - 8.5, ry)?.target).not.toBe('tree-node');
    // and stops at the label column: over a row's label the hit is the row's, whatever branch is near
    for (let k = 0; k < lay.nodes.length; k += 3)
      expect(hitView(f.layout, lay.labelX + 1, lay.nodes[k + 2])?.target).not.toBe('tree-node');
  });

  it('scrolls its rows in the band and zooms them; below 3 device px they draw as image tiles', async () => {
    const s = await settle(spec({}, { height: 40 }).build());
    const f = await s.frame({
      size: Size,
      dpr: 1,
      settled: true,
      viewports: { '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { bandTop: 1, rowScale: 20 } } } },
    });
    const lay = rowsLay(f as never);
    expect([lay.first, lay.last]).toEqual([1, 3]);
    expect(lay.labelText).toEqual(['HBB/2-7', 'MYG/1-6']);
    // the scroll is clamped to the rows
    const g = await s.frame({
      size: Size,
      dpr: 1,
      settled: true,
      viewports: { '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { bandTop: 999, rowScale: 20 } } } },
    });
    expect([rowsLay(g as never).top, rowsLay(g as never).bandTop]).toEqual([3, 60]);
    const h = await s.frame({
      size: Size,
      dpr: 1,
      settled: true,
      viewports: { '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { bandTop: 0, rowScale: 2 } } } },
    });
    const l3 = rowsLay(h as never);
    expect(l3.images.length).toBe(1);
    expect(l3.n).toBe(0);
  });

  it('pins a row in a band at the top (kept in place, marked) and hides rows behind a chip', async () => {
    const s = await settle(spec({ pinned: ['MYG/1-6'], hidden_rows: ['GLB/3-8'] }).build());
    const f = await s.frame({ size: Size, dpr: 1, settled: true });
    const lay = rowsLay(f as never);
    expect(lay.display.pinned).toBe(1);
    expect(lay.labelText).toEqual(['MYG/1-6', 'HBA/1-6', 'HBB/2-7', '· MYG/1-6', 'CYG/1-6']);
    expect(lay.chip).toBe('1 hidden');
  });

  it('the canvas and the export draw one spec alike: the same lists, the folds the export is given', async () => {
    // a band of its own height (a figure takes no fill)
    const b = spec({ order: 'tree', collapsed: ['hb'] }, { height: 120 });
    const s = await settle(b.build());
    const f = await s.frame({ size: Size, dpr: 1, settled: true });
    const canvas = recordDrawLists(f.lists, ['data', 'labels']);
    const svg = await renderSvg(b.build(), { size: Size });
    expect(svg.settled).toBe(true);
    expect(recordDrawLists(svg.lists, ['data', 'labels'])).toBe(canvas);
    expect(svg.svg).toContain('HBA/1-6 +1');
    // the viewer's own folds (none folded) passed to the export
    const open = await renderSvg(b.build(), { size: Size, folds: { '/views/v/s/rows': [] } });
    expect(open.svg).not.toContain('HBA/1-6 +1');
    expect(open.svg).toContain('HBB/2-7');
  });

  it('the logo stacks each column’s letters by information, lowest first, and hits a column with its letter', async () => {
    const b = createSeqSpec();
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
    b.views.view('v', (v) =>
      v.section('s', { of: aln }, (s) => {
        s.track('logo', { label: 'logo' }, (t) => t.layer('l', { type: 'logo', data: aln } as never));
      }),
    );
    expect(b.check()).toEqual([]);
    const st = await settle(b.build());
    const f = await st.frame({ size: Size, dpr: 1, settled: true });
    let lay: LogoLay | undefined;
    for (const sec of f.layout.sections)
      for (const l of sec.lays) for (const x of l ?? []) if (x.type.name === 'logo') lay = x.lay as LogoLay;
    // column 0 is all M: one letter, the column's full information (log2 20)
    const first = Array.from({ length: lay!.n }, (_, i) => i).filter((i) => lay!.lx[i] === lay!.lx[0]);
    expect(first.length).toBe(1);
    expect(String.fromCharCode(lay!.code[first[0]])).toBe('M');
    expect(lay!.vs[first[0]] * lay!.fontSize * 0.72).toBeCloseTo(lay!.box.height, 3);
    const h = hitView(f.layout, lay!.lx[0], lay!.box.y + lay!.box.height - 2);
    expect(h).toMatchObject({ position: 0, state: 'M' });
  });

  it('in text (an alignment figure) every row is on every line, unfolded, with no tree; each cell hits its row', async () => {
    const b = createSeqSpec();
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
    const tree = b.data.tree('tree', { of: aln, newick: Newick });
    b.views.view('v', (v) =>
      v.section('s', { of: aln, layout: 'text', group: 3, max_width: 48 } as never, (s) =>
        s
          .track('rows', { label: '' })
          .layer('r', { type: 'rows', data: aln, tree, order: 'tree', collapsed: ['hb'] } as never),
      ),
    );
    expect(b.check()).toEqual([]);
    const st = await settle(b.build());
    const f = await st.frame({ size: Size, dpr: 1, settled: true });
    const lays: RowsLay[] = [];
    for (const sec of f.layout.sections)
      for (const l of sec.lays) for (const x of l ?? []) if (x.type.name === 'rows') lays.push(x.lay as RowsLay);
    expect(lays.length).toBeGreaterThan(1);
    for (const lay of lays) {
      expect(lay.figure).toBe(true);
      expect(lay.display.count).toBe(5);
      expect(lay.tree.length).toBe(0);
      expect(lay.labelText.length).toBe(5);
    }
    expect(f.layout.issues.some((i) => i.code === 'seq.render.rows.text.unfolded')).toBe(true);
  });

  describe('the tree’s look: one line width on device pixels, a light wedge to the deepest leaf', () => {
    const nums = (l: string) =>
      l
        .split(' ')
        .filter((x) => x !== '' && !Number.isNaN(Number(x)))
        .map(Number);
    // a tree width whose thirds fall between pixels
    const params = { order: 'tree', collapsed: ['hb'], tree_width: 100 };

    /** The branch rects of the labels pass: the muted ink's rects left of the labels. */
    async function branches(dpr: number) {
      const s = await settle(spec(params, { height: 120 }).build());
      const f = await s.frame({ size: Size, dpr, settled: true });
      const lay = rowsLay(f as never);
      const rects = recordDrawLists(f.lists, ['labels'])
        .split('\n')
        .filter((l) => l.startsWith('rect '))
        .map(nums)
        .filter(([x]) => x < lay.labelX);
      return { lay, rects, f };
    }

    for (const dpr of [1, 2])
      it(`every branch rect’s edges on device pixels, every branch one width (dpr ${dpr})`, async () => {
        const { lay, rects } = await branches(dpr);
        // a segment a rect: the tree's segments in view
        expect(rects.length).toBe(lay.tree.length / 4);
        const onGrid = (v: number) => Math.abs(v * dpr - Math.round(v * dpr)) < 0.02;
        for (const [x, y, w, h] of rects) {
          for (const v of [x, y, x + w, y + h]) expect(onGrid(v), `${[x, y, w, h]}`).toBe(true);
          expect(Math.min(w, h)).toBeCloseTo(1, 5);
        }
      });

    it('a folded clade’s wedge ends at its deepest leaf’s x, a light fill outlined in the branch ink, canvas and SVG alike', async () => {
      const { lay, f } = await branches(1);
      // hb's leaves lie 2 from the root of a tree 3 deep (GLB, CYG): two thirds across
      // (the tree inset by half a stroke each side, so its strokes lie inside the column)
      const far = lay.treeX0 + 0.5 + ((lay.treeX1 - lay.treeX0 - 1) * 2) / 3;
      expect(lay.wedges[2]).toBeCloseTo(far, 5);
      const paths = recordDrawLists(f.lists, ['labels'])
        .split('\n')
        .filter((l) => l.startsWith('path '));
      expect(paths).toHaveLength(1);
      const xs = nums(paths[0].split(' stroke')[0]).filter((_, i) => i % 2 === 0);
      expect(Math.max(...xs)).toBeCloseTo(far, 1);
      expect(Math.max(...xs)).toBeLessThan(lay.treeX1 - 1);
      // its far side about half the row high (its corners softened inside it)
      const ys = nums(paths[0].split(' stroke')[0]).filter((_, i) => i % 2 === 1);
      const high = Math.max(...ys) - Math.min(...ys);
      expect(high).toBeGreaterThan(lay.rowScale * 0.4);
      expect(high).toBeLessThanOrEqual(lay.rowScale / 2);
      expect(ys.length).toBeGreaterThan(3);
      // filled with a tint (not opaque), outlined in the branch ink at the branches' width
      expect(paths[0]).toMatch(/stroke #[0-9a-f]{8} 1\b/);
      const fill = / #([0-9a-f]{8}) stroke /.exec(paths[0])?.[1];
      expect(fill).toBeTruthy();
      expect(Number.parseInt(fill!.slice(6), 16)).toBeLessThan(0x80);
      const svg = await renderSvg(spec(params, { height: 120 }).build(), { size: Size });
      expect(recordDrawLists(svg.lists, ['labels'])).toBe(recordDrawLists(f.lists, ['labels']));
    });
  });

  it('a folded clade whose node is at the tree’s right edge: its wedge still 6 px long, ending at the edge, and hit', async () => {
    const b = createSeqSpec();
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
    // hb's leaves at zero length below it, as deep as the tree goes: its node sits at the tree's right edge
    const tree = b.data.tree('tree', { of: aln, newick: '((HBA:0,HBB:0)hb:3,(MYG:1,(GLB:1,CYG:1):1):1)root;' });
    b.views.view('v', (v) =>
      v.section('s', { of: aln }, (s) =>
        s.track('rows', { label: '' }, (t) =>
          t.layer('r', { type: 'rows', data: aln, tree, order: 'tree', collapsed: ['hb'] } as never),
        ),
      ),
    );
    const s = await settle(b.build());
    const f = await s.frame({ size: Size, dpr: 1, settled: true });
    const lay = rowsLay(f as never);
    const [x0, y0, x1, y1] = lay.wedges;
    expect(x1).toBeCloseTo(lay.treeX1 - 0.5, 6);
    expect(x1 - x0).toBeGreaterThanOrEqual(6 - 1e-9);
    expect(hitView(f.layout, x0 + 2, (y0 + y1) / 2)).toMatchObject({ target: 'tree-node', element: { key: 'hb' } });
  });

  // 2,000 rows quick, the 10,000 medium (their settle a third of a second)
  for (const [Rows, tags] of [
    [2_000, []],
    [10_000, ['medium']],
  ] as [number, string[]][])
    describe(`the header below \`rows.labelMinPx\` on ${Rows} rows: the label column empty, the tree down to what resolves`, {
      tags,
    }, () => {
      // the rows and their tree settled once for the describe's frames
      let settled: ReturnType<typeof settle> | undefined;
      const scaled = () => (settled ??= settle(rowsSpec(Rows, true, { tree_width: 120 }), [Synthetic]));
      const at = (rowScale: number) =>
        ({
          size: { width: 900, height: 700 },
          dpr: 2,
          settled: true,
          viewports: { [Section]: { start: 0, end: Columns, bands: { [Track]: { bandTop: 0, rowScale } } } },
        }) as never;

      /** The header's draw items (the labels pass, clipped to the label column): kind and x. */
      const header = (f: { lists: Parameters<typeof recordDrawLists>[0] }, lay: RowsLay) =>
        recordDrawLists(f.lists, ['labels'])
          .split('\n')
          .map((l) => l.split(' '))
          // the track's own label (empty here) aside
          .filter(([kind, x]) => (kind === 'rect' || kind === 'text' || kind === 'path') && x !== '""')
          .filter(([, x]) => Number(x) < lay.box.x);

      it('at the band’s whole zoom: nothing in the label column, the tree’s few resolved clades only', async () => {
        const s = await scaled();
        const f = await s.frame(at(600 / Rows));
        const lay = rowsLay(f as never);
        expect(lay.rowScale).toBeLessThan(1);
        const items = header(f as never, lay);
        // the label column (past the tree): no label, number, chip or icon
        expect(items.filter(([, x]) => Number(x) >= lay.treeX1)).toEqual([]);
        expect(items.filter(([kind]) => kind === 'text')).toEqual([]);
        // the tree still drawn, down to the clades whose rows take 4 px: two segments a clade at most, a clade at least
        // 4 px of the band's 600 for each level, never one a row
        expect(lay.tree.length).toBeGreaterThan(0);
        expect(lay.tree.length / 4).toBeLessThan(2 * 2 * (600 / 4));
        expect(items.length).toBe(lay.tree.length / 4);
      });

      it('zoomed out, a node’s vertical spans its drawn branches only: each end meets a branch (none with no child drawn)', async () => {
        const s = await scaled();
        const f = await s.frame(at(600 / Rows));
        const lay = rowsLay(f as never);
        const t = lay.tree;
        const top = lay.box.y + lay.pinnedHeight;
        const bottom = lay.box.y + lay.box.height;
        // the branches' ys by the x they touch (either end)
        const meets = new Map<number, Set<number>>();
        const touching = (x: number) => meets.get(x) ?? meets.set(x, new Set()).get(x)!;
        for (let k = 0; k + 3 < t.length; k += 4)
          if (t[k + 1] === t[k + 3]) {
            touching(t[k]).add(t[k + 1]);
            touching(t[k + 2]).add(t[k + 1]);
          }
        let verticals = 0;
        for (let k = 0; k + 3 < t.length; k += 4) {
          if (t[k] !== t[k + 2] || t[k + 1] === t[k + 3]) continue;
          verticals++;
          const ys = meets.get(t[k]) ?? new Set();
          for (const y of [t[k + 1], t[k + 3]]) expect(ys.has(y) || y <= top || y >= bottom, `${t[k]} ${y}`).toBe(true);
        }
        expect(verticals).toBeGreaterThan(10);
      });

      it('at 8 px a row: the labels, the tree to every leaf', async () => {
        const s = await scaled();
        const f = await s.frame(at(8));
        const lay = rowsLay(f as never);
        const shown = Math.ceil(600 / 8);
        expect(lay.labelText.length).toBeGreaterThanOrEqual(shown - 1);
        expect(lay.labelText[0]).toMatch(/^R\d+$/);
        // every leaf in view has its branch
        const leaves = lay.treeNode.filter((n) => n < 0).length;
        expect(leaves).toBeGreaterThanOrEqual(shown - 1);
      });
    });
});
