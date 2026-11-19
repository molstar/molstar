/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView } from '../index.js';
import type { RowsLay } from '../rows.js';
import { settle } from './helpers.js';

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

const Size = { width: 400, height: 300 };

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

  it('scrolls its rows in the band and zooms them; below 3 device px they draw as image tiles', async () => {
    const s = await settle(spec({}, { height: 40 }).build());
    const f = await s.frame({
      size: Size,
      dpr: 1,
      settled: true,
      viewports: { '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { bandTop: 20, rowScale: 20 } } } },
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
    expect(rowsLay(g as never).bandTop).toBe(60);
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
});
