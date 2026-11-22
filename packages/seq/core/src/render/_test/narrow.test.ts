/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView } from '../index.js';
import { DefaultRenderConfig, renderConfig } from '../config.js';
import { narrowColumns, overviewWidth } from '../narrow.js';
import type { RowsLay } from '../rows.js';
import { settle } from './helpers.js';

// A narrow frame (_spec/README.md "Render", A narrow frame): the vertical overview narrows 96 → 48 px and hides below
// a 160 px plot, then the label column shrinks to its floor and the tree collapses, then the plot keeps its floor and
// the frame scrolls (`minWidth`); a wide frame keeps every column.

const Labels = ['HBA/1-6', 'HBB/2-7', 'MYG/1-6', 'GLB/3-8', 'CYG/1-6'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLSDG', 'MKAVLG', 'MEKVPG'];
const Newick = '((HBA:1,HBB:1)hb:1,(MYG:1,(GLB:1,CYG:1):1):1)root;';

function spec() {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick: Newick });
  b.views.view('v', (v) =>
    v.section('s', { of: aln }, (s) =>
      s.track('rows', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln, tree, order: 'tree' } as never)),
    ),
  );
  return b.build();
}

async function laid(width: number, treeOpen = false) {
  const s = await settle(spec());
  const f = await s.frame({
    size: { width, height: 200 },
    dpr: 1,
    settled: true,
    ...(treeOpen
      ? { viewports: { '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { bandTop: 0, tree: true } } } } }
      : {}),
  });
  let lay: RowsLay | undefined;
  for (const l of f.layout.sections[0].lays)
    for (const x of l ?? []) if (x.type.name === 'rows') lay = x.lay as RowsLay;
  const box = f.layout.sections[0].tracks[0].box;
  return { layout: f.layout, lay: lay!, box };
}

describe('a narrow frame', () => {
  const config = DefaultRenderConfig;

  it('the overview: its size while the plot keeps 160 px beside it, narrowed to 48, hidden below', () => {
    expect(overviewWidth(800, 96, config)).toBe(96);
    expect(overviewWidth(160 + 96, 96, config)).toBe(96);
    expect(overviewWidth(160 + 70, 96, config)).toBe(70);
    expect(overviewWidth(160 + 48, 96, config)).toBe(48);
    expect(overviewWidth(160 + 47, 96, config)).toBe(0);
    // `minPlotPx: 0` never hides; a narrower `minPx` keeps it longer
    expect(overviewWidth(10, 96, renderConfig({ navigator: { overview: { minPlotPx: 0 } } as never }))).toBe(96);
    expect(overviewWidth(160 + 30, 96, renderConfig({ navigator: { overview: { minPx: 24 } } as never }))).toBe(30);
  });

  it('the columns: the label column to its floor first, then the tree; the least frame width keeps the plot', () => {
    // a plot of 300 px: nothing gives way
    expect(narrowColumns(600, 300, 112, 120, config)).toMatchObject({ labelWidth: 112, tree: 120, plotPx: 300 });
    // 20 px short: the label column gives them
    expect(narrowColumns(440, 300, 112, 120, config)).toMatchObject({ labelWidth: 92, tree: 120 });
    // 50 px short: the label column to 80, the tree collapses
    const n = narrowColumns(410, 300, 112, 120, config);
    expect(n).toMatchObject({ labelWidth: 80, tree: 0, plotPx: 110, columnsPx: 300 });
    expect(n.minWidth).toBe(300 - 32 - 120 + 160);
  });

  it('at 1280 px an alignment with its tree keeps every column', async () => {
    const { layout, lay, box } = await laid(1280);
    expect(layout.narrow).toMatchObject({ labelWidth: 112 });
    expect(layout.narrow!.plotPx).toBeGreaterThan(160);
    expect(lay.treeX1 - lay.treeX0).toBe(120);
    expect(box.width).toBeCloseTo(layout.narrow!.plotPx, 6);
  });

  it('at 480 px the label column shrinks toward its floor, narrower the tree collapses; the least width keeps a 160 px plot', async () => {
    const wide = await laid(1280);
    const columns = wide.layout.narrow!.columnsPx;
    // at 480 px the columns (112 px labels, a 120 px tree, the numbers) leave the plot under 160 px
    const { layout, lay, box } = await laid(480);
    expect(480 - columns).toBeLessThan(160);
    expect(layout.narrow!.labelWidth).toBeLessThan(112);
    expect(layout.narrow!.labelWidth).toBeGreaterThanOrEqual(80);
    expect(layout.narrow!.labelWidth).toBe(Math.max(80, Math.floor(112 - (160 - (480 - columns)))));
    // the data box is what the narrowed columns leave
    expect(box.width).toBeGreaterThan(480 - columns);
    // narrower, past the label column's floor: the tree collapses
    const narrower = await laid(columns - 32 + 160 - 40);
    expect(narrower.layout.narrow).toMatchObject({ labelWidth: 80, tree: 0 });
    expect(narrower.lay.treeX1 - narrower.lay.treeX0).toBe(0);
    expect(lay.treeX1 - lay.treeX0).toBe(layout.narrow!.tree ? 120 : 0);
    // the least width keeps the plot at its floor with the columns at theirs
    expect(layout.narrow!.minWidth).toBe(columns - (112 - 80) - 120 + 160);
    const least = await laid(Math.ceil(layout.narrow!.minWidth));
    expect(least.box.width).toBeGreaterThanOrEqual(160 - 1);
  });

  it('a collapsed tree leaves its toggle in its place: a small tree glyph, hit as `tree-toggle`; opened (a band’s `tree`), whole', async () => {
    const wide = await laid(1280);
    const columns = wide.layout.narrow!.columnsPx;
    expect(wide.lay.toggleX).toBe(-1);
    const narrow = await laid(columns - 32 + 160 - 40);
    expect(narrow.layout.narrow!.tree).toBe(0);
    const { lay } = narrow;
    expect(lay.toggleX).toBeGreaterThanOrEqual(0);
    // the labels after it, never under it
    expect(lay.labelX).toBeGreaterThanOrEqual(lay.toggleX + 12);
    expect(hitView(narrow.layout, lay.toggleX + 6, lay.toggleY + 6)).toMatchObject({ target: 'tree-toggle' });
    // the reader opened it: the tree whole, no toggle, the least width holding it (the frame scrolls)
    const open = await laid(columns - 32 + 160 - 40, true);
    expect(open.layout.narrow!.tree).toBe(120);
    // its toggle stays at the band's top, pressed, the tree beside it: the way back on screen
    expect(open.lay.toggleX).toBeGreaterThanOrEqual(0);
    expect(open.lay.toggleOpen).toBe(true);
    expect(open.lay.treeX0).toBeGreaterThanOrEqual(open.lay.toggleX + 12);
    expect(open.lay.treeX1 - open.lay.treeX0).toBeGreaterThan(100);
    expect(hitView(open.layout, open.lay.toggleX + 6, open.lay.toggleY + 6)).toMatchObject({ target: 'tree-toggle' });
    expect(narrow.lay.toggleOpen).toBe(false);
    expect(open.layout.narrow!.minWidth).toBe(narrow.layout.narrow!.minWidth + 120);
  });
});
