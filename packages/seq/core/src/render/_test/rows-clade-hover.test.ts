/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A hovered clade (_spec/README.md "Render", The rows layer, Hover): a tree node's or a summary row's hit lights the
// clade's branches in the hover ink and one band over its rows; another clade's branches and rows stay dark.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type Frame, type FrameInput, recordDrawLists } from '../index.js';
import type { RowsLay } from '../rows.js';
import { settle } from './helpers.js';

const Labels = ['HBA/1-6', 'HBB/2-7', 'MYG/1-6', 'GLB/3-8', 'CYG/1-6'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLS--', 'MKAVLG', 'MEKVPG'];
const Newick = '((HBA:1,HBB:1)hb:1,(MYG:1,(GLB:1,CYG:1)gc:1)rest:1)root;';
/** Support values as inner labels: both clades `95`, so each is named by its first and last leaves. */
const Support = '((HBA:1,HBB:1)95:1,(MYG:1,(GLB:1,CYG:1)95:1)rest:1)root;';
const Chain = '(((HBA:1,HBB:1)x:1)y:1,(MYG:1,(GLB:1,CYG:1)gc:1)rest:1)root;';
const Rows = '/data/aln/(rows)';
const T = '/views/v/s/rows';

function spec(params: Record<string, unknown> = {}, newick = Newick, height?: number) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick });
  b.views.view('v', (v) =>
    v.section('s', { of: aln } as never, (s) =>
      s.track('rows', { label: '', ...(height ? { height } : {}) }, (t) =>
        t.layer('r', { type: 'rows', data: aln, tree, ...params } as never),
      ),
    ),
  );
  return b;
}

const at = (band: Record<string, unknown>, extra: Partial<FrameInput> = {}): FrameInput => ({
  size: { width: 600, height: 300 },
  dpr: 1,
  settled: true,
  viewports: { '/views/v/s': { start: 0, end: 6, bands: { [T]: { bandTop: 0, ...band } as never } } },
  ...extra,
});

function rowsLay(f: Frame): RowsLay {
  for (const s of f.layout.sections)
    for (const l of s.lays) for (const x of l ?? []) if (x.type.name === 'rows') return x.lay as RowsLay;
  throw new Error('no rows layer');
}

/** The hover pass's rects `[x, y, w, h, fill]`. */
function rects(f: Frame): number[][] {
  return recordDrawLists(f.lists, ['hover'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) =>
      l
        .slice(5)
        .split(' ')
        .map((v) => Number(v)),
    );
}

interface Options {
  readonly newick?: string;
  readonly height?: number;
  readonly pinInPlace?: boolean;
  readonly layer?: string;
}

async function hovered(
  params: Record<string, unknown>,
  band: Record<string, unknown>,
  key: string,
  target: string,
  o: Options = {},
) {
  const s = await settle(
    spec(params, o.newick, o.height).build(),
    [],
    o.pinInPlace === undefined ? {} : { rows: { pinInPlace: o.pinInPlace } },
  );
  const lay = rowsLay(await s.frame(at(band)));
  const hover = {
    layer: o.layer ?? lay.key,
    section: 's',
    system: Rows,
    position: undefined,
    target: target as 'tree-node',
    element: { owner: lay.owner, key },
  };
  const f = s.renderer.frame({ ...at(band), emphasis: { hover } });
  const all = rects(f);
  const bands = all.filter((r) => r[0] === lay.box.x && r[2] === lay.box.width);
  const branches = all.filter((r) => r[0] < lay.box.x);
  return { lay, bands, branches };
}

describe('a hovered clade', () => {
  it('a tree node: its branches and one band over exactly its rows, nothing for another clade', async () => {
    const { lay, bands, branches } = await hovered({ order: 'tree' }, { rowScale: 20 }, 'gc', 'tree-node');
    expect(bands.map((r) => [r[1], r[3]])).toEqual([[lay.rowY(3), 2 * lay.rowScale]]);
    // gc's own stem, its two leaf branches and its vertical: four segments, all left of the data
    expect(branches).toHaveLength(4);
    // none above gc's rows (the HBA, HBB and MYG branches are not lit)
    for (const r of branches) expect(r[1] + r[3] / 2).toBeGreaterThanOrEqual(lay.rowY(2) + lay.rowScale / 2 - 1);
    const other = await hovered({ order: 'tree' }, { rowScale: 20 }, 'hb', 'tree-node');
    expect(other.bands.map((r) => [r[1], r[3]])).toEqual([[other.lay.rowY(0), 2 * other.lay.rowScale]]);
  });

  it('a summary row: the folded clade’s row band and its wedge', async () => {
    const { lay, bands, branches } = await hovered(
      { order: 'tree' },
      { rowScale: 20, folds: ['hb'] },
      'hb',
      'summary-row',
    );
    expect(bands.map((r) => [r[1], r[3]])).toEqual([[lay.rowY(0), lay.rowScale]]);
    expect(branches.length).toBeGreaterThan(0);
  });

  it('a pinned row of the clade: lit in the pinned band and in its place', async () => {
    const { lay, bands } = await hovered({ order: 'tree', pinned: ['CYG/1-6'] }, { rowScale: 20 }, 'gc', 'tree-node');
    const rs = lay.rowScale;
    // the pinned slot, and one run over GLB and CYG in the order
    expect(bands.map((r) => [r[1], r[3]]).sort((a, b) => a[0] - b[0])).toEqual([
      [lay.rowY(lay.display.displayRow(4)), rs],
      [lay.rowY(lay.display.placeRow(3)), 2 * rs],
    ]);
    expect(lay.display.placeRow(4)).toBe(lay.display.placeRow(3) + 1);
  });

  it('`pinInPlace: false`: a pinned row of the clade is tinted once, never joined to the scrolled rows', async () => {
    const { lay, bands } = await hovered(
      { order: 'tree', pinned: ['CYG/1-6'] },
      { rowScale: 20 },
      'rest',
      'tree-node',
      { pinInPlace: false },
    );
    const rs = lay.rowScale;
    expect(bands.map((r) => [r[1], r[3]]).sort((a, b) => a[0] - b[0])).toEqual([
      [lay.rowY(0), rs],
      [lay.rowY(lay.display.placeRow(2)), 2 * rs],
    ]);
  });

  it('a clade scrolled under the pinned band never covers it', async () => {
    const { lay, bands } = await hovered(
      { order: 'tree', pinned: ['CYG/1-6'] },
      { rowScale: 20, bandTop: 1.5 },
      'hb',
      'tree-node',
      { height: 60 },
    );
    expect(lay.pinnedHeight).toBeGreaterThan(0);
    expect(bands.length).toBeGreaterThan(0);
    for (const r of bands) expect(r[1]).toBeGreaterThanOrEqual(lay.box.y + lay.pinnedHeight);
  });

  it('a unary chain: hovering the inner node leaves its parent’s stem dark', async () => {
    const x = await hovered({ order: 'tree' }, { rowScale: 20 }, 'x', 'tree-node', { newick: Chain });
    const y = await hovered({ order: 'tree' }, { rowScale: 20 }, 'y', 'tree-node', { newick: Chain });
    expect(x.branches.length).toBeLessThan(y.branches.length);
  });

  it('two clades labelled alike (support values): each hover reaches its own rows, a fold its own clade', async () => {
    const o = { newick: Support };
    const gc = await hovered({ order: 'tree' }, { rowScale: 20 }, 'GLB..CYG', 'tree-node', o);
    expect(gc.bands.map((r) => [r[1], r[3]])).toEqual([[gc.lay.rowY(3), 2 * gc.lay.rowScale]]);
    const hb = await hovered({ order: 'tree' }, { rowScale: 20 }, 'HBA..HBB', 'tree-node', o);
    expect(hb.bands.map((r) => [r[1], r[3]])).toEqual([[hb.lay.rowY(0), 2 * hb.lay.rowScale]]);
    const folded = await hovered(
      { order: 'tree' },
      { rowScale: 20, folds: ['GLB..CYG'] },
      'GLB..CYG',
      'summary-row',
      o,
    );
    // HBA, HBB, MYG, then GLB and CYG in one summary row
    expect(folded.lay.display.count).toBe(4);
    expect(folded.bands.map((r) => [r[1], r[3]])).toEqual([[folded.lay.rowY(3), folded.lay.rowScale]]);
  });

  it('`collapsed` naming no clade (a support value): said once, an info, the rest folded', async () => {
    const s = await settle(spec({ order: 'tree', collapsed: ['95', 'GLB..CYG'] }, Support).build());
    const f = await s.frame(at({ rowScale: 20 }));
    const said = f.layout.issues.filter((i) => i.code === 'seq.render.rows.clade');
    expect(said).toHaveLength(1);
    expect(said[0].severity).toBe('info');
    expect(said[0].message).toContain('95');
    expect(said[0].message).not.toContain('GLB..CYG');
    // GLB and CYG folded into one row
    expect(rowsLay(f).display.count).toBe(4);
    const named = await settle(spec({ order: 'tree', collapsed: ['gc'] }).build());
    expect(
      (await named.frame(at({ rowScale: 20 }))).layout.issues.some((i) => i.code === 'seq.render.rows.clade'),
    ).toBe(false);
  });

  it('a hover from another layer lights no clade', async () => {
    const { bands, branches } = await hovered({ order: 'tree' }, { rowScale: 20 }, 'gc', 'tree-node', {
      layer: '/views/v/s/rows/other',
    });
    expect([bands.length, branches.length]).toEqual([0, 0]);
  });
});
