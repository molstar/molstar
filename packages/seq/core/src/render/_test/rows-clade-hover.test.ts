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
const Rows = '/data/aln/(rows)';
const T = '/views/v/s/rows';

function spec(params: Record<string, unknown> = {}) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick: Newick });
  b.views.view('v', (v) =>
    v.section('s', { of: aln } as never, (s) =>
      s.track('rows', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln, tree, ...params } as never)),
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

async function hovered(params: Record<string, unknown>, band: Record<string, unknown>, key: string, target: string) {
  const s = await settle(spec(params).build());
  const lay = rowsLay(await s.frame(at(band)));
  const hover = {
    layer: lay.key,
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
    const ys = bands.map((r) => r[1]).sort((a, b) => a - b);
    expect(ys).toContain(lay.rowY(lay.display.displayRow(4)));
    expect(ys).toContain(lay.rowY(lay.display.placeRow(3)));
  });
});
