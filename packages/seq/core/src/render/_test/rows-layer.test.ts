/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rows layer's places and marks (_spec/README.md "Render", The rows layer): a band's top in display rows and a
// kept window by its top row's key; a pinned row's place in the order (its tree leaf, its hover and selection there
// too); the pinned rows drawn on their own over image tiles offset from the pinned band; the selection's cost bounded
// by the rows in view and a summary row's partial mark; the tree's geometry kept a pan and a branch hit by its
// segment; one scratch hit; the measures merged over the theme's; the header (the tree's own width, labels cut with an
// ellipsis, the last numbers in a right gutter); the text figure only in text; one rows layer a section; a lazy
// alignment's pins said, its sub-pixel band a neutral fill.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it, vi } from 'vitest';
import { type Frame, type FrameInput, hitView, Renderer, recordDrawLists } from '../index.js';
import type { RowsLay } from '../rows.js';
import { settle } from './helpers.js';
import { rowsSpec, Section, Synthetic, Track } from './rows-fixtures.js';

const Labels = ['HBA/1-6', 'HBB/2-7', 'MYG/1-6', 'GLB/3-8', 'CYG/1-6'];
const Letters = ['MVLSPA', 'MVHLTP', 'MGLS--', 'MKAVLG', 'MEKVPG'];
const Newick = '((HBA:1,HBB:1)hb:1,(MYG:1,(GLB:1,CYG:1)gc:1)rest:1)root;';
const Rows = '/data/aln/(rows)';
const T = '/views/v/s/rows';

function spec(params: Record<string, unknown> = {}, track: Record<string, unknown> = {}, section = {}) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
  const tree = b.data.tree('tree', { of: aln, newick: Newick });
  b.views.view('v', (v) =>
    v.section('s', { of: aln, ...section } as never, (s) =>
      s.track('rows', { label: '', ...track }, (t) =>
        t.layer('r', { type: 'rows', data: aln, tree, ...params } as never),
      ),
    ),
  );
  return b;
}

const Size = { width: 600, height: 300 };

function rowsLays(f: Frame): RowsLay[] {
  const out: RowsLay[] = [];
  for (const s of f.layout.sections)
    for (const l of s.lays) for (const x of l ?? []) if (x.type.name === 'rows') out.push(x.lay as RowsLay);
  return out;
}
const rowsLay = (f: Frame) => rowsLays(f)[0];

const at = (band: Record<string, unknown>, extra: Partial<FrameInput> = {}): FrameInput => ({
  size: Size,
  dpr: 1,
  settled: true,
  viewports: { '/views/v/s': { start: 0, end: 6, bands: { [T]: { bandTop: 0, ...band } as never } } },
  ...extra,
});

/** The rects of a pass as `[x, y, w, h, fill, stroke]` numbers. */
function rects(f: Frame, pass: 'hover' | 'selection' | 'data') {
  return recordDrawLists(f.lists, [pass])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => l.slice(5).split(' '));
}

describe('a band’s top in display rows and a kept window by key', () => {
  it('`bandTop` is display rows below the pinned band, the layout’s px its times the row height, clamped', async () => {
    const s = await settle(spec({}, { height: 40 }).build());
    const lay = rowsLay(await s.frame(at({ bandTop: 1.5, rowScale: 20 })));
    expect([lay.top, lay.bandTop, lay.first]).toEqual([1.5, 30, 1]);
    // the last rows fill the band: 5 rows of 20 px in 40 px, at most 3 rows down
    expect(rowsLay(await s.frame(at({ bandTop: 99, rowScale: 20 }))).top).toBe(3);
  });

  it('a key at the top and back: a folded row at its summary row, a hidden one at the nearest row shown', async () => {
    const s = await settle(spec({ order: 'tree' }, { height: 40 }).build());
    const lay = rowsLay(await s.frame(at({ bandTop: 1.25, rowScale: 20 })));
    expect(lay.keyAt(lay.top)).toEqual({ key: 'HBB/2-7', offset: 0.25 });
    expect(lay.topOf('HBB/2-7', 0.25)).toBe(1.25);
    // under the fold of `hb`, HBB is in its summary row (display row 0)
    expect(lay.topOf('HBB/2-7', 0.25, ['hb'])).toBe(0.25);
    // a fold above the window keeps the window's top row: GLB is row 3, then row 2
    expect(lay.topOf('GLB/3-8', 0, ['hb'])).toBe(2);
    expect(lay.topOf('nobody', 0)).toBeUndefined();
    const h = await settle(spec({ order: 'tree', hidden_rows: ['MYG/1-6'] }, { height: 40 }).build());
    const lh = rowsLay(await h.frame(at({ bandTop: 0, rowScale: 20 })));
    // MYG hidden: the nearest row shown after it in the order, GLB
    expect(lh.topOf('MYG/1-6', 0)).toBe(lh.display.placeRow(3));
    expect(lh.display.nearestRow(2)).toBe(2);
  });
});

describe('a pinned row’s place in the order', () => {
  it('its tree leaf in place, its hover band and its selection tint both in the band and in place', async () => {
    const s = await settle(spec({ order: 'tree', pinned: ['CYG/1-6'] }).build());
    const f = await s.frame(at({ rowScale: 20 }));
    const lay = rowsLay(f);
    const band = lay.display.displayRow(4);
    const place = lay.display.placeRow(4);
    expect([band, place]).toEqual([0, 5]);
    // the tree's leaf of CYG at its place in the order: a branch at its row's middle
    const mid = lay.rowY(place) + lay.rowScale / 2;
    const ys = [];
    for (let k = 0; k + 3 < lay.tree.length; k += 4) if (lay.tree[k + 1] === lay.tree[k + 3]) ys.push(lay.tree[k + 1]);
    expect(ys).toContain(mid);
    expect(ys).not.toContain(lay.rowY(band) + lay.rowScale / 2);
    // the hover of CYG: its band at both places
    const hover = { layer: lay.key, section: 's', system: Rows, position: 4, target: 'row' as const };
    const fh = s.renderer.frame({ ...at({ rowScale: 20 }), emphasis: { hover } });
    expect(
      rects(fh, 'hover')
        .map((r) => Number(r[1]))
        .sort((a, b) => a - b),
    ).toEqual([lay.rowY(band), lay.rowY(place)]);
    // the selection of CYG: tinted at both places
    const selection = rangesLocation(Rows, [4, 5]);
    const fs = s.renderer.frame({ ...at({ rowScale: 20 }), emphasis: { selection } });
    expect(
      rects(fs, 'selection')
        .map((r) => Number(r[1]))
        .sort((a, b) => a - b),
    ).toEqual([lay.rowY(band), lay.rowY(place)]);
  });

  it('`pinInPlace: false`: the row in the band only, the tree and the hover there', async () => {
    const s = await settle(spec({ pinned: ['CYG/1-6'] }).build(), [], { rows: { pinInPlace: false } });
    const lay = rowsLay(await s.frame(at({ rowScale: 20 })));
    expect(lay.display.count).toBe(5);
    expect([lay.display.displayRow(4), lay.display.placeRow(4)]).toEqual([0, 0]);
    expect(lay.labelText).toEqual(['CYG/1-6', 'HBA/1-6', 'HBB/2-7', 'MYG/1-6', 'GLB/3-8']);
  });

  it('`rowAtY` in the pinned band is the pinned row, below it the scrolled rows', async () => {
    const s = await settle(spec({ pinned: ['CYG/1-6'] }, { height: 60 }).build());
    const lay = rowsLay(await s.frame(at({ bandTop: 2, rowScale: 20 })));
    expect(lay.rowAtY(lay.box.y + 5)).toBe(0);
    // below the pinned band: display row 1 + the top's 2 rows
    expect(lay.rowAtY(lay.box.y + 25)).toBe(3);
  });
});

describe('the pinned rows over image tiles', () => {
  it('drawn on their own, the tiles offset from the pinned band and spanning the scrolled rows only', async () => {
    const s = await settle(spec({ pinned: ['CYG/1-6'] }, { height: 60 }).build());
    const f = await s.frame(at({ bandTop: 0, rowScale: 2 }));
    const lay = rowsLay(f);
    expect(lay.images.length).toBe(1);
    const [, y0, , h] = lay.ib;
    // the tiles start below the pinned row and hold the 5 scrolled rows (the pinned one in place among them)
    expect(y0).toBe(lay.box.y + lay.pinnedHeight);
    expect(h).toBe(5 * lay.rowScale);
    // the pinned row's cells are rects of their own, past `pinnedFrom`
    expect(lay.n - lay.pinnedFrom).toBeGreaterThan(0);
    for (let k = lay.pinnedFrom; k < lay.n; k++) expect(lay.y[k]).toBe(lay.box.y);
  });
});

describe('the selection', () => {
  it('reads the rows in view only, whatever its size (200,000 lazy rows selected)', async () => {
    const s = await settle(rowsSpec(200_000, false), [Synthetic]);
    const input: FrameInput = {
      size: { width: 900, height: 400 },
      dpr: 1,
      settled: true,
      viewports: { [Section]: { start: 0, end: 60, bands: { [Track]: { bandTop: 100_000, rowScale: 20 } } } },
    };
    const lay = rowsLay(await s.frame(input));
    const shown = lay.last - lay.first;
    const rowAt = vi.spyOn(lay.display, 'rowAt');
    const displayRow = vi.spyOn(lay.display, 'displayRow');
    const f = s.renderer.frame({ ...input, emphasis: { selection: rangesLocation(lay.rowsSystem, [0, 200_000]) } });
    expect(rects(f, 'selection').length).toBe(shown);
    expect(rowAt.mock.calls.length + displayRow.mock.calls.length).toBeLessThanOrEqual(2 * shown);
  });

  it('a summary row draws the partial mark: the tint over the share selected, the outline when all are', async () => {
    const s = await settle(spec({ order: 'tree', collapsed: ['hb'] }).build());
    const lay = rowsLay(await s.frame(at({ rowScale: 20, folds: ['hb'] })));
    const sel = (r: number[]) =>
      rects(
        s.renderer.frame({ ...at({ rowScale: 20, folds: ['hb'] }), emphasis: { selection: rangesLocation(Rows, r) } }),
        'selection',
      );
    const half = sel([0, 1]);
    expect(half.length).toBe(1);
    expect(Number(half[0][2])).toBeCloseTo(lay.box.width / 2, 6);
    expect(half[0].join(' ')).not.toContain('stroke');
    const all = sel([0, 2]);
    expect(Number(all[0][2])).toBeCloseTo(lay.box.width, 6);
  });
});

describe('the tree and the hits', () => {
  it('keeps its geometry a pan (nothing sized by the nodes made), and hits a branch by its segment', async () => {
    const s = await settle(spec({ order: 'tree' }, { height: 60 }).build());
    const l0 = rowsLay(await s.frame(at({ bandTop: 0, rowScale: 20 })));
    const coord = l0.coord;
    const l1 = rowsLay(await s.frame(at({ bandTop: 1, rowScale: 20 })));
    expect(l1.coord).toBe(coord);
    // a point on `gc`'s branch (from `rest`'s x to `gc`'s, at gc's y), far from both nodes' points
    const f = await (await settle(spec({ order: 'tree' }).build())).frame(at({ bandTop: 0, rowScale: 20 }));
    const lay = rowsLay(f);
    const gc = lay.clades!.byName.get('gc')!;
    const s0 = lay.treeNode.indexOf(gc);
    const [x0, y, x1] = [lay.tree[4 * s0], lay.tree[4 * s0 + 1], lay.tree[4 * s0 + 2]];
    expect(Math.abs(x1 - x0)).toBeGreaterThan(20);
    const h = hitView(f.layout, (x0 + x1) / 2, y + 4);
    expect(h).toMatchObject({ target: 'tree-node', element: { key: 'gc' } });
    // a leaf's branch hits no node: the row's label there is past the tree
    expect(lay.treeNode.filter((n) => n >= 0).every((n) => lay.clades!.hi[n] - lay.clades!.lo[n] > 1)).toBe(true);
  });

  it('one scratch hit: a hover allocates no hit object', async () => {
    const s = await settle(spec().build());
    const f = await s.frame(at({ rowScale: 20 }));
    const lay = rowsLay(f);
    lay.axis.cell(1);
    const x = (lay.axis.sx0 + lay.axis.sx1) / 2;
    const h1 = hitView(f.layout, x, lay.rowY(0) + 5);
    const h2 = hitView(f.layout, x, lay.rowY(2) + 5);
    expect(h1).toBe(lay.scratch);
    expect(h2).toBe(lay.scratch);
    expect(h2).toMatchObject({ target: 'cell', position: 1, yPosition: 2 });
  });
});

describe('the colors and the header', () => {
  it('the rows’ measures are merged over the theme’s, never replacing them', async () => {
    const custom = { name: 'custom', read: () => 1 };
    const s = await settle(spec().build());
    await s.frame(at({ rowScale: 20 }));
    const r = new Renderer({
      ...s.renderer.options,
      theme: { measures: new Map([['custom', custom]]) },
    });
    const lay = rowsLay(r.frame(at({ rowScale: 20 })));
    expect([...lay.measureEnv!.measures.keys()].sort()).toEqual(['custom', 'identity', 'matches-reference', 'share']);
  });

  it('the label column widened by the tree’s width; a label past its room cut with an ellipsis', async () => {
    const plain = rowsLay(await (await settle(spec().build())).frame(at({ rowScale: 20 })));
    const tree = rowsLay(
      await (await settle(spec({ order: 'tree', tree_width: 120 }).build())).frame(at({ rowScale: 20 })),
    );
    expect(tree.box.x - plain.box.x).toBe(120);
    expect(tree.treeX1 - tree.treeX0).toBe(120);
    expect(tree.labelX).toBe(plain.labelX + 124);
    expect(tree.labelText).toContain('HBA/1-6');
    const narrow = rowsLay(
      await (await settle(spec({}, {}).build(), [], { labelWidth: 20 })).frame(at({ rowScale: 20 })),
    );
    expect(narrow.labelText.some((t) => t.endsWith('…'))).toBe(true);
    for (const t of narrow.labelText) expect(Labels.includes(t) || t.endsWith('…')).toBe(true);
    for (let k = 0; k < narrow.labelW.length; k++) expect(narrow.labelW[k]).toBeLessThanOrEqual(narrow.labelMax);
  });

  it('the last residue numbers in a right gutter, none for a row of gaps', async () => {
    const lay = rowsLay(await (await settle(spec().build())).frame(at({ rowScale: 20 })));
    expect(lay.numberText).toEqual(['1', '2', '1', '3', '1']);
    // MYG has 4 letters from 1: its last is 4
    expect(lay.lastText).toEqual(['6', '7', '4', '8', '6']);
    expect(lay.numberRightX).toBeGreaterThan(lay.box.x + lay.box.width);
    expect(lay.box.x + lay.box.width).toBeLessThan(Size.width - 30);
  });
});

describe('where the rows draw', () => {
  it('a figure only in text; in another layout not drawn, said', async () => {
    const s = await settle(spec({}, {}, { layout: 'blocks' }).build());
    const f = await s.frame({ size: Size, dpr: 1, settled: true });
    expect(f.layout.issues.some((i) => i.code === 'seq.render.rows.layout')).toBe(true);
    for (const l of rowsLays(f)) expect(l.n + l.images.length).toBe(0);
  });

  it('a text figure’s rows share its band: the display’s count (hidden rows out), not the alignment’s', async () => {
    const b = spec({ hidden_rows: ['GLB/3-8'] }, {}, { layout: 'text', group: 3, max_width: 48 });
    const s = await settle(b.build());
    const lays = rowsLays(await s.frame({ size: Size, dpr: 1, settled: true }));
    expect(lays.length).toBeGreaterThan(0);
    for (const l of lays) {
      expect(l.figure).toBe(true);
      expect(l.display.count).toBe(4);
      expect(l.rowScale).toBeCloseTo(l.box.height / 4, 6);
    }
  });

  it('one rows layer a section: a second one said and not drawn', async () => {
    const b = createSeqSpec();
    const aln = b.data.alignment('aln', { alphabet: 'protein', rows: { label: Labels, letters: Letters } });
    b.views.view('v', (v) =>
      v.section('s', { of: aln }, (s) => {
        s.track('a', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln } as never));
        s.track('b', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln } as never));
      }),
    );
    const f = await (await settle(b.build())).frame({ size: Size, dpr: 1, settled: true });
    expect(rowsLays(f).length).toBe(1);
    expect(f.layout.issues.filter((i) => i.code === 'seq.render.rows.second').map((i) => i.address)).toEqual([
      '/views/v/s/b/r',
    ]);
  });

  it('a lazy alignment: its pins and hidden rows said, its band below a pixel a neutral fill', async () => {
    const s = await settle(rowsSpec(200_000, false, { pinned: ['row 3'] }), [Synthetic]);
    const input: FrameInput = {
      size: { width: 900, height: 400 },
      dpr: 1,
      settled: true,
      viewports: { [Section]: { start: 0, end: 60, bands: { [Track]: { bandTop: 0, rowScale: 0.001 } } } },
    };
    const f = await s.frame(input);
    const lay = rowsLay(f);
    expect(f.layout.issues.some((i) => i.code === 'seq.render.rows.lazy')).toBe(true);
    expect(lay.neutral).toBe(true);
    expect(lay.pending.length).toBe(0);
  });
});

describe('the order and the cells', () => {
  it('a field order: by `first`, descending, ties in data order', async () => {
    const lay = rowsLay(
      await (await settle(spec({ order: { field: 'first', direction: 'descending' } }).build())).frame(
        at({ rowScale: 20 }),
      ),
    );
    // firsts 1, 2, 1, 3, 1
    expect(lay.labelText).toEqual(['GLB/3-8', 'HBB/2-7', 'HBA/1-6', 'MYG/1-6', 'CYG/1-6']);
  });

  it('a summary row’s cells: the consensus opaque by its share; `bar` the non-gap share, no letters', async () => {
    const s = await settle(spec({ order: 'tree' }).build());
    const lay = rowsLay(await s.frame(at({ rowScale: 20, folds: ['hb'] })));
    // hb: MVLSPA and MVHLTP; columns 0 and 1 agree (whole), column 2 half (L or H)
    const y = lay.rowY(0);
    const alphas: number[] = [];
    for (let k = 0; k < lay.n; k++) if (lay.y[k] === y) alphas.push(lay.fill[k] & 0xff);
    expect(alphas.length).toBe(6);
    const full = alphas[0];
    expect(alphas[1]).toBe(full);
    expect(alphas[2]).toBe(Math.round((full * 128) / 255));
    const b = await settle(spec({ order: 'tree', summary: 'bar' }).build());
    const lb = rowsLay(await b.frame(at({ rowScale: 20, folds: ['rest'] })));
    // rest: MGLS--, MKAVLG, MEKVPG; column 4 is 2 of 3 letters, column 0 all
    const sum = lb.display.rowAt(2) < 0 ? 2 : -1;
    expect(sum).toBe(2);
    const cells: number[][] = [];
    for (let k = 0; k < lb.n; k++) if (lb.y[k] === lb.rowY(sum)) cells.push([lb.code[k], lb.fill[k] & 0xff]);
    expect(cells.every((c) => c[0] === 0)).toBe(true);
    expect(cells[4][1]).toBe(Math.round((cells[0][1] * 170) / 255));
  });

  it('at the cells step a run of one fill on a row is one rect', async () => {
    const b = createSeqSpec();
    const aln = b.data.alignment('aln', {
      alphabet: 'protein',
      rows: { label: ['a', 'b'], letters: ['AAAAAA', 'AWAWAW'] },
    });
    b.views.view('v', (v) =>
      v.section('s', { of: aln }, (s) =>
        s.track('rows', { label: '' }, (t) => t.layer('r', { type: 'rows', data: aln } as never)),
      ),
    );
    const lay = rowsLay(await (await settle(b.build())).frame(at({ rowScale: 5 })));
    expect(lay.letters).toBe(false);
    const per = (i: number) => Array.from(lay.y.subarray(0, lay.n)).filter((y) => y === lay.rowY(i)).length;
    expect([per(0), per(1)]).toEqual([1, 6]);
  });
});
