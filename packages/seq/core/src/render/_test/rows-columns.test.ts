/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rowsSystemOf } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { DarkTokens, hitView, recordDrawLists, renderSvg } from '../index.js';
import type { RowsLay } from '../rows.js';
import { fitRowColumns, rowColumnsOf } from '../rows-model.js';
import { settle } from './helpers.js';

// The rows layer's columns (_spec/README.md "Render", The rows layer, Row columns): the fields they read (the rows'
// own, the measures, a joined records source), how each draws on either side in the canvas and the SVG, their titles
// (a hit that orders), their yielding to a narrow view, and the text figure's.

const Labels = ['A/1-5', 'B/1-5', 'C/1-4'];
// identity to A over the columns where both have a residue: B 3 of 4, C 3 of 3; lengths 5, 5, 4; gaps 1, 1, 2
const Letters = ['MVLS-A', 'MV-LTA', 'M--SDA'];
const Species = ['Homo sapiens', 'Mus musculus', null];

type Col = { field: string; side?: string; as?: string; width?: number; title?: string };

const AllKinds: Col[] = [
  { field: 'OS', title: 'Species' },
  { field: 'OS', as: 'chip', title: 'Kind' },
  { field: 'identity', side: 'right', as: 'number' },
  { field: 'length', side: 'right', as: 'bar' },
];

function spec(columns: Col[], params: Record<string, unknown> = {}, records?: string) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: Labels, letters: Letters, OS: Species },
    reference: 'A/1-5',
  } as never);
  if (records !== undefined)
    b.data.source('extra', {
      format: 'tsv',
      value: records,
      system: rowsSystemOf('/data/aln'),
      types: { score: 'int32' },
    } as never);
  b.views.view('v', (v) =>
    v.section('s', { of: aln, ...(params.section as object) }, (s) => {
      s.track('ruler', { label: 'Columns' }, (t) => t.layer('ruler', { type: 'ruler' }));
      s.track('rows', { label: '' }, (t) =>
        t.layer('r', { type: 'rows', data: aln, columns, height: 3, ...(params.layer as object) } as never),
      );
    }),
  );
  return b;
}

const Size = { width: 900, height: 300 };

function rowsLays(f: {
  layout: { sections: readonly { lays: readonly (readonly { lay: unknown; type: { name: string } }[])[] }[] };
}) {
  const out: RowsLay[] = [];
  for (const s of f.layout.sections)
    for (const l of s.lays) for (const x of l ?? []) if (x.type.name === 'rows') out.push(x.lay as RowsLay);
  return out;
}

async function laid(columns: Col[], params: Record<string, unknown> = {}, size = Size, records?: string) {
  const b = spec(columns, params, records);
  expect(b.check()).toEqual([]);
  const s = await settle(b.build(), records !== undefined ? [RecordsRuntime] : []);
  const f = await s.frame({ size, dpr: 1, settled: true });
  return { b, s, f, lay: rowsLays(f as never)[0] };
}

/** Each display row's text of a column, by its index in `columns`. */
const cellsOf = (lay: RowsLay, k: number) =>
  Array.from({ length: lay.display.count }, (_, i) => {
    const q = lay.cellRow.findIndex((r, j) => r === i && lay.cellCol[j] === k);
    return q < 0 ? '' : lay.cellText[q];
  });

describe('row columns', () => {
  it('read the rows’ own fields and the measures: identity, length and gaps as computed by hand', async () => {
    const { lay } = await laid(AllKinds);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'identity'))).toEqual([1, 0.75, 1]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'length'))).toEqual([5, 5, 4]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'gaps'))).toEqual([1, 1, 2]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'OS'))).toEqual(['Homo sapiens', 'Mus musculus', undefined]);
    // a missing value is an empty cell; the number reads as a percent
    expect(cellsOf(lay, 0)).toEqual(['Homo sapiens', 'Mus musculus', '']);
    expect(cellsOf(lay, 2)).toEqual(['100%', '75%', '100%']);
    // a bar is its share of the column's range (from 0): 5, 5 and 4 of 5
    const bars = lay.cellShare.filter((_, q) => lay.cellCol[q] === 3);
    expect(bars).toEqual([1, 1, 0.8]);
    // a chip's color is its category's (two categories, two colors)
    const chips = new Set(lay.cellColor.filter((_, q) => lay.cellCol[q] === 1));
    expect(chips.size).toBe(2);
  });

  it('draw left in the header and right after the last numbers, alike in the canvas and the SVG, in both themes', async () => {
    const { b, f, lay } = await laid(AllKinds);
    // left: between the label and the first numbers; right: past the data and its last numbers
    expect(lay.colX[0]).toBeGreaterThan(lay.labelX);
    expect(lay.colX[1] + lay.colW[1]).toBeLessThanOrEqual(lay.numberX);
    expect(lay.colX[2]).toBeGreaterThanOrEqual(lay.box.x + lay.box.width);
    expect(lay.colX[3]).toBeGreaterThan(lay.colX[2]);
    expect(Math.max(...lay.lastW.map((w) => lay.numberRightX + w))).toBeLessThanOrEqual(lay.colX[2]);
    for (const tokens of [undefined, DarkTokens] as const) {
      const s = await settle(b.build());
      if (tokens) s.renderer.tokens = tokens;
      const g = await s.frame({ size: Size, dpr: 1, settled: true });
      const svg = await renderSvg(b.build(), { size: Size, ...(tokens ? { tokens } : {}) });
      expect(svg.issues).toEqual([]);
      expect(recordDrawLists(svg.lists, ['data', 'labels'])).toBe(recordDrawLists(g.lists, ['data', 'labels']));
      for (const t of ['Species', 'Kind', 'identity', 'length', 'Homo sapiens', 'Mus musculus', '75%'])
        expect(svg.svg).toContain(t);
    }
    expect(f.layout.issues).toEqual([]);
  });

  it('hit as their row’s header, naming the value; a title hits as itself, its field the state', async () => {
    const { f, lay } = await laid(AllKinds);
    const y = lay.rowY(1) + lay.rowScale / 2;
    expect(hitView(f.layout, lay.colX[0] + 2, y)).toMatchObject({
      target: 'row',
      position: 1,
      state: 'Species: Mus musculus',
    });
    expect(hitView(f.layout, lay.colX[2] + lay.colW[2] - 2, y)).toMatchObject({
      target: 'row',
      position: 1,
      state: 'identity: 75%',
    });
    // the label's own hit names nothing
    expect(hitView(f.layout, lay.labelX + 2, y)?.state).toBeUndefined();
    const t = lay.titleCol.indexOf(2);
    expect(hitView(f.layout, lay.colX[2] + 2, lay.titleY - 3)).toMatchObject({
      target: 'column-title',
      state: 'identity',
      layer: '/views/v/s/rows/r',
    });
    expect(lay.titleText[t]).toBe('identity');
  });

  it('an order by a column’s field orders the rows, its title marked; a field no row has is said once', async () => {
    const { lay } = await laid(AllKinds, { layer: { order: { field: 'identity', direction: 'ascending' } } });
    expect(lay.labelText).toEqual(['B/1-5', 'A/1-5', 'C/1-4']);
    expect(lay.titleText[lay.titleCol.indexOf(2)]).toBe('identity ↑');
    const { lay: byOs } = await laid(AllKinds, { layer: { order: { field: 'OS', direction: 'descending' } } });
    // missing last in either direction
    expect(byOs.labelText).toEqual(['B/1-5', 'A/1-5', 'C/1-4']);
    const { f } = await laid([{ field: 'nothing' }, { field: 'OS' }]);
    const said = f.layout.issues.filter((i) => i.code === 'seq.render.rows.column');
    expect(said.length).toBe(1);
    expect(said[0].message).toContain('nothing');
  });

  it('join a records source whose system is the rows system, by row key (an id) or by position', async () => {
    // a score of 10 above 9: numbers by value, not by their text
    const tsv = 'id\tclade\tscore\nC/1-4\tfish\t10\nA/1-5\tprimate\t9\nZ\tnone\t1\n';
    const { lay, f } = await laid(
      [
        { field: 'clade', side: 'right' },
        { field: 'score', side: 'right', as: 'number' },
      ],
      { layer: { order: { field: 'score', direction: 'descending' } } },
      Size,
      tsv,
    );
    expect(f.layout.issues).toEqual([]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'clade'))).toEqual(['primate', undefined, 'fish']);
    // ordered by the joined field: B has none (last)
    expect(lay.labelText).toEqual(['C/1-4', 'A/1-5', 'B/1-5']);
    expect(cellsOf(lay, 1)).toEqual(['10', '9', '']);
    // by position: a record a row, from the first
    const { lay: byPosition } = await laid([{ field: 'note' }], {}, Size, 'position\tnote\n2\tsecond\n');
    expect([0, 1, 2].map((r) => byPosition.fieldValue(r, 'note'))).toEqual([undefined, 'second', undefined]);
  });

  it('yield to a narrow view: right ones first, each down to a floor and then hidden, a chip counting them', async () => {
    const cols = rowColumnsOf([
      { field: 'a', width: 80 },
      { field: 'b', width: 80 },
      { field: 'c', side: 'right', width: 80 },
      { field: 'd', side: 'right', width: 80 },
    ]);
    expect(fitRowColumns(cols, 0, 32, 6, [])).toEqual([80, 80, 80, 80]);
    expect(fitRowColumns(cols, 40, 32, 6, [])).toEqual([80, 80, 80, 40]);
    expect(fitRowColumns(cols, 60, 32, 6, [])).toEqual([80, 80, 80, 0]);
    // the right ones hidden (the last first), then the left ones shrink from the last
    expect(fitRowColumns(cols, 180, 32, 6, [])).toEqual([80, 72, 0, 0]);
    expect(fitRowColumns(cols, 250, 32, 6, [])).toEqual([80, 0, 0, 0]);
    expect(fitRowColumns(cols, 1000, 32, 6, [])).toEqual([0, 0, 0, 0]);
    // laid out: the wide view shows all, a narrower one hides the right ones first, the narrowest all of them
    const widths = async (width: number) => {
      const { lay } = await laid(AllKinds, {}, { width, height: 300 });
      return { w: [...lay.colW], chip: lay.colChip, x: lay.box.x, data: lay.box.width };
    };
    const wide = await widths(1000);
    expect(wide.w.every((w) => w > 0)).toBe(true);
    expect(wide.chip).toBe('');
    const mid = await widths(640);
    expect(mid.w[3]).toBe(0);
    expect(mid.w[0]).toBeGreaterThan(0);
    expect(mid.chip).toMatch(/columns? hidden$/);
    const narrow = await widths(400);
    expect(narrow.w).toEqual([0, 0, 0, 0]);
    expect(narrow.chip).toBe('4 columns hidden');
    // the label keeps its room: the data starts no later than without columns
    const { lay: none } = await laid([], {}, { width: 400, height: 300 });
    expect(narrow.x).toBe(none.box.x);
  });

  it('a text figure draws the same columns on every line, their titles over the first', async () => {
    const { f } = await laid(AllKinds, { section: { layout: 'text', group: 3, max_width: 48 } });
    const lays = rowsLays(f as never);
    expect(lays.length).toBeGreaterThan(1);
    for (const lay of lays) {
      expect(lay.figure).toBe(true);
      expect(cellsOf(lay, 0)).toEqual(['Homo sapiens', 'Mus musculus', '']);
      expect(cellsOf(lay, 2)).toEqual(['100%', '75%', '100%']);
    }
    expect(lays[0].titleText.length).toBe(4);
    expect(lays[1].titleText.length).toBe(0);
    const svg = await renderSvg(spec(AllKinds, { section: { layout: 'text', group: 3, max_width: 48 } }).build(), {
      size: Size,
    });
    expect(recordDrawLists(svg.lists, ['labels'])).toBe(recordDrawLists(f.lists, ['labels']));
  });
});
