/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { openAlignment } from '@molstar/seqx-core/data';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rowsSystemOf } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import { DarkTokens, hitView, Passes, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import type { RowsLay } from '../rows.js';
import { fitRowColumns, numericValues, rowColumnsOf } from '../rows-model.js';
import { settle } from './helpers.js';

// The rows layer's columns (_spec/README.md "Render", The rows layer, Row columns): the fields they read (the rows'
// own, the measures, a joined records source), how each draws on either side in the canvas and the SVG, their titles
// (a hit that orders) in the layer's title row, their yielding to a narrow view, and the text figure's.

const Labels = ['A/1-5', 'B/1-5', 'C/1-4'];
// identity to A over the columns where either has a residue (a gap against a residue a mismatch, gap against gap
// left out): B 3 of 6, C 3 of 6; lengths 5, 5, 4; gaps 1, 1, 2
const Letters = ['MVLS-A', 'MV-LTA', 'M--SDA'];
const Species = ['Homo sapiens', 'Mus musculus', null];
// a `#=GS` TAXID: text that reads as numbers (by text "10090" sorts before "9606")
const Taxid = ['9606', '10090', '7955'];

type Col = { field: string; side?: string; as?: string; width?: number; title?: string };

const AllKinds: Col[] = [
  { field: 'OS', title: 'Species' },
  { field: 'OS', as: 'chip', title: 'Kind' },
  { field: 'identity', side: 'right', as: 'number' },
  { field: 'length', side: 'right', as: 'bar' },
];

interface SpecOptions {
  readonly section?: Record<string, unknown>;
  readonly layer?: Record<string, unknown>;
  readonly records?: string;
  /** Another records source, a second joined one. */
  readonly records2?: string;
  readonly rows?: Record<string, unknown>;
  /** The rows track first in its section (no ruler above). */
  readonly first?: boolean;
}

function spec(columns: Col[], o: SpecOptions = {}) {
  const b = createSeqSpec();
  const aln = b.data.alignment('aln', {
    alphabet: 'protein',
    rows: { label: Labels, letters: Letters, OS: Species, TAXID: Taxid, ...o.rows },
    reference: 'A/1-5',
  } as never);
  for (const [id, value] of [
    ['extra', o.records],
    ['more', o.records2],
  ] as const)
    if (value !== undefined)
      b.data.source(id, {
        format: 'tsv',
        value,
        system: rowsSystemOf('/data/aln'),
        types: { score: 'int32' },
      } as never);
  b.views.view('v', (v) =>
    v.section('s', { of: aln, ...o.section }, (s) => {
      if (!o.first) s.track('ruler', { label: 'Columns' }, (t) => t.layer('ruler', { type: 'ruler' }));
      s.track('rows', { label: '' }, (t) =>
        t.layer('r', { type: 'rows', data: aln, columns, height: 3, ...o.layer } as never),
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

async function laid(columns: Col[], o: SpecOptions = {}, size = Size) {
  const b = spec(columns, o);
  const s = await settle(b.build(), o.records !== undefined ? [RecordsRuntime] : []);
  const f = await s.frame({ size, dpr: 1, settled: true });
  return { b, s, f, lay: rowsLays(f as never)[0] };
}

/** Each display row's cut text of a column, by its index in `columns`. */
const cellsOf = (lay: RowsLay, k: number) =>
  Array.from({ length: lay.display.count }, (_, i) => {
    const c = lay.cols;
    for (let q = 0; q < c.n; q++) if (c.display[q] === i && c.col[q] === k) return c.cellText(q);
    return '';
  });

const shares = (lay: RowsLay, k: number) =>
  Array.from({ length: lay.cols.n }, (_, q) => q)
    .filter((q) => lay.cols.col[q] === k)
    .map((q) => lay.cols.share[q]);

const titles = (lay: RowsLay) => lay.cols.titleText.slice(0, lay.cols.titleN);

describe('row columns', () => {
  it('read the rows’ own fields and the measures: identity, length and gaps as computed by hand', async () => {
    const { b, lay } = await laid(AllKinds);
    expect(b.check()).toEqual([]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'identity'))).toEqual([1, 0.5, 0.5]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'length'))).toEqual([5, 5, 4]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'gaps'))).toEqual([1, 1, 2]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'OS'))).toEqual(['Homo sapiens', 'Mus musculus', undefined]);
    // a missing value is an empty cell; the measure's identity reads as a percent
    expect(cellsOf(lay, 0)).toEqual(['Homo sapiens', 'Mus musculus', '']);
    expect(cellsOf(lay, 2)).toEqual(['100%', '50%', '50%']);
    // a bar is its share of the column's range (from 0): 5, 5 and 4 of 5
    expect(shares(lay, 3)).toEqual([1, 1, 0.8]);
    // a chip's color is its category's (two categories, two colors)
    expect(
      new Set(
        Array.from({ length: lay.cols.n }, (_, q) => q)
          .filter((q) => lay.cols.col[q] === 1)
          .map((q) => lay.cols.color[q]),
      ).size,
    ).toBe(2);
  });

  it('identity leaves out the inserts (lowercase, `.`) on either row; `gaps` counts a row’s `-` only', async () => {
    const b = createSeqSpec();
    // column 2 an insert in B (`d`) and C (`.`): out of both identities; B's `.` not a gap
    const aln = b.data.alignment('aln', {
      alphabet: 'protein',
      rows: { label: ['A', 'B', 'C'], letters: ['MV-SA', 'MVd-A', 'MV.SA'] },
      reference: 'A',
    });
    b.views.view('v', (v) =>
      v.section('s', { of: aln }, (s) =>
        s.track('rows', { label: '' }, (t) =>
          t.layer('r', { type: 'rows', data: aln, columns: [{ field: 'identity' }] } as never),
        ),
      ),
    );
    const s = await settle(b.build());
    const lay = rowsLays((await s.frame({ size: Size, dpr: 1, settled: true })) as never)[0];
    // B: M V (d out) -/S mismatch, A → 3 of 4; C: M V (. out) S A → 4 of 4
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'identity'))).toEqual([1, 0.75, 1]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'gaps'))).toEqual([1, 1, 0]);
  });

  it('a column whose values all read as numbers draws and orders as numbers; an own `identity` is no percent', async () => {
    expect(numericValues(['9606', null, ' 1e3 '])).toEqual([9606, null, 1000]);
    expect(numericValues(['9606', 'x'])).toBeUndefined();
    const { lay } = await laid(
      [
        { field: 'TAXID', as: 'bar' },
        { field: 'identity', as: 'number' },
      ],
      { layer: { order: { field: 'TAXID', direction: 'ascending' } }, rows: { identity: [0.5, 0.25, 1] } },
    );
    // by value: 7955 (C), 9606 (A), 10090 (B)
    expect(lay.labelText).toEqual(['C/1-4', 'A/1-5', 'B/1-5']);
    expect(shares(lay, 0)).toEqual([7955 / 10090, 9606 / 10090, 1]);
    // the rows' own `identity` shadows the measure, and reads as it is
    expect(cellsOf(lay, 1)).toEqual(['1', '0.5', '0.25']);
  });

  it('draw left in the header and right after the last numbers, alike in the canvas and the SVG, in both themes', async () => {
    const { b, f, lay } = await laid(AllKinds);
    const c = lay.cols;
    expect(c.x[0]).toBeGreaterThan(lay.labelX);
    expect(c.x[1] + c.width[1]).toBeLessThanOrEqual(lay.numberX);
    expect(c.x[2]).toBeGreaterThanOrEqual(lay.box.x + lay.box.width);
    expect(c.x[3]).toBeGreaterThan(c.x[2]);
    expect(Math.max(...lay.lastW.map((w) => lay.numberRightX + w))).toBeLessThanOrEqual(c.x[2]);
    for (const tokens of [undefined, DarkTokens] as const) {
      const s = await settle(b.build());
      if (tokens) s.renderer.tokens = tokens;
      const g = await s.frame({ size: Size, dpr: 1, settled: true });
      const svg = await renderSvg(b.build(), { size: Size, ...(tokens ? { tokens } : {}) });
      expect(svg.issues).toEqual([]);
      expect(recordDrawLists(svg.lists, ['data', 'labels'])).toBe(recordDrawLists(g.lists, ['data', 'labels']));
      for (const t of ['Species', 'Kind', 'identity', 'length', 'Homo sapiens', 'Mus musculus', '50%'])
        expect(svg.svg).toContain(t);
    }
    expect(f.layout.issues).toEqual([]);
  });

  it('titles sit in the layer’s own title row: never over the track above, never off the page as the first track', async () => {
    const { f, lay } = await laid(AllKinds);
    const tracks = f.layout.sections[0].tracks;
    const ruler = tracks[0].box;
    const rowsTrack = tracks[1].box;
    const size = 10;
    // the title row is the rows track's top, the rows below it
    expect(lay.box.y - lay.cols.titleRow).toBe(rowsTrack.y);
    expect(lay.cols.titleY - size).toBeGreaterThanOrEqual(ruler.y + ruler.height);
    // three rows and the title row: the track grew by it
    expect(rowsTrack.height).toBe(3 * 20 + lay.cols.titleRow);
    const { lay: top } = await laid(AllKinds, { first: true });
    expect(top.cols.titleY - size).toBeGreaterThanOrEqual(0);
    expect(titles(top).length).toBe(4);
  });

  it('hit as their row’s header, naming the value; a title hits by its drawn text, the chip is no title', async () => {
    const { f, lay } = await laid(AllKinds);
    const c = lay.cols;
    const y = lay.rowY(1) + lay.rowScale / 2;
    expect(hitView(f.layout, c.x[0] + 2, y)).toMatchObject({
      target: 'row',
      position: 1,
      state: 'Species: Mus musculus',
    });
    expect(hitView(f.layout, c.x[2] + c.width[2] - 2, y)).toMatchObject({
      target: 'row',
      position: 1,
      state: 'identity: 50%',
    });
    expect(hitView(f.layout, lay.labelX + 2, y)?.state).toBeUndefined();
    const t = c.titleCol.indexOf(2);
    expect(hitView(f.layout, c.titleX[t] + 2, c.titleY - 3)).toMatchObject({
      target: 'column-title',
      state: 'identity',
      layer: '/views/v/s/rows/r',
    });
    // past the title's text, within its column: no title
    expect(hitView(f.layout, c.titleX[t] + c.titleW[t] + 2, c.titleY - 3)?.target).not.toBe('column-title');
    // narrowed: the chip in the title row is no title hit
    const { f: g, lay: n } = await laid(AllKinds, {}, { width: 560, height: 300 });
    expect(n.cols.chip).toMatch(/hidden$/);
    expect(hitView(g.layout, n.cols.chipX + 2, n.cols.titleY - 3)?.target).not.toBe('column-title');
  });

  it('an order by a column’s field orders the rows, its title marked; a field no row has is said once', async () => {
    const { lay } = await laid(AllKinds, { layer: { order: { field: 'identity', direction: 'ascending' } } });
    expect(lay.labelText).toEqual(['B/1-5', 'C/1-4', 'A/1-5']);
    expect(titles(lay)[lay.cols.titleCol.indexOf(2)]).toBe('identity ↑');
    const { lay: byOs } = await laid(AllKinds, { layer: { order: { field: 'OS', direction: 'descending' } } });
    expect(byOs.labelText).toEqual(['B/1-5', 'A/1-5', 'C/1-4']);
    const { f } = await laid([{ field: 'nothing' }, { field: 'OS' }]);
    const said = f.layout.issues.filter((i) => i.code === 'seq.render.rows.column');
    expect(said.length).toBe(1);
    expect(said[0].message).toContain('nothing');
  });

  it('join a records source whose system is the rows system, by row key (an id) or by position', async () => {
    const tsv = 'id\tclade\tscore\nC/1-4\tfish\t10\nA/1-5\tprimate\t9\nZ\tnone\t1\n';
    const { lay, f } = await laid(
      [
        { field: 'clade', side: 'right' },
        { field: 'score', side: 'right', as: 'number' },
      ],
      { layer: { order: { field: 'score', direction: 'descending' } }, records: tsv },
    );
    expect(f.layout.issues).toEqual([]);
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'clade'))).toEqual(['primate', undefined, 'fish']);
    expect(lay.labelText).toEqual(['C/1-4', 'A/1-5', 'B/1-5']);
    expect(cellsOf(lay, 1)).toEqual(['10', '9', '']);
    // only the fields a column or the order names are joined
    expect([...(lay.joined?.keys() ?? [])].sort()).toEqual(['clade', 'score']);
    const { lay: byPosition } = await laid([{ field: 'note' }], { records: 'position\tnote\n2\tsecond\n' });
    expect([0, 1, 2].map((r) => byPosition.fieldValue(r, 'note'))).toEqual([undefined, 'second', undefined]);
  });

  it('the joined fields follow the alignment read again (its rows in another order), never the old rows', async () => {
    const tsv = 'id\tclade\nA/1-5\tprimate\nC/1-4\tfish\n';
    const { s, lay } = await laid([{ field: 'clade' }], { records: tsv });
    expect([0, 1, 2].map((r) => lay.fieldValue(r, 'clade'))).toEqual(['primate', undefined, 'fish']);
    // the alignment's source answers again with its rows reversed (the joined source untouched)
    s.store.detach('/data/aln');
    s.store.attach(
      '/data/aln',
      openAlignment({
        label: [...Labels].reverse(),
        letters: [...Letters].reverse(),
        fields: { OS: [...Species].reverse() },
      }),
    );
    const g = await s.frame({ size: Size, dpr: 1, settled: true });
    const again = rowsLays(g as never)[0];
    expect(again).toBe(lay);
    expect(again.label(0)).toBe('C/1-4');
    expect([0, 1, 2].map((r) => again.fieldValue(r, 'clade'))).toEqual(['fish', undefined, 'primate']);
  });

  it('a joined field an own field or a measure shadows, or two sources share, is said once; the first read', async () => {
    const { lay, f } = await laid([{ field: 'identity', as: 'number' }, { field: 'OS' }, { field: 'clade' }], {
      records: 'id\tidentity\tOS\tclade\nA/1-5\t0.3\tX\tprimate\n',
      records2: 'id\tclade\nA/1-5\tape\n',
    });
    const said = f.layout.issues.filter((i) => i.code === 'seq.render.rows.field');
    expect(said.length).toBe(1);
    expect(said[0].message).toContain('identity');
    expect(said[0].message).toContain('OS');
    expect(said[0].message).toContain('clade');
    // the measure and the own field win; the first source's clade
    expect(cellsOf(lay, 0)[0]).toBe('100%');
    expect(cellsOf(lay, 1)[0]).toBe('Homo sapiens');
    expect(lay.fieldValue(0, 'clade')).toBe('primate');
  });

  it('a source that could join and fails holds nothing: the field no row has is said once the view settles', async () => {
    // keyed records need an id: none here, the source refused
    const { f } = await laid([{ field: 'clade' }], { records: 'name\tclade\nA/1-5\tprimate\n' });
    expect(f.layout.issues.some((i) => i.code === 'seq.render.rows.column')).toBe(true);
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
    const widths = async (width: number) => {
      const { lay } = await laid(AllKinds, {}, { width, height: 300 });
      return { w: [...lay.cols.width], chip: lay.cols.chip, x: lay.box.x, chipEnd: lay.cols.chipX + lay.cols.chipW };
    };
    const wide = await widths(1000);
    expect(wide.w.every((w) => w > 0)).toBe(true);
    expect(wide.chip).toBe('');
    const mid = await widths(640);
    expect(mid.w[3]).toBe(0);
    expect(mid.w[0]).toBeGreaterThan(0);
    expect(mid.chip).toMatch(/columns? hidden$/);
    // the chip ends on the page, before the right inset
    expect(mid.chipEnd).toBeLessThanOrEqual(640);
    const narrow = await widths(400);
    expect(narrow.w).toEqual([0, 0, 0, 0]);
    expect(narrow.chip).toBe('4 columns hidden');
    const { lay: none } = await laid([], {}, { width: 400, height: 300 });
    expect(narrow.x).toBe(none.box.x);
  });

  it('a text figure draws the same columns on every line, their titles over the first (and hit there)', async () => {
    const section = { layout: 'text', group: 3, max_width: 48 };
    const { f } = await laid(AllKinds, { section });
    const lays = rowsLays(f as never);
    expect(lays.length).toBeGreaterThan(1);
    for (const lay of lays) {
      expect(lay.figure).toBe(true);
      expect(cellsOf(lay, 0)).toEqual(['Homo sapiens', 'Mus musculus', '']);
      expect(cellsOf(lay, 2)).toEqual(['100%', '50%', '50%']);
    }
    expect(titles(lays[0]).length).toBe(4);
    expect(titles(lays[1]).length).toBe(0);
    const first = lays[0].cols;
    const strips = f.layout.sections[0].placement.strips;
    const t = first.titleCol.indexOf(0);
    expect(hitView(f.layout, first.titleX[t] + 2 + strips.dx[0], first.titleY - 3 + strips.dy[0])).toMatchObject({
      target: 'column-title',
      state: 'OS',
    });
    const svg = await renderSvg(spec(AllKinds, { section }).build(), { size: Size });
    expect(recordDrawLists(svg.lists, ['labels'])).toBe(recordDrawLists(f.lists, ['labels']));
  });

  describe('a hovered cell read whole (the tip pass)', () => {
    const Narrow: Col[] = [
      { field: 'OS', title: 'Species', width: 40 },
      { field: 'identity', side: 'right', as: 'number' },
    ];

    /** The tip pass's lines (`rect x y w h …`, `text "…" x y …`) of the row 0 hovered. */
    async function tip(rows: Record<string, unknown>, width = Size.width) {
      const { s, lay } = await laid(Narrow, { rows }, { width, height: 300 });
      // the hover as the pointer makes it: the hit on the first row's species cell
      const first = await s.frame({ size: { width, height: 300 }, dpr: 1, settled: true });
      const hover = hitView(first.layout, lay.cols.x[0] + 4, lay.rowY(0) + lay.rowScale / 2) as never;
      expect(hover).toBeTruthy();
      const f = s.renderer.frame({ size: { width, height: 300 }, dpr: 1, settled: true, emphasis: { hover } });
      const lines = recordDrawLists(f.lists, ['tip']).split('\n');
      return { lay, f, lines, hoverPass: recordDrawLists(f.lists, ['hover']) };
    }
    const nums = (l: string) =>
      l
        .split(' ')
        .filter((x) => x !== '' && !Number.isNaN(Number(x)))
        .map(Number);

    it('is drawn after every label, on a surface over the columns it overlaps, in the hover pass nowhere', async () => {
      expect(Passes.indexOf('tip')).toBeGreaterThan(Passes.indexOf('labels'));
      const { lay, lines, hoverPass } = await tip({ OS: ['Homo sapiens (human)', 'Mus musculus', null] });
      const texts = lines.filter((l) => l.startsWith('text "'));
      expect(texts).toHaveLength(1);
      expect(texts[0]).toContain('"Homo sapiens (human)"');
      expect(hoverPass).not.toContain('Homo sapiens');
      // the surface: the page's background and then the track's well, both over the cell's x to the text's end
      const rects = lines.filter((l) => l.startsWith('rect '));
      expect(rects).toHaveLength(2);
      const [x, , w] = nums(rects[1]);
      const next = lay.cols.x[0] + lay.cols.width[0];
      expect(x).toBeLessThanOrEqual(lay.cols.x[0]);
      expect(x + w).toBeGreaterThan(next);
      expect(nums(rects[0]).slice(0, 4)).toEqual(nums(rects[1]).slice(0, 4));
      // the text is the last thing drawn
      expect(lines.filter((l) => l !== '').at(-2)).toBe(texts[0]);
    });

    it('is clipped to the band’s scrolled area: a half-scrolled row never paints over the pinned band', async () => {
      const rows = { OS: ['Homo sapiens (human)', 'Mus musculus (house mouse)', 'Danio rerio (zebrafish)'] };
      const { s, lay } = await laid(Narrow, { rows, layer: { pinned: ['A/1-5'], height: 2 } });
      const viewports = { '/views/v/s': { start: 0, end: 6, bands: { '/views/v/s/rows': { bandTop: 1.5 } } } };
      const size = Size;
      const first = await s.frame({ size, dpr: 1, settled: true, viewports } as never);
      // data row 1 (display row 2, the pinned row standing in its place before it), half of it (10 px) scrolled under the
      // pinned band
      const l2 = rowsLays(first as never)[0];
      expect(l2.rowY(2)).toBe(l2.box.y + l2.pinnedHeight - 10);
      const hover = { layer: l2.key, section: 's', system: l2.rowsSystem, position: 1 } as never;
      const f = s.renderer.frame({ size, dpr: 1, settled: true, viewports, emphasis: { hover } } as never);
      const lines = recordDrawLists(f.lists, ['tip']).split('\n');
      const clip = lines.find((l) => l.startsWith('clip '))!;
      expect(clip, lines.join('\n')).toBeTruthy();
      const [, cy, , ch] = nums(clip);
      expect(cy).toBe(lay.box.y + l2.pinnedHeight);
      expect(cy + ch).toBe(lay.box.y + lay.box.height);
      // the surface starts above the clip (the row is half out), so the clip is what keeps it off the pinned band
      expect(nums(lines.find((l) => l.startsWith('rect '))!)[1]).toBeLessThan(cy);
    });

    it('a label longer than the view ends in an ellipsis at the view’s edge', async () => {
      const width = 640;
      const { lines } = await tip({ OS: [`Homo sapiens ${'very long '.repeat(40)}`, 'Mus musculus', null] }, width);
      const text = lines.find((l) => l.startsWith('text "'))!;
      expect(/^text "[^"]*…"/.test(text)).toBe(true);
      const rect = lines.filter((l) => l.startsWith('rect '))[1];
      const [x, , w] = nums(rect);
      expect(x + w).toBeLessThanOrEqual(width);
      expect(x + w).toBeGreaterThan(width - 12);
    });
  });
});
