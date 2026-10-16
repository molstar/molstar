/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// LY3a (seqspec:viewer#extension-points; seqspec:format#themes): each built-in theme against hand-computed colors,
// palette parity through the recorded MVS outputs, Clustal X against Jalview's thresholds (now a `rules` table),
// conservation shading's alpha as a byte, nested bases, theme and glyph fallbacks reported; then `rules` itself (every
// measure, the first rule that holds, alpha, labels in the legend, granularity), a size channel, a continuous palette's
// ramp, the element table and the theme node on case 3.

import { readFileSync } from 'node:fs';
import { AminoGroupColors, Color, LetterDictionaries, type Palette, type ThemeRef } from '@molstar/seq-core/format';
import { type KeyTarget, type SeqSet, SeqSets, type SetInput } from '@molstar/seq-core/model';
import { describe, expect, it } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import { trackPositions } from '../src/actions.js';
import { GlyphShape, NO_COLOR } from '../src/draw/backend.js';
import { createSeqViewer, type Theme, ThemeLocation } from '../src/index.js';
import { DefaultSeqViewerSpec } from '../src/default-spec.js';
import {
  CellInk,
  cellInk,
  createTheme,
  GlyphTable,
  Granularity,
  type LegendSink,
  NO_GLYPH,
  NO_SIZE,
  OPAQUE,
  SIZE_ONE,
  themeContext,
} from '../src/themes.js';
import { layer, opened } from './layout-fixtures.js';
import { at, channel, context } from './theme-fixtures.js';

const t: KeyTarget = { kind: -1, id: -1, row: -1 };
const idOf = (set: SeqSet, key: string) => (set.byName(key, t) ? t.id : -1);

const INPUT: SetInput = {
  molecules: [
    { key: 'p', residues: 'ACDEFGHIKL' },
    { key: 'd', alphabet: 'dna', residues: 'ACGTACGT' },
  ],
  units: [
    { key: 'A', molecule: 'p' },
    { key: 'B', molecule: 'p' },
    { key: 'C', molecule: 'd' },
  ],
  alignments: [{ key: 'aln', gapped: ['ACD-E', 'ACDKE', 'WCD-E', 'ACKKE'], flavor: 'aligned', alphabet: 'protein' }],
  annotations: [
    {
      key: 'f',
      geometry: 'ranges',
      target: 'p',
      columns: {
        start: [0, 2, 4, 6],
        stop: [2, 4, 6, 8],
        category: ['x', 'y', 'x', 'z'],
        color: ['#ff0000', 'nope', 'blue', null],
        glyph: ['diamond', '★', 'square', null],
        score: [1, 5, 9, 2],
      },
    },
    { key: 'v', geometry: 'values', target: 'p', columns: { index: [0, 1, 1], value: [0, 10, 5] } },
    {
      key: 'cons',
      geometry: 'values',
      target: 'aln',
      columns: { start: [0, 1, 2, 3, 4], value: [1, 0.2, 0.75, 0.5, 2] },
    },
  ],
};

const C = (s: string) => Color.decodePrimary(s);
const code = (ch: string) => ch.charCodeAt(0);
/** A legend's writes, as data. */
function legendOf(theme: Theme | null) {
  const out: unknown[] = [];
  const sink: LegendSink = {
    item: (label, color, glyph) => out.push(['item', label, color, glyph]),
    ramp: (stops, lo, hi) => out.push(['ramp', Array.from(stops), lo, hi]),
    size: (label, value, unit) => out.push(['size', label, value, unit]),
  };
  theme?.legend?.(sink);
  return out;
}

describe('color themes, against hand-computed colors', () => {
  const set = SeqSets.build(INPUT).set;
  const f = idOf(set, 'f');

  it('uniform, field_color, field_palette (category by first occurrence; a row read at its data row)', () => {
    const x = context(set, { data: ['f'] });
    expect(x.color('uniform', { color: '#123456' }).value(at({}))).toBe(0x123456);
    const fc = x.color('field_color');
    expect([0, 1, 2, 3].map((row) => fc.value(at({ annotation: f, row })))).toEqual([
      0xff0000,
      NO_COLOR,
      0x0000ff,
      NO_COLOR,
    ]);
    const fp = x.color('field_palette');
    const c10 = ['#1f77b4', '#ff7f0e', '#2ca02c'].map(C);
    expect([0, 1, 2, 3].map((row) => fp.value(at({ annotation: f, row })))).toEqual([c10[0], c10[1], c10[0], c10[2]]);
    expect(fp.value(at({ annotation: f + 99, row: 0 }))).toBe(NO_COLOR);
    expect([fp.granularity, fp.reduction]).toEqual([Granularity.item, 'classes']);
    // A values layer: `value` through a continuous palette over the layer's values (0..10); a `mean` theme.
    const v = context(set, { data: ['v'] }).color('field_palette', {
      palette: { kind: 'continuous', colors: ['#000000', '#ffffff'] } satisfies Palette,
    });
    expect([0, 5, 10].map((value) => v.value(at({ value })))).toEqual([0x000000, 0x7f7f7f, 0xffffff]);
    expect([v.granularity, v.reduction]).toEqual([Granularity.value, 'mean']);
  });

  it("letter: the alphabet's dictionary (either case), a named one; atoms by the model's element table", () => {
    const x = context(set, { blocks: [{ unit: 'A' }, { unit: 'C' }] });
    const l = x.color('letter');
    expect(l.value(at({ block: 0, letter: code('A') }))).toBe(LetterDictionaries.Zappo.A);
    expect(l.value(at({ block: 0, letter: code('a') }))).toBe(LetterDictionaries.Zappo.A);
    expect(l.value(at({ block: 1, letter: code('G') }))).toBe(LetterDictionaries.Nucleotide.G);
    expect(l.value(at({ block: 0, letter: code('-') }))).toBe(NO_COLOR);
    expect(x.color('letter', { colors: 'Taylor' }).value(at({ block: 1, letter: code('C') }))).toBe(
      LetterDictionaries.Taylor.C,
    );
    // Atoms by atomic number past barium, with the default or given dictionaries alike.
    const lig = SeqSets.build({
      molecules: [{ key: 'l', kind: 'ligand', elements: ['C', 'N', 'Fe', 'W', 'Pt', 'Au'], ownUnit: true }],
    }).set;
    const atoms = (colors?: unknown) => context(lig, { blocks: [{ unit: 'l/unit' }] }).color('letter', { colors });
    const Z = [6, 7, 26, 74, 78, 79];
    expect(Z.map((z) => atoms().value(at({ block: 0, letter: z })))).toEqual(
      ['#909090', '#3050f8', '#e06633', '#2194d6', '#d0d0e0', '#ffd123'].map(C),
    );
    expect(Z.map((z) => atoms('ElementSymbol').value(at({ block: 0, letter: z })))).toEqual(
      Z.map((z) => atoms().value(at({ block: 0, letter: z }))),
    );
    const given = atoms({ c: '#010101', Fe: '#020202', PT: '#030303' });
    expect([6, 26, 78, 7].map((z) => given.value(at({ block: 0, letter: z })))).toEqual([
      0x010101, 0x020202, 0x030303, -1,
    ]);
  });

  it('by_unit and by_copy_group: one color per unit, per group by its first unit', () => {
    const x = context(set);
    const u = x.color('by_unit');
    const g = x.color('by_copy_group');
    const [a, b, c] = ['A', 'B', 'C'].map((k) => idOf(set, k));
    const c10 = ['#1f77b4', '#ff7f0e', '#2ca02c'].map(C);
    expect([a, b, c].map((unit) => u.value(at({ unit })))).toEqual([c10[0], c10[1], c10[2]]);
    expect([a, b, c].map((unit) => g.value(at({ unit })))).toEqual([c10[0], c10[0], c10[1]]);
  });

  it("identity_to_row over a base; conservation_shading's alpha is the column's value as a byte, none below", () => {
    const x = context(set, { blocks: [{ alignment: 'aln' }] });
    const id = x.color('identity_to_row', { row: 0, base: { theme: 'uniform', params: { color: '#00ff00' } } });
    // Column 0: row 0 has A; W differs. Column 3: row 0 has a gap.
    expect(id.value(at({ column: 0, letter: code('A') }))).toBe(0x00ff00);
    expect(id.value(at({ column: 0, letter: code('a') }))).toBe(0x00ff00);
    expect(id.value(at({ column: 0, letter: code('W') }))).toBe(NO_COLOR);
    expect(id.value(at({ column: 3, letter: code('K') }))).toBe(NO_COLOR);
    // It varies by column, so the letters' cell classes are never keyed by the letter alone.
    expect(id.granularity).toBe(Granularity.column);
    expect(x.color('identity_to_row', { row: 0 }).granularity).toBe(Granularity.column);
    const cs = x.color('conservation_shading', { data: 'cons', threshold: 0.5 });
    const columns = [0, 1, 2, 3, 4];
    expect(columns.map((column) => cs.value(at({ column, letter: code('C') })))).toEqual(
      columns.map((k) => (k === 1 ? NO_COLOR : LetterDictionaries.Zappo.C)),
    );
    expect(columns.map((column) => cs.alpha?.(at({ column, letter: code('C') })))).toEqual([255, 0, 191, 128, 255]);
  });

  it('nested bases resolve once, at most 4 deep (deeper: svs.color.invalid)', () => {
    const x = context(set, { blocks: [{ alignment: 'aln' }] });
    const nest = (n: number): ThemeRef =>
      n === 0
        ? { theme: 'uniform', params: { color: '#0000ff' } }
        : { theme: 'identity_to_row', params: { base: nest(n - 1) } };
    expect(x.color('identity_to_row', { base: nest(1) }).value(at({ column: 1, letter: code('C') }))).toBe(0x0000ff);
    expect(x.issues).toEqual([]);
    x.color('identity_to_row', { base: nest(4) });
    expect(x.issues.map((i) => i.code)).toEqual(['svs.color.invalid']);
  });

  it("last row wins where values overlap: the position's value, then its color", async () => {
    const { Annotations } = await import('@molstar/seq-core/model');
    const v = set.annotations[idOf(set, 'v')] as NonNullable<SeqSet['annotations'][number]>;
    const value = Annotations.valueAt(v, 0, 1);
    expect(value).toBe(5);
    const fp = context(set, { data: ['v'] }).color('field_palette', {
      palette: { kind: 'categorical', colors: ['red', 'green', 'blue'] } satisfies Palette,
    });
    expect(fp.value(at({ value }))).toBe(fp.value(at({ annotation: v.id, row: 2 })));
  });

  it("a theme without an implementation falls back (its own, then its entry's) and says so; none: null", () => {
    const x = context(set, { blocks: [{ unit: 'C' }] });
    const red = x.color('lab:heat', {}, { theme: 'uniform', params: { color: 'red' } });
    expect(red.value(at({}))).toBe(0xff0000);
    expect(x.color('lab:none')).toBeNull();
    expect(x.issues.map((i) => i.code)).toEqual(['svs.render.theme', 'svs.render.theme']);
  });
});

describe('glyph themes', () => {
  const set = SeqSets.build(INPUT).set;
  const f = idOf(set, 'f');
  it("a track's positions list a data ref and filter once, the filter's keys in any order (VS-20)", () => {
    const layer = (filter: unknown) => ({ kind: 'layer', params: { type: 'ranges', data: 'f', filter } });
    const track = {
      kind: 'track',
      children: [layer({ field: 'kind', values: ['a'] }), layer({ values: ['a'], field: 'kind' }), layer(undefined)],
    };
    expect(trackPositions(track as never)).toEqual([{ data: 'f', filter: { field: 'kind', values: ['a'] } }, 'f']);
  });

  it('a glyph object interned again with its keys in another order is the same glyph (VS-20)', () => {
    const g = new GlyphTable();
    const a = g.intern({ char: 'x', fallback: 'square' });
    expect(g.intern({ fallback: 'square', char: 'x' })).toBe(a);
    expect(g.intern({ char: 'y', fallback: 'square' })).not.toBe(a);
  });

  it('uniform, field_map with a default, field_glyph (shapes and one character), a fallback reported', () => {
    const x = context(set, { data: ['f'] });
    expect(x.glyph('uniform', { glyph: 'star' }).value(at({}))).toBe(GlyphShape.star);
    const m = x.glyph('field_map', { map: { x: 'diamond', y: { char: '✶' } }, default: 'cross' });
    const ids = [0, 1, 2, 3].map((row) => m.value(at({ annotation: f, row })));
    expect([ids[0], ids[2], ids[3]]).toEqual([1, 1, GlyphShape.cross]);
    expect([x.tc.glyphs.glyph(ids[1])?.char, x.tc.glyphs.shape(ids[1])]).toEqual(['✶', 0]);
    const g = x.glyph('field_glyph');
    const gs = [0, 1, 2, 3].map((row) => g.value(at({ annotation: f, row })));
    expect([gs[0], x.tc.glyphs.glyph(gs[1])?.char, gs[2], gs[3]]).toEqual([1, '★', 3, NO_GLYPH]);
    expect(x.glyph('lab:marks')).toBeNull();
    expect(x.issues.map((i) => i.code)).toEqual(['svs.render.theme']);
  });

  it('field_glyph takes any one grapheme as a character, lower case and accented too (VS-16)', () => {
    const x = context(set);
    const g = x.glyph('field_glyph');
    const chars = ['x', 'X', 'é', '★', '👍🏽'].map(
      (c) => x.tc.glyphs.glyph(g.raw?.({ glyph: [c] }, 0) ?? NO_GLYPH)?.char,
    );
    expect(chars).toEqual(['x', 'X', 'é', '★', '👍🏽']);
    expect(g.raw?.({ glyph: ['star'] }, 0)).toBe(GlyphShape.star);
    expect(g.raw?.({ glyph: ['xy'] }, 0)).toBe(NO_GLYPH);
  });
});

describe('palette parity through the recorded MVS outputs', () => {
  const recorded = JSON.parse(
    readFileSync(new URL('../../core/_test/format/palettes/mvs-recorded.json', import.meta.url), 'utf8'),
  ) as {
    cases: {
      name: string;
      palette: Palette;
      values: (string | number)[];
      queries: (string | number)[];
      expected: (string | null)[];
    }[];
  };
  for (const c of recorded.cases) {
    it(c.name, () => {
      const numeric = c.palette.kind !== 'categorical';
      const values = c.values.map(String);
      const set = SeqSets.build({
        molecules: [{ key: 'p', residues: 'A'.repeat(values.length) }],
        annotations: [
          {
            key: 'a',
            geometry: 'ranges',
            target: 'p',
            columns: { start: values.map((_, k) => k), stop: values.map((_, k) => k + 1), v: c.values },
          },
        ],
      }).set;
      // An MVS palette in SVS: MVS's default `repeat_color_list` (false) written out, SVS's being true
      // (seqspec:format#palettes).
      const palette = numeric ? c.palette : { repeat_color_list: false, ...c.palette };
      const theme = context(set, { data: ['a'] }).color('field_palette', { field: 'v', palette });
      const hex = (x: number) => (x < 0 ? null : Color.toHex(x));
      const a = idOf(set, 'a');
      c.queries.forEach((q, k) => {
        if (numeric && typeof q === 'number')
          expect([q, hex(theme.value(at({ value: q })))]).toEqual([q, c.expected[k]]);
        const row = values.indexOf(String(q));
        if (!numeric && row >= 0) expect([q, hex(theme.value(at({ annotation: a, row })))]).toEqual([q, c.expected[k]]);
      });
    });
  }
});

describe("clustal_x, a rule table, against Jalview's thresholds (a fixture column set)", () => {
  // Ten rows, five columns: A 70% (hydrophobic > 60%); A 60% (not above); C 90% (> 85%: cysteine; the A beside it
  // hydrophobic, C being in the group); K+R 70% (positive, and E negative by KR > 60%); G with L 90% (glycine).
  const cols = ['AAAAAAAKKK', 'AAAAAAKKKK', 'CCCCCCCCCA', 'KKKKRRREEE', 'GLLLLLLLLL'];
  const rows = Array.from({ length: 10 }, (_, r) => cols.map((c) => c[r]).join(''));
  const set = SeqSets.build({
    molecules: [{ key: 'p', residues: 'MKLAW' }],
    units: [{ key: 'A', molecule: 'p' }],
    alignments: [{ key: 'j', gapped: rows, flavor: 'aligned', alphabet: 'protein' }],
  }).set;
  const G = AminoGroupColors;
  it('each residue by the first rule its column meets', () => {
    const cx = context(set, { blocks: [{ alignment: 'j' }] }).color('clustal_x');
    const of = (column: number, ch: string) => cx.value(at({ column, letter: code(ch) }));
    expect([of(0, 'A'), of(0, 'K'), of(1, 'A'), of(2, 'C'), of(2, 'A')]).toEqual([
      G.hydrophobic,
      NO_COLOR,
      NO_COLOR,
      G.cysteine,
      G.hydrophobic,
    ]);
    expect([of(3, 'K'), of(3, 'R'), of(3, 'E'), of(4, 'G'), of(4, 'L'), of(4, '-')]).toEqual([
      G.positive,
      G.positive,
      G.negative,
      G.glycine,
      G.hydrophobic,
      NO_COLOR,
    ]);
    expect(cx.granularity).toBe(Granularity.column);
  });
  it('on a protein sequence a position is a column of its own letter (as the hand-written theme)', () => {
    const cx = context(set, { blocks: [{ unit: 'A' }] }).color('clustal_x');
    expect([...'MKLAWGPC'].map((ch) => cx.value(at({ block: 0, column: -1, letter: code(ch) })))).toEqual([
      G.hydrophobic,
      G.positive,
      G.hydrophobic,
      G.hydrophobic,
      G.hydrophobic,
      G.glycine,
      G.proline,
      G.cysteine,
    ]);
  });
});

describe('rules', () => {
  const set = SeqSets.build(INPUT).set;
  const f = idOf(set, 'f');
  const x = context(set, { data: ['f'], blocks: [{ alignment: 'aln' }] });
  const rules = (rs: unknown[], more: Record<string, unknown> = {}) => x.color('rules', { rules: rs, ...more });

  it('the first rule that holds; `value`, a number and a string `field` (with `in`), `column`; default', () => {
    const th = rules(
      [
        { when: [{ measure: { field: 'score' }, min: 5 }], color: '#ff0000', label: 'high' },
        { when: [{ measure: { field: 'category' }, in: ['z'] }], color: '#00ff00', label: 'z' },
        { when: [{ measure: 'value', max: 0.5 }], color: '#0000ff' },
        { when: [{ measure: { column: 'cons' }, min: 0.6, max: 1 }], color: '#ffff00' },
      ],
      { default: { theme: 'uniform', params: { color: '#808080' } } },
    );
    expect([0, 1, 2, 3].map((row) => th.value(at({ annotation: f, row })))).toEqual([
      0x808080, 0xff0000, 0xff0000, 0x00ff00,
    ]);
    expect([th.value(at({ value: 0.25 })), th.value(at({ value: 0.75 }))]).toEqual([0x0000ff, 0x808080]);
    expect([0, 1, 2].map((column) => th.value(at({ column })))).toEqual([0xffff00, 0x808080, 0xffff00]);
    expect(th.granularity).toBe(Granularity.position);
    expect(legendOf(th)).toEqual([
      ['item', 'high', 0xff0000, NO_GLYPH],
      ['item', 'z', 0x00ff00, NO_GLYPH],
    ]);
  });

  it("a `column` measure over data on another alignment is NaN: the default, never that alignment's values (VS-04)", () => {
    const two = SeqSets.build({
      ...INPUT,
      alignments: [
        ...(INPUT.alignments ?? []),
        { key: 'aln2', gapped: ['ACDEF', 'ACDKF'], flavor: 'aligned', alphabet: 'protein' },
      ],
      annotations: [
        {
          key: 'cons2',
          geometry: 'values',
          target: 'aln2',
          columns: { start: [0, 1, 2, 3, 4], value: [1, 1, 1, 1, 1] },
        },
      ],
    }).set;
    const y = context(two, { blocks: [{ alignment: 'aln' }] });
    const red = [{ when: [{ measure: { column: 'cons2' }, min: 0.5 }], color: '#ff0000' }];
    const th = y.color('rules', { rules: red, default: { theme: 'uniform', params: { color: '#808080' } } });
    expect([0, 1, 2].map((column) => th.value(at({ column })))).toEqual([0x808080, 0x808080, 0x808080]);
    const on2 = context(two, { blocks: [{ alignment: 'aln2' }] }).color('rules', { rules: red });
    expect(on2.value(at({ column: 1 }))).toBe(0xff0000);
  });

  it('`share` of the letter and of a set, `rank`, `equals_row`, `in` on the letter; NaN fails a clause', () => {
    // Column 2: D D D K (D 0.75, K 0.25); column 3: - K - K (K 1.0 of letters).
    const share = (m: unknown, min: number) => rules([{ when: [{ measure: m, min }], color: '#010101' }]);
    const s = share({ share: 'letter' }, 0.75);
    expect([s.value(at({ column: 2, letter: code('D') })), s.value(at({ column: 2, letter: code('K') }))]).toEqual([
      0x010101, -1,
    ]);
    expect(s.value(at({ column: 3, letter: code('k') }))).toBe(0x010101);
    expect(s.value(at({ column: 3, letter: code('-') }))).toBe(-1);
    const set2 = share({ share: 'DK' }, 1);
    expect([
      set2.value(at({ column: 2, letter: code('D') })),
      set2.value(at({ column: 0, letter: code('A') })),
    ]).toEqual([0x010101, -1]);
    const rank = rules([{ when: [{ measure: { rank: 'letter' }, max: 0 }], color: '#020202' }]);
    expect([
      rank.value(at({ column: 2, letter: code('D') })),
      rank.value(at({ column: 2, letter: code('K') })),
    ]).toEqual([0x020202, -1]);
    const eq = rules([{ when: [{ measure: { equals_row: 2 }, min: 1 }, { in: ['W', 'C'] }], color: '#030303' }]);
    // Row 2 is WCD-E: W in column 0, C in column 1 (`in` lets only W and C through).
    expect([0, 1, 2].map((column) => eq.value(at({ column, letter: code('WCD'[column]) })))).toEqual([
      0x030303, 0x030303, -1,
    ]);
    expect(eq.granularity).toBe(Granularity.column);
  });

  it("alpha as a byte by its measure; a nested theme's value; glyph and size channels; `rules` is `classes`", () => {
    const fade = rules([{ when: [], theme: { theme: 'letter' }, alpha: { column: 'cons' } }]);
    expect([0, 1, 2].map((column) => fade.alpha?.(at({ column, letter: code('C') })))).toEqual([255, 51, 191]);
    expect(fade.value(at({ column: 1, letter: code('C') }))).toBe(LetterDictionaries.Zappo.C);
    expect([fade.reduction, fade.granularity]).toEqual(['classes', Granularity.column]);
    const g = x.glyph('rules', {
      rules: [{ when: [{ measure: { field: 'score' }, min: 5 }], glyph: 'star', label: 'strong' }],
      default: { theme: 'uniform', params: { glyph: 'tick' } },
    });
    expect([0, 1].map((row) => g.value(at({ annotation: f, row })))).toEqual([GlyphShape.tick, GlyphShape.star]);
    expect(legendOf(g)).toEqual([['item', 'strong', NO_COLOR, GlyphShape.star]]);
  });

  it("a size channel: `uniform` and `rules` serve it, in 1/256 of the unit (an int), none the viewer's", () => {
    const uni = x.size('uniform', { size: 2.5 });
    expect([uni.value(at({})), uni.kind]).toEqual([640, 'size']);
    const th = x.size('rules', {
      rules: [{ when: [{ measure: { field: 'score' }, min: 5 }], size: 2.5, label: '≥ 5' }],
    });
    expect([0, 1].map((row) => th.value(at({ annotation: f, row })))).toEqual([NO_SIZE, 640]);
    expect(legendOf(th)).toEqual([['size', '≥ 5', 2.5, 0]]);
    expect(Number.isInteger(th.value(at({ annotation: f, row: 1 })))).toBe(true);
  });

  it('`field_scale`: a number through linear, sqrt and log into the range, clamped; a string through a map', () => {
    const scale = (params: Record<string, unknown>) => x.size('field_scale', { range: [2, 6], ...params });
    const at4 = (th: Theme) => [0, 1, 2, 3].map((row) => th.value(at({ annotation: f, row })) / SIZE_ONE);
    // score is [1, 5, 9, 2]; the domain from the data (1..9)
    expect(at4(scale({ field: 'score' }))).toEqual([2, 4, 6, 2.5]);
    expect(at4(scale({ field: 'score', domain: [3, 7] }))).toEqual([2, 4, 6, 2]);
    expect(at4(scale({ field: 'score', scale: 'sqrt' }))[1]).toBeCloseTo(2 + 4 * Math.sqrt(0.5), 1);
    expect(at4(scale({ field: 'score', scale: 'log' }))[1]).toBeCloseTo(2 + (4 * Math.log1p(4)) / Math.log1p(8), 1);
    const mapped = scale({ field: 'category', map: { x: 3, y: 5 }, default: 1 });
    expect(at4(mapped)).toEqual([3, 5, 3, 1]);
    expect(at4(scale({ field: 'category', map: { x: 3 } }))[1]).toBe(NO_SIZE / SIZE_ONE);
    // Far below the domain (VS-04): clamped to its low end on every scale, never no size.
    expect(at4(scale({ field: 'score', domain: [10, 100], scale: 'log' }))).toEqual([2, 2, 2, 2]);
    expect(at4(scale({ field: 'score', domain: [10, 100], scale: 'sqrt' }))).toEqual([2, 2, 2, 2]);
    // The value of a position with no row (a values layer's), over the layer's domain.
    const v = context(set, { data: ['v'], domain: [0, 10] }).size('field_scale', { range: [0, 10] });
    expect([0, 5, 20].map((value) => v.value(at({ value })) / SIZE_ONE)).toEqual([0, 5, 10]);
  });

  it("EN1: a string field provider's name in an `in` clause; the rule found for a value serves its alpha", () => {
    const registry = createSeqViewer(DefaultSeqViewerSpec()).registry as ReturnType<
      typeof createSeqViewer
    >['registry'] & {
      register(e: unknown): () => void;
    };
    registry.register({
      svs: {
        fields: [
          {
            name: 'lab:parity',
            type: 'string',
            applies: () => true,
            read: (a: { count: number }, out: Int32Array, strings: { intern(s: string): number }) => {
              for (let r = 0; r < a.count; r++) out[r] = strings.intern(r % 2 ? 'odd' : 'even');
            },
          },
        ],
      },
    } as never);
    const tc = themeContext(x.tc, false, registry, () => {});
    const th = createTheme(
      tc,
      {
        theme: 'rules',
        params: {
          rules: [
            { when: [{ measure: { field: 'lab:parity' }, in: ['odd'] }], color: '#ff0000', alpha: { column: 'cons' } },
            { when: [{ measure: { column: 'cons' }, min: 0.6 }], color: '#0000ff', alpha: { column: 'cons' } },
          ],
        },
      },
      channel('color'),
      registry,
      () => {},
    ) as Theme;
    expect([0, 1, 2, 3].map((row) => th.value(at({ annotation: f, row })))).toEqual([
      NO_COLOR,
      0xff0000,
      NO_COLOR,
      0xff0000,
    ]);
    // Value then alpha of one location share the search; the next location, the next search (never a stale rule).
    const trace = [
      { annotation: f, row: 1, column: 1 },
      { annotation: f, row: 0, column: 1 },
      { annotation: f, row: 0, column: 0 },
      { annotation: f, row: 1, column: 0 },
    ];
    expect(
      trace.map((l) => [th.value(at(l)), th.alpha?.(at(l))]).map(([c, a]) => [c === NO_COLOR ? 'none' : c, a]),
    ).toEqual([
      [0xff0000, 51],
      ['none', 0],
      [0x0000ff, 255],
      [0xff0000, 255],
    ]);
    // A different location asked of `alpha` alone, after another's `value`.
    th.value(at({ annotation: f, row: 1, column: 0 }));
    expect(th.alpha?.(at({ annotation: f, row: 0, column: 1 }))).toBe(0);
  });

  it('`equals_row` on a lazy alignment reads its row where it is held; else not applicable, reported (VS-16)', () => {
    const build = (held: boolean) =>
      SeqSets.build({
        molecules: [],
        alignments: [
          {
            key: 'z',
            gapped: [],
            flavor: 'aligned',
            alphabet: 'protein',
            lazy: {
              rowCount: 4,
              columnCount: 5,
              rows: new Map(held ? [[0, { label: 'r0', gapped: 'ACDEF' }]] : []),
              labelOf: () => null,
            },
          },
        ],
      } as never).set;
    const rule = [{ when: [{ measure: { equals_row: 0 }, min: 1 }], color: '#ff0000' }];
    const off = context(build(false), { blocks: [{ alignment: 'z' }] });
    off.color('rules', { rules: rule });
    expect(off.issues.map((i) => i.code)).toEqual(['svs.render.theme']);
    const on = context(build(true), { blocks: [{ alignment: 'z' }] });
    const th = on.color('rules', { rules: rule });
    expect(on.issues).toEqual([]);
    expect([
      th.value(at({ block: 0, column: 1, letter: code('C') })),
      th.value(at({ block: 0, column: 1, letter: code('K') })),
    ]).toEqual([0xff0000, NO_COLOR]);
  });

  it('`share` on a lazy alignment is not applicable: `letter`, reported; with a profile it is', () => {
    const lazy = SeqSets.build({
      molecules: [],
      alignments: [
        {
          key: 'z',
          gapped: [],
          flavor: 'aligned',
          alphabet: 'protein',
          lazy: { rowCount: 4, columnCount: 5, rows: new Map(), labelOf: () => null },
        },
      ],
    }).set;
    const y = context(lazy, { blocks: [{ alignment: 'z' }] });
    const th = y.color('clustal_x');
    expect(th.value(at({ block: 0, letter: code('A') }))).toBe(LetterDictionaries.Zappo.A);
    expect(y.issues.map((i) => [i.code, i.message])).toEqual([
      ['svs.render.theme', '"clustal_x" is drawn as "letter"'],
    ]);
  });
});

describe("legends: a continuous palette's ramp; categorical keys", () => {
  it("a matrix 0..20 under Viridis: a ramp over the layer's domain, 15 at t = 0.75", () => {
    const set = SeqSets.build({
      molecules: [{ key: 'p', residues: 'ACDEF' }],
      annotations: [
        {
          key: 'm',
          geometry: 'matrix',
          target: 'p',
          states: ['a', 'b'],
          columns: {
            index: [0, 1, 2],
            value: [
              [0, 5],
              [10, 15],
              [20, 1],
            ],
          },
        } as never,
      ],
    }).set;
    const th = context(set, { data: ['m'], domain: [0, 20] }).color('field_palette', {
      palette: { kind: 'continuous', colors: 'Viridis' },
    });
    const [ramp] = legendOf(th) as [string, number[], number, number][];
    expect([ramp[0], ramp[2], ramp[3], ramp[1].length]).toEqual(['ramp', 0, 20, 32]);
    expect(th.value(at({ value: 15 }))).toBe(th.value(at({ value: 0.75 * 20 })));
    expect(ramp[1][2 * 15 + 1]).toBe(th.value(at({ value: 20 })));
    expect(th.value(at({ value: 15 }))).not.toBe(th.value(at({ value: 20 })));
  });
  it("a categorical palette on a grid (a dot plot's cells) colors each value as its legend says (VS-03)", () => {
    const set = SeqSets.build({
      molecules: [{ key: 's', residues: 'MK' }],
      units: [],
      alignments: [],
      annotations: [{ key: 'g', geometry: 'grid', target: 's', columns: { b_target: ['s'], value: [[1, 2, 3, 1]] } }],
    } as unknown as SetInput).set;
    const x = context(set, { data: ['g'] });
    const palette = { kind: 'categorical', colors: ['red', 'green', 'blue'] };
    const th = x.color('field_palette', { field: 'value', palette });
    const rgb = [C('red'), C('green'), C('blue')];
    expect([1, 2, 3].map((value) => th.value(at({ value })))).toEqual(rgb);
    expect(legendOf(th)).toEqual(['1', '2', '3'].map((k, i) => ['item', k, rgb[i], NO_GLYPH]));
    // A grid is numeric: the palette's own default colors by the value, as on `values` and `matrix`.
    expect(x.color('field_palette').granularity).toBe(Granularity.value);
  });

  it('a streamed item with no category column is colored by its id, as a held row is (VS-16)', () => {
    const ids = SeqSets.build({
      molecules: [{ key: 'p', residues: 'ACDEFGHIKL' }],
      units: [{ key: 'A', molecule: 'p' }],
      annotations: [
        { key: 'g', geometry: 'ranges', target: 'p', columns: { start: [0, 2], stop: [1, 3], id: ['a', 'b'] } },
      ],
    }).set;
    const th = context(ids, { data: ['g'] }).color('field_palette');
    const g = idOf(ids, 'g');
    expect(th.raw?.({ id: ['a', 'b'] }, 1)).toBe(th.value(at({ annotation: g, row: 1 })));
    expect(th.raw?.({ id: ['a', 'b'] }, 1)).not.toBe(NO_COLOR);
  });

  it("a streamed annotation with no category: its manifest's vocabulary of ids colors them (LB-08)", () => {
    const streamed = (vocabulary: Record<string, string[]>) => {
      const set = SeqSets.build({
        molecules: [{ key: 'p', residues: 'ACDEFGHIKL' }],
        units: [{ key: 'A', molecule: 'p' }],
        annotations: [{ key: 'g', geometry: 'ranges', streamed: ['p'], columns: {} }],
      } as never).set;
      set.annotations[0]?.custom.props.set('svs:stream', { manifest: { vocabulary } });
      return context(set, { data: ['g'] }).color('field_palette');
    };
    const byId = streamed({ id: ['a', 'b'] });
    const [a, b] = [byId.raw?.({ id: ['a', 'b'] }, 0), byId.raw?.({ id: ['a', 'b'] }, 1)];
    expect([a === NO_COLOR, b === NO_COLOR, a === b]).toEqual([false, false, false]);
    // With a category vocabulary the category colors, as before; a label vocabulary colors labels.
    const byCategory = streamed({ category: ['x', 'y'], id: ['a', 'b'] });
    expect(byCategory.raw?.({ category: ['x', 'y'] }, 1)).not.toBe(NO_COLOR);
    expect(streamed({ label: ['L1'] }).raw?.({ label: ['L1'] }, 0)).not.toBe(NO_COLOR);
  });

  it('a streamed annotation with no vocabulary: each value a stable color by its hash, one note (LB-08)', () => {
    const open = () => {
      const set = SeqSets.build({
        molecules: [{ key: 'p', residues: 'ACDEFGHIKL' }],
        units: [{ key: 'A', molecule: 'p' }],
        annotations: [{ key: 'g', geometry: 'ranges', streamed: ['p'], columns: {} }],
      } as never).set;
      const x = context(set, { data: ['g'] });
      return { x, th: x.color('field_palette') };
    };
    const { x, th } = open();
    // Two regions' tables: "a" in both, at other rows; ids by FNV-1a into Category10 ("a" 0, "b" 7, "c" 8).
    const one = { id: ['a', 'b'] };
    const two = { id: ['c', 'x', 'a'] };
    const [a1, b, c, a2] = [th.raw?.(one, 0), th.raw?.(one, 1), th.raw?.(two, 0), th.raw?.(two, 2)];
    expect([a1, b, c]).toEqual(['#1f77b4', '#7f7f7f', '#bcbd22'].map(C));
    expect(a2).toBe(a1);
    // A reload (a new context and theme): the same colors.
    const again = open().th;
    expect([again.raw?.(one, 0), again.raw?.(two, 0)]).toEqual([a1, c]);
    // The legend lists none of them, and one note says so.
    expect(legendOf(th)).toEqual([]);
    expect(x.issues.map((i) => [i.code, i.severity])).toEqual([['svs.legend.unlisted', 'info']]);
  });

  it("categorical: the palette's keys in order, with their colors", () => {
    const set = SeqSets.build(INPUT).set;
    const th = context(set, { data: ['f'] }).color('field_palette');
    expect(legendOf(th)).toEqual(
      ['x', 'y', 'z'].map((k, i) => ['item', k, C(['#1f77b4', '#ff7f0e', '#2ca02c'][i]), NO_GLYPH]),
    );
  });
});

describe('locations and ink', () => {
  it("`ThemeLocation.at` fills the block's unit, so `by_unit` colors features and positions", async () => {
    const tt = await opened(case3());
    const m = tt.ctx.graph.assembled.value?.view;
    const axis = tt.ctx.layout.axis(tt.ctx.viewports.sections.value[0])?.axis;
    if (!m || !axis) throw new Error('no view');
    const l = ThemeLocation.at(ThemeLocation.create(), axis, 2, 0, 17);
    expect([l.layer, l.block, l.unit, l.index, l.alignment >= 0, l.column]).toEqual([2, 0, axis.unit[0], 17, true, 17]);
  });
  it('dark ink on light cells, light ink on dark ones, whatever the skin', () => {
    expect([0xffff00, 0xe0e0e0, 0x80a0f0].map(cellInk)).toEqual([CellInk.dark, CellInk.dark, CellInk.dark]);
    expect([0x000080, 0x202020, 0x8b0000].map(cellInk)).toEqual([CellInk.light, CellInk.light, CellInk.light]);
    expect(cellInk(0xffff00)).toBe(CellInk.dark);
  });
});

describe('the theme node', () => {
  it("case 3: the rows' conservation shading is a compiled `rules` theme per channel, fading by the column", async () => {
    const tt = await opened(case3());
    const r = layer(tt.ctx, 'Members');
    const [color] = r.m.prepared[r.k].themes as (Theme | null)[];
    expect(r.m.resolved.layers.themes[r.k][0]).toMatchObject({ theme: 'rules', table: 'conservation_shading' });
    expect(color?.kind).toBe('color');
    expect(color?.granularity).toBe(Granularity.column);
    const a = r.m.set.alignments[r.m.resolved.blocks.alignment[0]];
    let faded = 0;
    for (let column = 0; a && column < a.columnCount; column++) {
      const alpha = color?.alpha?.(at({ block: 0, column, letter: code('L') })) ?? OPAQUE;
      if (alpha > 0 && alpha < OPAQUE) faded++;
    }
    expect(faded).toBeGreaterThan(10);
  });
});
