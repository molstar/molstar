/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The theme sublanguage, `seq.theme` (_spec/README.md "Themes", Forms), described as SeqViewSpec is: the forms a
// channel slot holds, their palettes and rules' clauses. Its nodes sit in slots (a layer's channels, a scale's
// palette) and have no ids: a slot's body is replaced whole. What the description cannot say (a color in a size, a
// comparison across types, nesting past `themeDepth`) is the language's checks (`seq.lang.theme.*`).

import { defineLanguage, type ParamType } from '@molstar/core/repr/schema';

export const ThemeLanguageId = 'seq.theme';

const node = (kinds?: readonly string[]): ParamType => ({
  type: 'node',
  language: ThemeLanguageId,
  ...(kinds ? { kinds } : {}),
});

/** A channel's value: a color, a number, a word, or a nested form. */
const channelValue: ParamType = { type: 'union', of: [{ type: 'string' }, { type: 'number' }, node()] };

const PaletteKinds = ['categorical', 'direct', 'binned', 'continuous'];
const MeasureKinds = ['field', 'value', 'letter', 'share', 'identity', 'matches-reference'];

const colors: ParamType = {
  type: 'union',
  of: [
    { type: 'string' },
    { type: 'list', of: { type: 'union', of: [{ type: 'string' }, { type: 'list', of: { type: 'any' } }] } },
    { type: 'map', of: { type: 'string' } },
  ],
};

/** A domain's end: a number, or `null` (open). */
const endpoint: ParamType = { type: 'any' };

/** The scale params a `field` form takes inline and a `scale` node holds. */
export const ScaleParamsDescription = {
  type: {
    type: 'string',
    enum: ['linear', 'log', 'sqrt'],
    doc: 'How values are read: `linear` (when absent; a radius: `sqrt`), `log`, `sqrt`.',
  },
  domain: {
    type: 'list',
    of: endpoint,
    'min-length': 2,
    'max-length': 2,
    doc: '`[a, b]`; an end `null` (or absent) is open: the union over the layers’ data in view, at settle.',
  },
  palette: { ...node(PaletteKinds), doc: 'A color scale’s palette.' },
  range: {
    type: 'list',
    of: { type: 'number', doc: 'An end.' },
    'min-length': 2,
    'max-length': 2,
    doc: 'A size scale’s output, `[from, to]` in the channel’s unit.',
  },
  mid: { type: 'number', doc: 'A diverging center.' },
  nice: { type: 'boolean', doc: 'Legend ticks at round numbers.' },
  label: { type: 'string', doc: 'A label for people (a legend’s title).' },
  unit: { type: 'string', doc: 'The values’ unit, for people.' },
  non_positive: {
    type: 'string',
    enum: ['missing', 'clamp'],
    doc: 'On a log scale, a value at or below zero: `missing` (when absent) or `clamp` to the smallest positive.',
  },
} as const;

const clause = (doc: string, rest: boolean) => ({
  doc,
  parents: [] as string[],
  data: {
    args: [
      {
        name: 'measure',
        ...node(),
        doc: 'What it reads: `(field :name n)`, `(value)`, `(letter)`, a registered measure.',
      },
      { name: 'value', type: 'any' as const, doc: 'What it compares with.' },
    ],
    ...(rest ? { rest: true, max: 256 } : { min: 2 }),
  },
});

export const ThemeLanguage = defineLanguage({
  id: ThemeLanguageId,
  version: '0.4.0',
  root: 'uniform',
  doc: 'The theme sublanguage: what colors and sizes a mark, as data in a channel slot.',
  kinds: {
    uniform: {
      doc: 'One value for every mark (a plain string or number in a slot is this form).',
      parents: [],
      ids: 'forbidden',
      params: { value: { ...channelValue, required: true, doc: 'The value.' } },
    },
    field: {
      doc: 'A field’s value through a scale: inline (`palette`, `domain`, `type`, `mid`, `range`) or a `scale` node by address.',
      parents: [],
      ids: 'forbidden',
      params: {
        name: { type: 'string', doc: 'The column (`value`: a values layer’s value; `letter`: the cell’s letter).' },
        data: { type: 'ref', kinds: ['values', 'ranges'], doc: 'Another data node, read per position.' },
        scale: { type: 'ref', kinds: ['scale'], doc: 'A scale node: layers naming one share its domain and palette.' },
        ...ScaleParamsDescription,
      },
    },
    rules: {
      doc: 'The first rule whose clauses all hold gives its value; else `else` (absent: the layer type’s default).',
      parents: [],
      ids: 'forbidden',
      params: { else: { ...channelValue, doc: 'The value where no rule holds.' } },
      data: { args: [{ name: 'rule', ...node(['rule']), doc: 'A rule.' }], rest: true, min: 0, max: 256 },
    },
    rule: {
      doc: 'Clauses (an `and`) and the value they give.',
      parents: [],
      ids: 'forbidden',
      params: {
        when: { type: 'list', of: node(), required: true, doc: 'The clauses, all of which hold.' },
        value: { ...channelValue, required: true, doc: 'A channel value or a nested form.' },
        label: { type: 'string', doc: 'A legend label (its clauses when absent).' },
      },
    },
    scheme: {
      doc: 'A named scheme from the registry, with overrides.',
      parents: [],
      ids: 'forbidden',
      params: {
        name: { type: 'string', required: true, doc: 'The scheme’s name.' },
        overrides: { type: 'map', of: { type: 'string' }, doc: 'Colors by value or letter over the scheme’s.' },
        field: { type: 'string', doc: 'The field a numbers or categories scheme reads.' },
        mid: { type: 'number', doc: 'The `diverging` scheme’s center.' },
        gaps: { type: 'string', enum: ['count', 'ignore'], doc: 'ESPript: gaps counted against identity (`count`).' },
        similarity: { type: 'number', min: 0, max: 1, doc: 'ESPript: the share a similar column holds (0.7).' },
        mode: { type: 'string', enum: ['color', 'bw'], doc: 'ESPript: `bw` is its black-and-white figure.' },
      },
    },
    letters: {
      doc: 'A value by the letter: `colors` first, then the groups, then `default`.',
      parents: [],
      ids: 'forbidden',
      params: {
        colors: { type: 'map', of: { type: 'string' }, doc: 'A color by letter.' },
        groups: { type: 'list', of: node(['group']), doc: 'Groups of letters with a value and a label.' },
        default: { type: 'string', doc: 'The value of any other letter (absent: the layer’s default).' },
      },
    },
    group: {
      doc: 'Letters that take one value, with a label.',
      parents: [],
      ids: 'forbidden',
      params: {
        letters: { type: 'string', required: true, doc: 'The letters.' },
        value: { type: 'string', required: true, doc: 'Their value.' },
        label: { type: 'string', doc: 'A label for people.' },
      },
    },
    categorical: {
      doc: 'A color per value: a list (by `values`’ order), a named list, or a dictionary.',
      parents: [],
      ids: 'forbidden',
      params: {
        colors: { ...colors, doc: 'A named list, colors, or `{value: color}`.' },
        values: { type: 'list', of: { type: 'string' }, doc: 'The listed values, in order.' },
        repeat: { type: 'boolean', doc: 'The list repeats past its end (true).' },
        unlisted: {
          type: 'string',
          enum: ['hash', 'missing'],
          doc: 'An unlisted value: `hash` (FNV-1a over UTF-8, the default), `missing`.',
        },
        others: { ...colors, doc: 'The list an unlisted value hashes into (a dictionary’s: `categories`).' },
      },
    },
    direct: { doc: 'The field holds colors.', parents: [], ids: 'forbidden' },
    binned: {
      doc: 'Bins: `edges` (n + 1, closed below, the last closed) and `colors` (n), or `bins` n over the domain.',
      parents: [],
      ids: 'forbidden',
      params: {
        edges: { type: 'list', of: { type: 'number' }, 'min-length': 2, doc: 'The edges, ascending.' },
        colors: { ...colors, doc: 'A color a bin, or a named list.' },
        bins: { type: 'integer', min: 1, doc: 'Bins over the domain, evenly through the scale.' },
        under: { type: 'string', doc: 'Below the first edge.' },
        over: { type: 'string', doc: 'Above the last edge.' },
        mid: { type: 'number', doc: 'A diverging center.' },
      },
    },
    continuous: {
      doc: 'Colors interpolated in Oklab: a named list or stops `[color, at]` (`at` 0..1, or values with `stops: absolute`).',
      parents: [],
      ids: 'forbidden',
      params: {
        colors: { ...colors, doc: 'A named list or stops.' },
        stops: { type: 'string', enum: ['relative', 'absolute'], doc: '`absolute` reads `at` as values.' },
        mid: { type: 'number', doc: 'A diverging center.' },
        under: { type: 'string', doc: 'Below the domain (the low end’s color).' },
        over: { type: 'string', doc: 'Above the domain (the high end’s color).' },
      },
    },
    eq: clause('The measure equals the value (text with text, numbers with numbers).', false),
    ne: clause('The measure differs from the value (of the same type).', false),
    lt: clause('The measure is below the number.', false),
    le: clause('The measure is at or below the number.', false),
    gt: clause('The measure is above the number.', false),
    ge: clause('The measure is at or above the number.', false),
    'one-of': clause('The measure equals one of the values.', true),
    in: {
      doc: 'The mark overlaps a location by one position or more.',
      parents: [],
      data: { args: [{ name: 'location', type: 'ref', doc: 'An address of a location.' }] },
    },
    within: {
      doc: 'The mark’s extent lies inside a location.',
      parents: [],
      data: { args: [{ name: 'location', type: 'ref', doc: 'An address of a location.' }] },
    },
    value: { doc: 'The measure of a values layer’s value.', parents: [], ids: 'forbidden' },
    share: {
      doc: 'A column measure (the rows layer’s): the share of the column’s rows with the cell’s letter (or its group).',
      parents: [],
      ids: 'forbidden',
      params: {
        gaps: { type: 'string', enum: ['count', 'ignore'], doc: 'Gaps counted against the share (`count`).' },
        by: { type: 'string', enum: ['letter', 'group'], doc: 'The letter itself (`letter`) or its residue group.' },
      },
    },
    identity: { doc: 'A column measure: the column’s identity, 0..1.', parents: [], ids: 'forbidden' },
    'matches-reference': {
      doc: 'A column measure: 1 where the cell’s letter is the reference row’s.',
      parents: [],
      ids: 'forbidden',
    },
    letter: { doc: 'The measure of a cell’s letter.', parents: [], ids: 'forbidden' },
  },
});

/** The kinds a clause's measure may be (a registered measure's kind besides). */
export const ThemeMeasureKinds: readonly string[] = MeasureKinds;
