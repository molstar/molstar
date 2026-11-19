/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The theme entry: palettes (properties: a continuous ramp monotone in lightness, log monotone, bins closed below, the
// hash's known vectors and its stability, a dictionary code's color its value's), forms compiled to channels, schemes
// (fit, substitution, dark variants), legends (merge, Mol*'s types) and MolViewSpec's palettes both ways.

import { categoryColumn } from '@molstar/core/repr/tree';
import { Color } from '@molstar/core/util/color/color';
import {
  hashValueUtf8,
  makeValuePaletteFunction,
  valuePalettePropsFromSpec,
} from '@molstar/core/util/color/value-palette';
import { Oklab } from '@molstar/core/util/color/spaces/oklab';
import { DarkTokens, LightTokens, RowDatum } from '@molstar/seqx-core/render';
import { describe, expect, it } from 'vitest';
import {
  builtInSchemes,
  type Channel,
  compileChannel,
  compilePalette,
  contrastRgba,
  fromMvsPalette,
  inkAgainst,
  legendOf,
  mergeLegends,
  niceTicks,
  resolveColor,
  type ThemeBody,
  type ThemeDatum,
  type ThemeEnv,
  type ThemeNode,
  toMolLegend,
  toMvsPalette,
  toRgb,
} from '../index.js';

const L = (c: number) => Oklab.lightness(toRgb(c));
const schemes = builtInSchemes();

function env(over: Partial<ThemeEnv> = {}): ThemeEnv & { issues: string[] } {
  const issues: string[] = [];
  return { tokens: LightTokens, schemes, report: (code) => issues.push(code), issues, ...over };
}

const datum = (over: Partial<ThemeDatum> & { fields?: Record<string, string | number> } = {}): ThemeDatum => ({
  value: Number.NaN,
  letter: '',
  start: 0,
  end: 1,
  ...over,
  field: (n) => over.fields?.[n],
});

const color = (body: ThemeBody | undefined, e = env(), fallback: ThemeBody = 'token:range'): Channel =>
  compileChannel('color', body, fallback, e);

describe('palettes: properties', () => {
  it('a continuous ramp is monotone in lightness in both tables (the token ramp and its LUT)', () => {
    for (const tokens of [LightTokens, DarkTokens]) {
      const p = compilePalette({ kind: 'continuous', params: { colors: 'ramp' } }, { tokens, domain: [0, 100] });
      const ls = Array.from({ length: 21 }, (_, i) => L(p.num(i * 5)!));
      const up = ls[20] > ls[0];
      for (let i = 1; i < ls.length; i++) expect(up ? ls[i] >= ls[i - 1] - 1e-3 : ls[i] <= ls[i - 1] + 1e-3).toBe(true);
    }
  });

  it('a log scale is monotone, and a value at or below zero is missing', () => {
    const p = compilePalette(
      { kind: 'continuous', params: { colors: 'ramp' } },
      {
        tokens: LightTokens,
        domain: [1, 1e5],
        scale: { type: 'log' },
      },
    );
    const ls = [1, 10, 100, 1e3, 1e4, 1e5].map((x) => L(p.num(x)!));
    for (let i = 1; i < ls.length; i++) expect(ls[i]).toBeLessThan(ls[i - 1]);
    expect(p.num(0)).toBeUndefined();
    expect(p.num(-3)).toBeUndefined();
    // 10^2.5 sits at the middle of the ramp
    expect(p.num(10 ** 2.5)).toBe(p.num(10 ** 2.5));
  });

  it('bins are closed below and open above, the last closed; under and over take the ends', () => {
    const p = compilePalette(
      { kind: 'binned', params: { edges: [0, 50, 70, 90, 100], colors: ['#ff0000', '#00ff00', '#0000ff', '#000000'] } },
      { tokens: LightTokens },
    );
    expect(p.num(0)).toBe(0xff0000ff);
    expect(p.num(49.99)).toBe(0xff0000ff);
    expect(p.num(50)).toBe(0x00ff00ff);
    expect(p.num(90)).toBe(0x000000ff);
    expect(p.num(100)).toBe(0x000000ff);
    expect(p.num(-1)).toBe(0xff0000ff);
    expect(p.num(101)).toBe(0x000000ff);
    expect(p.num(Number.NaN)).toBeUndefined();
  });

  it('an unlisted value takes the color at its FNV-1a hash over UTF-8, the same in every palette instance', () => {
    const cats = LightTokens.colors;
    const list = Array.from({ length: 8 }, (_, i) => cats[`seq-cat-${i}` as keyof typeof cats]);
    const a = compilePalette({ kind: 'categorical', params: {} }, { tokens: LightTokens });
    const b = compilePalette({ kind: 'categorical', params: { values: ['x', 'y'] } }, { tokens: LightTokens });
    for (const v of ['foobar', 'Ω', 'kinase', '']) {
      expect(a.color(v)).toBe(list[hashValueUtf8(v) % 8]);
      // another reader (a new palette, another tile) agrees
      expect(compilePalette({ kind: 'categorical', params: {} }, { tokens: LightTokens }).color(v)).toBe(a.color(v));
    }
    expect(b.color('x')).toBe(list[0]);
    expect(b.color('y')).toBe(list[1]);
  });

  it('a dictionary-coded value’s color is its value’s', () => {
    const column = categoryColumn(['b', 'a', 'b', 'c', 'a']);
    const d = new RowDatum();
    d.bind([{ name: 'category', column }]);
    const ch = color({ kind: 'field', params: { name: 'category' } });
    for (let i = 0; i < 5; i++) {
      d.row = i;
      expect(ch.at(d)).toBe(ch.at(datum({ fields: { category: column.get(i)! } })));
    }
  });

  it('the LUT agrees with Mol*’s value palette at each stop (Oklab between, so only the stops are equal)', () => {
    const stops = ['#0d0887', '#cc4778', '#f0f921'];
    const p = compilePalette(
      { kind: 'continuous', params: { colors: stops } },
      { tokens: LightTokens, domain: [0, 1] },
    );
    const mol = makeValuePaletteFunction(
      valuePalettePropsFromSpec(
        { kind: 'continuous', colors: stops },
        {
          categorical: {
            kind: 'categorical',
            colors: [],
            repeat_color_list: false,
            sort: 'none',
            sort_direction: 'ascending',
            case_insensitive: false,
            missing_color: null,
          },
          discrete: { kind: 'discrete', colors: [], reverse: false, mode: 'normalized', value_domain: [null, null] },
          continuous: {
            kind: 'continuous',
            colors: [],
            reverse: false,
            mode: 'normalized',
            value_domain: [0, 1],
            underflow_color: null,
            overflow_color: null,
          },
        },
      ),
      [],
    );
    // the middle stop falls between two of the 256 entries: within a step there
    for (const x of [0, 1]) expect(toRgb(p.num(x)!)).toBe(mol(x));
    const [a, b] = [toRgb(p.num(0.5)!), mol(0.5)!];
    for (const shift of [16, 8, 0])
      expect(Math.abs(((a >> shift) & 255) - ((b >> shift) & 255))).toBeLessThanOrEqual(3);
  });
});

describe('forms', () => {
  it('a string is uniform (a constant), a token follows the table, light-dark picks by theme', () => {
    expect(color('#d62728').constant).toBe(0xd62728ff);
    expect(color('token:range').constant).toBe(LightTokens.colors['seq-range']);
    expect(color('light-dark(#000000, #ffffff)').constant).toBe(0x000000ff);
    expect(color('light-dark(#000000, #ffffff)', env({ tokens: DarkTokens })).constant).toBe(0xffffffff);
    expect(compileChannel('size', 0.5, 1, env()).constant).toBe(0.5);
    const e = env();
    expect(compileChannel('line', 'squiggly', 'solid', e).constant).toBeUndefined();
    expect(e.issues).toEqual(['seq.render.line-style']);
  });

  it('rules: the first that holds gives its value; a comparison across types never holds; else, then the default', () => {
    const plddt = (v: number | string) => datum({ fields: { plddt: v } });
    const rules: ThemeNode = {
      kind: 'rules',
      params: { else: '#cccccc' },
      args: [
        {
          kind: 'rule',
          params: {
            when: [{ kind: 'ge', args: [{ kind: 'field', params: { name: 'plddt' } }, 90] }],
            value: '#0053d6',
          },
        },
        {
          kind: 'rule',
          params: {
            when: [{ kind: 'ge', args: [{ kind: 'field', params: { name: 'plddt' } }, 70] }],
            value: '#65cbf3',
          },
        },
        {
          kind: 'rule',
          params: {
            when: [{ kind: 'eq', args: [{ kind: 'field', params: { name: 'plddt' } }, '50'] }],
            value: '#ff0000',
          },
        },
      ],
    };
    const ch = color(rules);
    expect(ch.at(plddt(95))).toBe(0x0053d6ff);
    expect(ch.at(plddt(75))).toBe(0x65cbf3ff);
    expect(ch.at(plddt(50))).toBe(0xccccccff);
    expect(ch.at(plddt('50'))).toBe(0xff0000ff);
    expect(ch.at(plddt(Number.NaN))).toBe(0xccccccff);
    expect(ch.legend).toMatchObject({
      kind: 'swatches',
      items: [{ label: 'plddt ≥ 90' }, { label: 'plddt ≥ 70' }, { label: 'plddt = 50' }, { label: 'otherwise' }],
    });
    // no else: the layer type's default
    const { params: _, ...noElse } = rules;
    expect(color(noElse).at(plddt(10))).toBe(LightTokens.colors['seq-range']);
  });

  it('a scheme that does not fit draws the default and says what it needs', () => {
    const e = env({ alphabet: 'protein' });
    const ch = compileChannel('color', { kind: 'scheme', params: { name: 'clustal-x' } }, 'token:range', e);
    expect(ch.constant).toBe(LightTokens.colors['seq-range']);
    expect(e.issues).toEqual(['seq.render.theme.scheme']);
    const n = env({ alphabet: 'dna' });
    compileChannel('color', { kind: 'scheme', params: { name: 'zappo' } }, 'token:range', n);
    expect(n.issues).toEqual(['seq.render.theme.scheme']);
  });

  it('letters: colors first, then groups, then the default; the reader’s substitution draws a safe scheme', () => {
    const zappo = color({ kind: 'scheme', params: { name: 'zappo' } }, env({ alphabet: 'protein' }));
    expect(zappo.at(datum({ letter: 'C' }))).toBe(0xffff00ff);
    expect(zappo.at(datum({ letter: 'c' }))).toBe(0xffff00ff);
    expect(zappo.at(datum({ letter: 'X' }))).toBe(LightTokens.colors['seq-missing']);
    const safe = color({ kind: 'scheme', params: { name: 'zappo' } }, env({ alphabet: 'protein', cvdSafe: true }));
    expect(safe.scheme).toEqual({ name: 'classes-safe', substituted: 'zappo' });
    expect(safe.at(datum({ letter: 'C' }))).toBe(0x56b4e9ff);
    const classes = color({ kind: 'scheme', params: { name: 'classes' } }, env({ alphabet: 'protein', cvdSafe: true }));
    expect(classes.scheme).toEqual({ name: 'classes-safe', substituted: 'classes' });
    const over = color(
      { kind: 'scheme', params: { name: 'zappo', overrides: { C: '#123456' } } },
      env({ alphabet: 'protein' }),
    );
    expect(over.at(datum({ letter: 'C' }))).toBe(0x123456ff);
  });

  it('a scheme’s literal colors take a derived dark variant: hue kept, lightness in the dark band', () => {
    const light = color({ kind: 'scheme', params: { name: 'taylor' } }, env({ alphabet: 'protein' }));
    const dark = color(
      { kind: 'scheme', params: { name: 'taylor' } },
      env({ alphabet: 'protein', tokens: DarkTokens }),
    );
    for (const l of 'ACDEFGHIKLMNPQRSTVWY') {
      const a = light.at(datum({ letter: l })) as number;
      const b = dark.at(datum({ letter: l })) as number;
      expect(L(b)).toBeLessThanOrEqual(0.63);
      expect(L(b)).toBeLessThan(L(a));
    }
  });

  it('an inline open domain is held by key and read back; a shared scale by its address', () => {
    const held = new Map<string, readonly [number, number]>([['/views/s', [0, 10]]]);
    const ch = color(
      { kind: 'field', params: { name: 'score', scale: '/views/s' } },
      env({
        scale: () => ({ palette: { kind: 'continuous', params: { colors: ['#000000', '#ffffff'] } } }),
        domain: (k) => held.get(k),
      }),
    );
    expect(ch.open.map((o) => o.key)).toEqual(['/views/s']);
    expect(ch.at(datum({ fields: { score: 10 } }))).toBe(0xffffffff);
    expect(ch.at(datum({ fields: { score: 0 } }))).toBe(0x000000ff);
  });

  it('a size channel maps through its scale’s range; a radius reads sqrt', () => {
    const size = compileChannel(
      'size',
      { kind: 'field', params: { name: 'n', domain: [0, 100], range: [2, 12] } },
      4,
      env({ radius: true }),
    );
    expect(size.at(datum({ fields: { n: 0 } }))).toBe(2);
    expect(size.at(datum({ fields: { n: 25 } }))).toBe(7);
    expect(size.at(datum({ fields: { n: 100 } }))).toBe(12);
  });

  it('a themed ink holds 3:1 against its fill', () => {
    for (const [ink, fill] of [
      [0x777777ff, 0x888888ff],
      [0xffff00ff, 0xffffffff],
      [0x0000ffff, 0x000033ff],
    ])
      expect(contrastRgba(inkAgainst(ink, fill), fill)).toBeGreaterThanOrEqual(3);
  });
});

describe('legends', () => {
  const field = (name: string, palette?: ThemeNode) =>
    color({ kind: 'field', params: { name, ...(palette ? { palette } : {}) } }, env({ numeric: (n) => n === 'score' }));

  it('derives swatches (vocabulary or painted), a ramp with ticks, steps; says hashed values once', () => {
    const listed = legendOf(field('kind', { kind: 'categorical', params: { values: ['a', 'b'] } }), {
      title: 'Kinds',
      layer: 'l1',
    })!;
    expect(listed.source).toBe('vocabulary');
    expect(listed.notes).toEqual([]);
    expect(listed.entries).toMatchObject({ kind: 'swatches', items: [{ label: 'a' }, { label: 'b' }] });
    const painted = legendOf(field('kind'), { title: 'Kinds', layer: 'l1', painted: ['x', 'y'] })!;
    expect(painted.source).toBe('painted');
    expect(painted.notes).toEqual(['unlisted']);
    const ramp = legendOf(field('score'), { title: 'Score', layer: 'l2' })!;
    expect(ramp.entries.kind).toBe('ramp');
    expect(ramp.open).toBe(true);
    expect(niceTicks(0, 1)).toEqual([0, 0.2, 0.4, 0.6, 0.8, 1]);
    expect(legendOf(color('#ff0000'), { title: 'x', layer: 'l' })).toBeUndefined();
  });

  it('merges legends that say the same, titled by what their titles share; converts to Mol*’s legends', () => {
    const ch = field('kind', { kind: 'categorical', params: { values: ['a', 'b'] } });
    const merged = mergeLegends([
      legendOf(ch, { title: 'Domains of p53', layer: 'l1' })!,
      legendOf(ch, { title: 'Domains of p63', layer: 'l2' })!,
    ]);
    expect(merged).toHaveLength(1);
    expect(merged[0].title).toBe('Domains of');
    expect(merged[0].layers).toEqual(['l1', 'l2']);
    expect(toMolLegend(merged[0])).toMatchObject({
      kind: 'table-legend',
      table: [
        ['a', expect.any(Number)],
        ['b', expect.any(Number)],
      ],
    });
    const ramp = legendOf(
      color({ kind: 'field', params: { name: 'v', domain: [0, 10] } }, env({ numeric: () => true })),
      { title: 'v', layer: 'l' },
    )!;
    expect(toMolLegend(ramp)).toMatchObject({ kind: 'scale-legend', minLabel: '0', maxLabel: '10' });
  });
});

describe('MolViewSpec’s palettes', () => {
  it('converts both ways, a value having one color in each', () => {
    const node: ThemeNode = { kind: 'binned', params: { edges: [0, 50, 100], colors: ['#ff0000', '#0000ff'] } };
    const mvs = toMvsPalette(node, LightTokens)!;
    expect(mvs.kind).toBe('discrete');
    const back = fromMvsPalette(mvs);
    expect(back.kind).toBe('binned');
    const ours = compilePalette(back, { tokens: LightTokens });
    expect(ours.num(25)).toBe(0xff0000ff);
    expect(ours.num(75)).toBe(0x0000ffff);
    const cat = toMvsPalette({ kind: 'categorical', params: { colors: { a: 'token:range' } } }, LightTokens)!;
    expect(cat).toMatchObject({ kind: 'categorical', unlisted: 'hash' });
    expect(resolveColor((cat as { colors: Record<string, string> }).colors.a, LightTokens)).toBe(
      LightTokens.colors['seq-range'],
    );
    expect(fromMvsPalette({ kind: 'continuous', colors: 'viridis', mode: 'absolute' })).toEqual({
      kind: 'continuous',
      params: { colors: 'viridis', stops: 'absolute' },
    });
    // the core's list order is Mol*'s: a value is the same color in MVS and here
    expect(
      Color.toHexStyle(
        toRgb(
          compilePalette(fromMvsPalette({ kind: 'categorical', colors: { q: '#112233' } }), {
            tokens: LightTokens,
          }).color('q')!,
        ),
      ),
    ).toBe('#112233');
  });
});
