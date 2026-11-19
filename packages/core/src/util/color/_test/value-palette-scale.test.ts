/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Color } from '../color.js';
import {
  hashValueUtf8,
  makeValuePaletteFunction,
  type ValuePaletteFunction,
  ValuePaletteLists,
  type ValuePaletteSpec,
  type ValuePaletteSpecDefaults,
  valuePalettePropsFromSpec,
} from '../value-palette.js';

// The optional fields: `scale`, `non_positive`, `mid` on the numeric palettes, `unlisted: 'hash'` and `others` on the
// categorical one, and the lazy lists. MVS's defaults (its records as they were) leave them out.
const Defaults: ValuePaletteSpecDefaults = {
  categorical: {
    kind: 'categorical',
    colors: 'Category10',
    repeat_color_list: false,
    sort: 'none',
    sort_direction: 'ascending',
    case_insensitive: false,
    missing_color: null,
  },
  discrete: { kind: 'discrete', colors: 'YlGn', reverse: false, mode: 'normalized', value_domain: [null, null] },
  continuous: {
    kind: 'continuous',
    colors: 'YlGn',
    reverse: false,
    mode: 'normalized',
    value_domain: [null, null],
    underflow_color: null,
    overflow_color: null,
  },
};

const hex = (c: Color | undefined) => (c === undefined ? undefined : Color.toHexStyle(c));

function palette(spec: ValuePaletteSpec, values: string[] = []): ValuePaletteFunction {
  return makeValuePaletteFunction(valuePalettePropsFromSpec(spec, Defaults), values);
}

const grey = (spec: Partial<ValuePaletteSpec> = {}) =>
  palette({ kind: 'continuous', colors: ['#000000', '#ffffff'], ...spec } as ValuePaletteSpec);

describe('value palettes: the optional fields', () => {
  it('converts a list when read, once, to the colors it held', () => {
    const a = ValuePaletteLists.Viridis;
    expect(ValuePaletteLists.Viridis).toBe(a);
    expect(a.length).toBeGreaterThan(10);
    expect(Object.keys(ValuePaletteLists)).toContain('Tableau10');
  });

  it('reads a log scale by the logarithm: non-positive missing, or clamped to the low end; monotone', () => {
    const spec = { value_domain: [1, 100] as [number, number], scale: 'log' as const };
    const f = grey(spec);
    expect(hex(f(1))).toBe('#000000');
    expect(hex(f(100))).toBe('#ffffff');
    expect(hex(f(10))).toBe(hex(grey()(0.5)));
    expect(f(0)).toBeUndefined();
    expect(f(-1)).toBeUndefined();
    expect(hex(grey({ ...spec, non_positive: 'clamp' })(0))).toBe('#000000');
    let prev = -1;
    for (const x of [1, 2, 5, 10, 20, 50, 100]) {
      const c = (f(x) as number) & 255;
      expect(c).toBeGreaterThanOrEqual(prev);
      prev = c;
    }
  });

  it('opens a log domain over the values through the scale; sqrt likewise', () => {
    const f = palette({ kind: 'continuous', colors: ['#000000', '#ffffff'], scale: 'log' }, ['1', '10', '100', '-5']);
    expect(hex(f(1))).toBe('#000000');
    expect(hex(f(100))).toBe('#ffffff');
    expect(hex(grey({ scale: 'sqrt', value_domain: [0, 100] })(25))).toBe(hex(grey()(0.5)));
  });

  it('reads a binned log palette’s absolute edges as values', () => {
    const f = palette({
      kind: 'discrete',
      mode: 'absolute',
      colors: [
        ['red', 1, 10],
        ['blue', 10, 1000],
      ],
      scale: 'log',
    });
    expect(hex(f(5))).toBe('#ff0000');
    expect(hex(f(500))).toBe('#0000ff');
  });

  it('puts mid at the middle color, each half scaled on its own', () => {
    const f = palette({
      kind: 'continuous',
      colors: ['#0000ff', '#ffffff', '#ff0000'],
      value_domain: [-1, 10],
      mid: 0,
    });
    expect(hex(f(0))).toBe('#ffffff');
    expect(hex(f(-1))).toBe('#0000ff');
    expect(hex(f(10))).toBe('#ff0000');
  });

  it('hashes an unlisted value with FNV-1a over UTF-8 into the list or `others`', () => {
    // FNV-1a's known vectors
    expect(hashValueUtf8('')).toBe(0x811c9dc5);
    expect(hashValueUtf8('a')).toBe(0xe40c292c);
    expect(hashValueUtf8('foobar')).toBe(0xbf9cf968);
    const list = ['#000001', '#000002', '#000003'];
    const f = palette({ kind: 'categorical', colors: list, unlisted: 'hash' }, ['x']);
    expect(hex(f('x'))).toBe('#000001');
    expect(hex(f('foobar'))).toBe(list[0xbf9cf968 % 3]);
    expect(hex(f('α'))).toBe(list[hashValueUtf8('α') % 3]);
    const g = palette({ kind: 'categorical', colors: { A: '#ff0000' }, unlisted: 'hash', others: list });
    expect(hex(g('A'))).toBe('#ff0000');
    expect(hex(g('foobar'))).toBe(list[0xbf9cf968 % 3]);
    // MVS's default unchanged: an unlisted value is missing
    expect(hex(palette({ kind: 'categorical', colors: { A: '#ff0000' }, missing_color: 'blue' })('B'))).toBe('#0000ff');
  });
});
