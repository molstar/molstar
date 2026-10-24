/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { ParamDefinition as PD } from '../../param-definition.js';
import { Color } from '../color.js';
import { ColorNames } from '../names.js';
import { SplitColor, SplitColorValue } from '../split-color.js';
import { decodeColor } from '../utils.js';
import {
  makeValuePaletteFunction,
  type ValuePaletteFunction,
  ValuePaletteLists,
  type ValuePaletteSpec,
  type ValuePaletteSpecDefaults,
  type ValuePaletteValue,
  valuePalettePropsFromSpec,
} from '../value-palette.js';
import { ValuePaletteParam } from '../value-palette-params.js';

/** MVS's defaults (`tree/mvs/param-types.ts`). */
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

function palette(spec: ValuePaletteSpec, values: string[] = [], cifMissing = false): ValuePaletteFunction {
  const props = valuePalettePropsFromSpec(spec, Defaults, { dictionaries: { Fruit: { apple: 'red', kiwi: '#0f0' } } });
  return makeValuePaletteFunction(props, values, { cifMissing });
}

describe('split colors', () => {
  it('decode both halves, or the one color twice', () => {
    expect(SplitColor.decode('red/white')).toEqual([ColorNames.red, ColorNames.white]);
    expect(SplitColor.decode('red')).toEqual([ColorNames.red, ColorNames.red]);
    expect(SplitColor.decode('nonsense')).toBeUndefined();
    expect(() => SplitColor.decodeStrict('nonsense')).toThrow();
    expect(SplitColor.toHexStyle(Color(0xff0000), Color(0xffffff))).toBe('#ff0000/#ffffff');
  });

  it('props: invalid halves are black, null is none', () => {
    expect(SplitColorValue.fromString('red/nonsense')).toEqual({
      name: 'splitColor',
      params: { color1: ColorNames.red, color2: ColorNames.black },
    });
    expect(SplitColorValue.fromString(null)).toEqual({ name: 'none', params: {} });
    expect(SplitColorValue.toTuple(SplitColorValue.fromString('blue'))).toEqual([ColorNames.blue, ColorNames.blue]);
  });

  it("names are only the table's own", () => {
    expect(decodeColor('constructor')).toBeUndefined();
    expect(decodeColor('RED')).toBe(ColorNames.red);
  });
});

describe('spec to props', () => {
  it('no palette is direct; the values are colors', () => {
    const props = valuePalettePropsFromSpec(undefined, Defaults);
    expect(props).toEqual({ name: 'direct', params: {} });
    const f = makeValuePaletteFunction(props, []);
    expect(hex(f('red/blue'))).toBe('#ff0000');
    expect(hex(f('red/blue', true))).toBe('#0000ff');
    expect(f('nonsense')).toBeUndefined();
  });

  it('names: lists, then dictionaries, unknown reported', () => {
    const onUnknownName = vi.fn();
    const list = valuePalettePropsFromSpec({ kind: 'categorical', colors: 'Set1' }, Defaults);
    expect(list.name === 'categorical' && list.params.colors.name).toBe('list');
    const dict = valuePalettePropsFromSpec({ kind: 'categorical', colors: 'Fruit' }, Defaults, {
      dictionaries: { Fruit: { apple: 'red' } },
    });
    expect(dict.name === 'categorical' && dict.params.colors).toEqual({
      name: 'dictionary',
      params: [{ value: 'apple', color: { name: 'oneColor', params: { color: ColorNames.red } } }],
    });
    valuePalettePropsFromSpec({ kind: 'continuous', colors: 'Fruit' }, Defaults, { onUnknownName });
    valuePalettePropsFromSpec({ kind: 'categorical', colors: 'toString' }, Defaults, { onUnknownName });
    expect(onUnknownName.mock.calls).toEqual([['Fruit'], ['toString']]);
  });

  it('`undefined` keeps a default, `null` is a value', () => {
    const props = valuePalettePropsFromSpec(
      { kind: 'continuous', colors: ['red', 'blue'], mode: undefined, overflow_color: null, value_domain: [0, null] },
      { ...Defaults, continuous: { ...Defaults.continuous, overflow_color: 'auto' } },
    );
    expect(props.name === 'continuous' && props.params).toMatchObject({ mode: 'normalized', xMin: 0, xMax: null });
    expect(props.name === 'continuous' && props.params.overflowColor).toEqual({ name: 'none', params: {} });
  });

  it("the named lists are Mol*'s, in MVS's names", () => {
    expect(Object.keys(ValuePaletteLists)).toHaveLength(50);
    expect(ValuePaletteLists.Category10[0]).toBe(0x1f77b4);
    expect(ValuePaletteLists.Chainbow.length).toBeGreaterThan(0);
  });

  it("the props are the param's values", () => {
    const param = ValuePaletteParam();
    expect(param.defaultValue).toEqual({ name: 'direct', params: {} });
    for (const spec of [
      { kind: 'categorical', colors: { a: 'red/white' }, missing_color: 'grey' },
      {
        kind: 'discrete',
        colors: [
          ['red', null, 1],
          [null, 1, null],
        ],
      },
      { kind: 'continuous', colors: 'Viridis', underflow_color: 'auto' },
    ] as ValuePaletteSpec[]) {
      const props = valuePalettePropsFromSpec(spec, Defaults);
      expect(PD.validate(PD.Group({ p: param }), { p: props }), spec.kind).toBeUndefined();
    }
  });
});

describe('categorical', () => {
  it('a list in order of first occurrence, then uncolored or repeated', () => {
    const values = ['b', 'a', 'b', 'c'];
    const f = palette({ kind: 'categorical', colors: ['red', 'green'] }, values);
    expect(['a', 'b', 'c'].map((v) => hex(f(v)))).toEqual(['#008000', '#ff0000', undefined]);
    expect(f.keys).toEqual(['b', 'a']);
    const g = palette({ kind: 'categorical', colors: ['red', 'green'], repeat_color_list: true }, values);
    expect(hex(g('c'))).toBe('#ff0000');
  });

  it('sort, direction, case and the missing color', () => {
    const f = palette(
      {
        kind: 'categorical',
        colors: ['red', 'green', 'blue'],
        sort: 'numeric',
        sort_direction: 'descending',
        case_insensitive: true,
        missing_color: 'black',
      },
      ['2', '10', '1'],
    );
    expect(f.keys).toEqual(['10', '2', '1']);
    expect(['10', '2', '1', '3'].map((v) => hex(f(v)))).toEqual(['#ff0000', '#008000', '#0000ff', '#000000']);
  });

  it('a dictionary, split colors by their halves', () => {
    const f = palette({ kind: 'categorical', colors: { A: 'red/white' }, case_insensitive: true }, ['a', 'b']);
    expect(hex(f('a'))).toBe('#ff0000');
    expect(hex(f('a', true))).toBe('#ffffff');
    expect(f('b')).toBeUndefined();
    expect(f.keys).toEqual(['A']);
  });

  it('`.` and `?` are values, or missing under `cifMissing`', () => {
    const spec: ValuePaletteSpec = { kind: 'categorical', colors: ['red', 'green'], missing_color: 'blue' };
    const values = ['.', 'a', '?'];
    expect(['.', 'a'].map((v) => hex(palette(spec, values)(v)))).toEqual(['#ff0000', '#008000']);
    const f = palette(spec, values, true);
    expect(['.', '?', 'a'].map((v) => hex(f(v)))).toEqual(['#0000ff', '#0000ff', '#ff0000']);
    expect(f.keys).toEqual(['a']);
  });

  it('lazy values: read only when needed', () => {
    const values = vi.fn(() => ['x']);
    const dict = valuePalettePropsFromSpec({ kind: 'categorical', colors: { x: 'red' } }, Defaults);
    const f = makeValuePaletteFunction(dict, values);
    expect(hex(f('x'))).toBe('#ff0000');
    expect(values).not.toHaveBeenCalled();
    expect(f.keys).toEqual(['x']);
    expect(f.keys).toEqual(['x']);
    expect(values).toHaveBeenCalledTimes(1);
  });
});

describe('discrete', () => {
  it('normalized over the values; the last matching bin', () => {
    const f = palette({ kind: 'discrete', colors: ['red', 'green'] }, ['10', '20']);
    expect([10, 14, 15, 20, 25].map((v) => hex(f(v)))).toEqual(['#ff0000', '#ff0000', '#008000', '#008000', undefined]);
    expect(f('abc')).toBeUndefined();
    expect(hex(f.at!({ value: 20 }))).toBe('#008000');
  });

  it('absolute ranges, open ends, no color', () => {
    const f = palette({
      kind: 'discrete',
      mode: 'absolute',
      colors: [
        ['red', null, 0],
        [null, 0, 1],
        ['blue', 1, null],
      ],
    });
    expect([-5, 0.5, 7].map((v) => hex(f(v)))).toEqual(['#ff0000', undefined, '#0000ff']);
  });

  it('one value: below, at and above it', () => {
    const f = palette({ kind: 'discrete', colors: ['red', 'green'] }, ['5']);
    expect([4, 5, 6].map((v) => hex(f(v)))).toEqual([undefined, '#008000', undefined]);
  });
});

describe('continuous', () => {
  it('interpolates between checkpoints, both halves', () => {
    const f = palette({
      kind: 'continuous',
      mode: 'absolute',
      colors: [
        ['#000000/#ffffff', 0],
        ['#ffffff/#000000', 1],
      ],
    });
    expect(hex(f(0.5))).toBe(Color.toHexStyle(Color.interpolate(Color(0), Color(0xffffff), 0.5)));
    expect(hex(f(0.25, true))).toBe(Color.toHexStyle(Color.interpolate(Color(0xffffff), Color(0), 0.25)));
    expect(hex(f(0))).toBe('#000000');
    expect(hex(f(1, true))).toBe('#000000');
    expect(hex(f.at!({ value: 0.5 }, true))).toBe(hex(f(0.5, true)));
  });

  it('under- and overflow: none, auto, a color', () => {
    const f = palette({ kind: 'continuous', mode: 'absolute', colors: ['red', 'blue'] });
    expect([f(-1), f(2)]).toEqual([undefined, undefined]);
    const g = palette({
      kind: 'continuous',
      mode: 'absolute',
      colors: ['red', 'blue'],
      underflow_color: 'auto',
      overflow_color: 'green',
    });
    expect([hex(g(-1)), hex(g(2))]).toEqual(['#ff0000', '#008000']);
  });

  it('normalized over the values and the given domain end; reversed lists', () => {
    const f = palette({ kind: 'continuous', colors: ['red', 'blue'], reverse: true, value_domain: [0, null] }, [
      '4',
      '8',
    ]);
    expect([hex(f(0)), hex(f(8))]).toEqual(['#0000ff', '#ff0000']);
    const empty = palette({ kind: 'continuous', colors: [] });
    expect(empty(0.5)).toBeUndefined();
  });

  it('normalizes over the same [min, max] as the distinct strings did (perf-audit-3 01/P14)', () => {
    const mixed: ValuePaletteValue[] = [
      -0,
      ' 5',
      '7abc',
      '1e1',
      '.',
      '?',
      null,
      undefined,
      true,
      Number.NaN,
      'nan',
      -2.5,
      '-0',
      3,
    ];
    const inputs: ValuePaletteValue[][] = [
      mixed,
      [-0, '-0'],
      ['.', '?'],
      [Number.POSITIVE_INFINITY, 1],
      ['x', null],
      [],
    ];
    for (const values of inputs) {
      for (const cifMissing of [false, true]) {
        // the former range: distinct strings, parsed, NaN dropped; 0..1 without any
        const seen = new Set<string>();
        for (const v of values) {
          if (v === null || v === undefined || (cifMissing && (v === '.' || v === '?'))) continue;
          seen.add(`${v}`);
        }
        const numbers = [...seen].map((v) => Number.parseFloat(v)).filter((x) => !Number.isNaN(x));
        const min = numbers.length ? numbers.reduce((a, b) => (a < b ? a : b)) : 0;
        const max = numbers.length ? numbers.reduce((a, b) => (a > b ? a : b)) : 1;

        const spec = (domain: [number | null, number | null]): ValuePaletteSpec => ({
          kind: 'continuous',
          colors: ['#000000', '#ffffff'],
          value_domain: domain,
          underflow_color: 'green',
          overflow_color: 'red',
        });
        const f = makeValuePaletteFunction(valuePalettePropsFromSpec(spec([null, null]), Defaults), values, {
          cifMissing,
        });
        const g = makeValuePaletteFunction(valuePalettePropsFromSpec(spec([min, max]), Defaults), [], { cifMissing });
        for (const x of [-3, -2.5, -1, -0, 0, 0.5, 3, 5, 7, 10, 11, 1e9]) expect(hex(f(x))).toBe(hex(g(x)));
      }
    }
  });

  it('does not sort the props in place', () => {
    const props = valuePalettePropsFromSpec(
      {
        kind: 'continuous',
        mode: 'absolute',
        colors: [
          ['blue', 1],
          ['red', 0],
        ],
      },
      Defaults,
    );
    const f = makeValuePaletteFunction(props, []);
    expect(hex(f(0))).toBe('#ff0000');
    expect(props.name === 'continuous' && props.params.colors.map((c) => c.checkpoint)).toEqual([1, 0]);
  });
});
