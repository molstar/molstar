/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Value palettes: colors for values (annotation fields, track values), as MolViewSpec and SeqViewSpec define them
 * (one spec: SVS's palettes are MVS's). Three steps: a palette as the spec writes it (`ValuePaletteSpec`, snake case),
 * its props (`ValuePaletteProps`, the values of `ValuePaletteParam` in `value-palette-params.ts`, what a color theme
 * saves), and the function of the props over the values it colors.
 *
 * No param definitions here, so a consumer that only evaluates palettes (Seq*'s core) loads no `ParamDefinition`.
 */

import { Color } from './color.js';
import { ColorLists } from './lists.js';
import { type OptionalSplitColorProp, SplitColor, type SplitColorProp, SplitColorValue } from './split-color.js';

// --- Props: the values of the param definitions in `value-palette-params.ts`.

export interface CategoricalValuePaletteProps {
  colors:
    | { name: 'list'; params: { color: SplitColorProp }[] }
    | { name: 'dictionary'; params: { value: string; color: SplitColorProp }[] };
  repeatColorList: boolean;
  sort: 'none' | 'lexical' | 'numeric';
  sortDirection: 'ascending' | 'descending';
  caseInsensitive: boolean;
  missingColor: OptionalSplitColorProp;
}

export interface DiscreteValuePaletteProps {
  colors: { color: OptionalSplitColorProp; fromValue: number; toValue: number }[];
  mode: 'normalized' | 'absolute';
  xMin: number | null;
  xMax: number | null;
}

export interface ContinuousValuePaletteProps {
  colors: { color: SplitColorProp; checkpoint: number }[];
  mode: 'normalized' | 'absolute';
  xMin: number | null;
  xMax: number | null;
  underflowColor: OptionalSplitColorProp;
  overflowColor: OptionalSplitColorProp;
}

/** `direct`: the values are colors themselves. */
export type ValuePaletteProps =
  | { name: 'direct'; params: {} }
  | { name: 'categorical'; params: CategoricalValuePaletteProps }
  | { name: 'discrete'; params: DiscreteValuePaletteProps }
  | { name: 'continuous'; params: ContinuousValuePaletteProps };

// --- Spec: a palette as MVS's tree (`tree/mvs/param-types.ts`) and SVS's documents write it.

/** A color as the spec writes it (`#rgb`, `#rrggbb`, a color name, a split color `a/b`), or a number-encoded color. */
export type ValuePaletteColor = string | number;

export interface CategoricalValuePaletteSpec {
  readonly kind: 'categorical';
  /** A list name, a dictionary name, colors, or `{value: Color}`. */
  readonly colors?: string | readonly ValuePaletteColor[] | Readonly<Record<string, ValuePaletteColor>>;
  readonly repeat_color_list?: boolean;
  readonly sort?: 'none' | 'lexical' | 'numeric';
  readonly sort_direction?: 'ascending' | 'descending';
  readonly case_insensitive?: boolean;
  readonly missing_color?: string | null;
}

export interface DiscreteValuePaletteSpec {
  readonly kind: 'discrete';
  /** A list name, colors, `[Color, checkpoint]` pairs, or `[Color | null, from | null, to | null]` ranges. */
  readonly colors?:
    | string
    | readonly ValuePaletteColor[]
    | readonly (readonly [ValuePaletteColor, number])[]
    | readonly (readonly [ValuePaletteColor | null, number | null, number | null])[];
  readonly reverse?: boolean;
  readonly mode?: 'normalized' | 'absolute';
  readonly value_domain?: readonly [number | null, number | null];
}

export interface ContinuousValuePaletteSpec {
  readonly kind: 'continuous';
  /** A list name, colors, or `[Color, checkpoint]` pairs. */
  readonly colors?: string | readonly ValuePaletteColor[] | readonly (readonly [ValuePaletteColor, number])[];
  readonly reverse?: boolean;
  readonly mode?: 'normalized' | 'absolute';
  readonly value_domain?: readonly [number | null, number | null];
  readonly underflow_color?: 'auto' | string | null;
  readonly overflow_color?: 'auto' | string | null;
}

export type ValuePaletteSpec = CategoricalValuePaletteSpec | DiscreteValuePaletteSpec | ContinuousValuePaletteSpec;

/** Every field of every kind: MVS's defaults (`CategoricalPaletteDefaults`, ...) or SVS's (`PaletteDefaults`). */
export interface ValuePaletteSpecDefaults {
  readonly categorical: Required<CategoricalValuePaletteSpec>;
  readonly discrete: Required<DiscreteValuePaletteSpec>;
  readonly continuous: Required<ContinuousValuePaletteSpec>;
}

/** MVS's color list names (`ColorListNameT`) and the lists of `ColorLists` they name, in MVS's order. */
const ListNames = {
  // Sequential single-hue
  Reds: 'reds',
  Oranges: 'oranges',
  Greens: 'greens',
  Blues: 'blues',
  Purples: 'purples',
  Greys: 'greys',

  // Sequential multi-hue
  OrRd: 'orange-red',
  BuGn: 'blue-green',
  PuBuGn: 'purple-blue-green',
  GnBu: 'green-blue',
  PuBu: 'purple-blue',
  BuPu: 'blue-purple',
  RdPu: 'red-purple',
  PuRd: 'purple-red',
  YlOrRd: 'yellow-orange-red',
  YlOrBr: 'yellow-orange-brown',
  YlGn: 'yellow-green',
  YlGnBu: 'yellow-green-blue',

  Magma: 'magma',
  Inferno: 'inferno',
  Plasma: 'plasma',
  Viridis: 'viridis',
  Cividis: 'cividis',
  Turbo: 'turbo',
  Warm: 'warm',
  Cool: 'cool',
  CubehelixDefault: 'cubehelix-default',

  // Cyclical
  Rainbow: 'rainbow',
  Sinebow: 'sinebow',

  // Diverging
  RdBu: 'red-blue',
  RdGy: 'red-grey',
  PiYG: 'pink-yellow-green',
  BrBG: 'brown-white-green',
  PRGn: 'purple-green',
  PuOr: 'purple-orange',
  RdYlGn: 'red-yellow-green',
  RdYlBu: 'red-yellow-blue',
  Spectral: 'spectral',

  // Categorical
  Category10: 'category-10',
  Observable10: 'observable-10',
  Tableau10: 'tableau-10',

  Set1: 'set-1',
  Set2: 'set-2',
  Set3: 'set-3',
  Pastel1: 'pastel-1',
  Pastel2: 'pastel-2',
  Dark2: 'dark-2',
  Paired: 'paired',
  Accent: 'accent',

  // Additional lists, not standard for visualization in general, but commonly used for structures
  Chainbow: 'turbo-no-black',
} as const satisfies Record<string, keyof typeof ColorLists>;

export type ValuePaletteListName = keyof typeof ListNames;

/** The named color lists a palette's `colors` may name (MVS's `ColorListNameT`), as colors. */
export const ValuePaletteLists = Object.fromEntries(
  Object.entries(ListNames).map(([name, list]) => [name, ColorLists[list].list.map(Color.fromColorListEntry)]),
) as unknown as Readonly<Record<ValuePaletteListName, readonly Color[]>>;

/** Named color dictionaries (MVS's `ColorDictNameT`): their tables live with their users (MVS's from Mol*'s themes). */
export type ValuePaletteDictionaries = Readonly<Record<string, Readonly<Record<string, ValuePaletteColor>>>>;

export interface ValuePalettePropsOptions {
  /** The dictionaries a categorical palette's `colors` may name; lists are `ValuePaletteLists`. */
  readonly dictionaries?: ValuePaletteDictionaries;
  /** Called with a `colors` name that is neither a list nor (for a categorical palette) a dictionary. */
  readonly onUnknownName?: (name: string) => void;
}

/** A palette's props: `palette` over `defaults` (`undefined` keeps a default, `null` is a value); none is `direct`. */
export function valuePalettePropsFromSpec(
  palette: ValuePaletteSpec | undefined | null,
  defaults: ValuePaletteSpecDefaults,
  options: ValuePalettePropsOptions = {},
): ValuePaletteProps {
  if (!palette) {
    return { name: 'direct', params: {} };
  }
  if (palette.kind === 'categorical') {
    const full = merge(defaults.categorical, palette);
    return {
      name: 'categorical',
      params: {
        colors: categoricalColors(full.colors, options),
        repeatColorList: full.repeat_color_list,
        sort: full.sort,
        sortDirection: full.sort_direction,
        caseInsensitive: full.case_insensitive,
        missingColor: SplitColorValue.fromString(full.missing_color),
      },
    };
  }
  if (palette.kind === 'discrete') {
    const full = merge(defaults.discrete, palette);
    return {
      name: 'discrete',
      params: {
        colors: discreteColors(full.colors, full.reverse, options),
        mode: full.mode,
        xMin: full.value_domain[0],
        xMax: full.value_domain[1],
      },
    };
  }
  if (palette.kind === 'continuous') {
    const full = merge(defaults.continuous, palette);
    const colors = continuousColors(full.colors, full.reverse, options);
    return {
      name: 'continuous',
      params: {
        colors,
        mode: full.mode,
        xMin: full.value_domain[0],
        xMax: full.value_domain[1],
        underflowColor:
          full.underflow_color === 'auto'
            ? extremeColor(colors, (a, b) => a < b)
            : SplitColorValue.fromString(full.underflow_color),
        overflowColor:
          full.overflow_color === 'auto'
            ? extremeColor(colors, (a, b) => a > b)
            : SplitColorValue.fromString(full.overflow_color),
      },
    };
  }
  throw new Error(
    `NotImplementedError: valuePalettePropsFromSpec is not implemented for palette kind "${(palette as any).kind}"`,
  );
}

/** Properties of `second` over `first`'s: `undefined` is a missing property, `null` a value. */
function merge<T extends object, U extends object>(first: T, second: U): T & U {
  const out: Partial<T & U> = { ...first };
  for (const key in second) {
    const value = second[key];
    if (value !== undefined) out[key] = value as any;
  }
  return out as T & U;
}

function splitColor(color: ValuePaletteColor): SplitColorProp {
  return typeof color === 'number'
    ? { name: 'oneColor', params: { color: color as Color } }
    : SplitColorValue.fromString(color as string);
}

function optionalSplitColor(color: ValuePaletteColor | null): OptionalSplitColorProp {
  return color === null ? SplitColorValue.none() : splitColor(color);
}

function namedList(name: string): readonly Color[] | undefined {
  return Object.hasOwn(ValuePaletteLists, name) ? ValuePaletteLists[name as ValuePaletteListName] : undefined;
}

function categoricalColors(
  colors: Required<CategoricalValuePaletteSpec>['colors'],
  options: ValuePalettePropsOptions,
): CategoricalValuePaletteProps['colors'] {
  if (typeof colors === 'string') {
    // Named color list or dict
    const name = colors;
    const resolved =
      namedList(name) ??
      (options.dictionaries && Object.hasOwn(options.dictionaries, name) ? options.dictionaries[name] : undefined);
    if (!resolved) {
      options.onUnknownName?.(name);
      return { name: 'list', params: [] };
    }
    colors = resolved;
  }

  if (Array.isArray(colors)) {
    // Color list
    return { name: 'list', params: (colors as readonly ValuePaletteColor[]).map((c) => ({ color: splitColor(c) })) };
  } else {
    // Color dict
    return {
      name: 'dictionary',
      params: Object.entries(colors as Readonly<Record<string, ValuePaletteColor>>).map(([value, color]) => ({
        value,
        color: splitColor(color),
      })),
    };
  }
}

function discreteColors(
  colors: Required<DiscreteValuePaletteSpec>['colors'],
  reverse: boolean,
  options: ValuePalettePropsOptions,
): DiscreteValuePaletteProps['colors'] {
  if (typeof colors === 'string') {
    // Named color list
    const list = namedList(colors);
    if (!list) {
      options.onUnknownName?.(colors);
      return [];
    }
    colors = list;
  }

  if (colors.every((t) => !Array.isArray(t))) {
    // Color list without checkpoints
    let list = colors as readonly ValuePaletteColor[];
    if (reverse) list = list.slice().reverse();
    const sectionLength = 1 / list.length;
    return list.map((col, i) => ({
      color: splitColor(col),
      fromValue: i * sectionLength,
      toValue: (i + 1) * sectionLength,
    }));
  } else if (colors.every((t) => Array.isArray(t) && t.length === 2)) {
    // Color list with 1 checkpoint
    const pairs = colors as readonly (readonly [ValuePaletteColor, number])[];
    return pairs.map(([col, from], i) => ({
      color: splitColor(col),
      fromValue: from,
      toValue: pairs[i + 1]?.[1] ?? Infinity,
    }));
  } else if (colors.every((t) => Array.isArray(t) && t.length === 3)) {
    // Color list with 2 checkpoints
    const ranges = colors as readonly (readonly [ValuePaletteColor | null, number | null, number | null])[];
    return ranges.map(([col, from, to]) => ({
      color: optionalSplitColor(col),
      fromValue: from ?? -Infinity,
      toValue: to ?? Infinity,
    }));
  } else {
    return []; // should be unreachable
  }
}

function continuousColors(
  colors: Required<ContinuousValuePaletteSpec>['colors'],
  reverse: boolean,
  options: ValuePalettePropsOptions,
): ContinuousValuePaletteProps['colors'] {
  if (typeof colors === 'string') {
    // Named color list
    const list = namedList(colors);
    if (!list) {
      options.onUnknownName?.(colors);
      return [];
    }
    colors = list;
  }

  if (colors.every((t) => !Array.isArray(t))) {
    // Color list without checkpoints
    let list = colors as readonly ValuePaletteColor[];
    if (reverse) list = list.slice().reverse();
    const n = list.length - 1;
    return list.map((col, i) => ({ color: splitColor(col), checkpoint: i / n }));
  } else {
    // Color list with checkpoints
    // Not applying `reverse` here, as it would have no effect
    return (colors as readonly (readonly [ValuePaletteColor, number])[]).map(([col, checkpoint]) => ({
      color: splitColor(col),
      checkpoint,
    }));
  }
}

/** The color of the lowest (`a < b`) or the highest (`a > b`) checkpoint; on a tie the later one. */
function extremeColor(
  colors: ContinuousValuePaletteProps['colors'],
  pick: (a: number, b: number) => boolean,
): OptionalSplitColorProp {
  if (colors.length === 0) return SplitColorValue.none();
  return colors.reduce((a, b) => (pick(a.checkpoint, b.checkpoint) ? a : b)).color;
}

// --- Evaluation

/** A value as a palette reads it: compared as a string, or parsed as a number. */
export type ValuePaletteValue = string | number | boolean | null | undefined;

/**
 * The color of one value (a split color's second half when `isSecondary`); `undefined`: not colored by this palette.
 */
export type ValuePaletteFunction = ((value: ValuePaletteValue, isSecondary?: boolean) => Color | undefined) & {
  /**
   * A categorical palette's values that it colors, in its order (its `sort`, case-folded under `caseInsensitive`): what
   * a legend lists.
   */
  readonly keys?: readonly string[];
  /**
   * A continuous or discrete palette at `o.value`, the same as the call: per element on a frame no double crosses a
   * call (V8 boxes a double passed to one it does not inline), so a draw's `palette.at(o)` allocates nothing.
   */
  readonly at?: (o: { readonly value: number }, isSecondary?: boolean) => Color | undefined;
};

export interface ValuePaletteOptions {
  /**
   * `.` and `?` (CIF's inapplicable and unknown) are missing values, as `null` and `undefined` always are: they never
   * take a list color, add nothing to a normalized domain, and a categorical palette colors them by its `missingColor`
   * alone. SVS sets it (seqspec:format#palettes); MVS reads them as values (its CIF annotations skip them before).
   */
  readonly cifMissing?: boolean;
}

/**
 * The palette's function. `values`: every value it will color (a categorical list's order of first occurrence and a
 * normalized domain come from them), or a function returning them, called at most once and only when needed (a list,
 * a normalized domain without both ends, a dictionary's `keys` when read: then a getter).
 */
export function makeValuePaletteFunction(
  props: ValuePaletteProps,
  values: Iterable<ValuePaletteValue> | (() => Iterable<ValuePaletteValue>),
  options: ValuePaletteOptions = {},
): ValuePaletteFunction {
  let resolved: Iterable<ValuePaletteValue> | undefined;
  const getValues = () => (resolved ??= typeof values === 'function' ? values() : values);
  const cifMissing = !!options.cifMissing;
  if (props.name === 'direct') return paletteFunctionDirect;
  if (props.name === 'categorical')
    return categorical(props.params, getValues, cifMissing, typeof values === 'function');
  if (props.name === 'discrete') return discrete(props.params, getValues, cifMissing);
  if (props.name === 'continuous') return continuous(props.params, getValues, cifMissing);
  throw new Error(`NotImplementedError: makeValuePaletteFunction for ${(props as any).name}`);
}

type Values = () => Iterable<ValuePaletteValue>;

function isMissing(value: ValuePaletteValue, cifMissing: boolean): boolean {
  return value === null || value === undefined || (cifMissing && (value === '.' || value === '?'));
}

/** Distinct values in order of first occurrence, as strings (MVS's `getDistinctValuesInField`). */
function distinct(values: Iterable<ValuePaletteValue>, upper: boolean, cifMissing: boolean): string[] {
  const seen = new Set<string>();
  for (const value of values) {
    if (isMissing(value, cifMissing)) continue;
    seen.add(upper ? `${value}`.toUpperCase() : `${value}`);
  }
  return [...seen];
}

/** A value as a number (`parseFloat` of its string). */
function toNumber(value: ValuePaletteValue, cifMissing: boolean): number {
  if (typeof value === 'number') return value;
  return isMissing(value, cifMissing) ? Number.NaN : Number.parseFloat(`${value}`);
}

const _directColors: [Color | undefined, Color | undefined] = [undefined, undefined];

const paletteFunctionDirect: ValuePaletteFunction = (value, isSecondary) => {
  if (value === null || value === undefined) return undefined;
  SplitColor.decodeTo(`${value}`, _directColors);
  return _directColors[isSecondary ? 1 : 0];
};

function categorical(
  p: CategoricalValuePaletteProps,
  values: Values,
  cifMissing: boolean,
  lazy: boolean,
): ValuePaletteFunction {
  const map = new Map<string, readonly [Color | undefined, Color | undefined]>();
  const upper = p.caseInsensitive;
  let sortedKeys: string[] | undefined;
  const getSortedKeys = () => {
    if (sortedKeys) return sortedKeys;
    const keys = distinct(values(), upper, cifMissing);
    if (p.sort === 'lexical') keys.sort();
    else if (p.sort === 'numeric') keys.sort((a, b) => Number.parseFloat(a) - Number.parseFloat(b));
    if (p.sortDirection === 'descending') keys.reverse();
    return (sortedKeys = keys);
  };
  if (p.colors.name === 'dictionary') {
    for (const { value, color } of p.colors.params) {
      map.set(upper ? value.toUpperCase() : value, SplitColorValue.toTuple(color));
    }
  } else {
    const list = p.colors.params.map((item) => SplitColorValue.toTuple(item.color));
    let next = 0;
    for (const key of getSortedKeys()) {
      if (next < list.length) map.set(key, list[next]);
      next++;
      if (next >= list.length) {
        if (p.repeatColorList) next = 0;
        else break;
      }
    }
  }
  const missing = SplitColorValue.toTuple(p.missingColor);
  const f = (value: ValuePaletteValue, isSecondary?: boolean) => {
    // A missing value takes `missingColor` too: a kernel fallback draws "no category".
    if (isMissing(value, cifMissing)) return missing[isSecondary ? 1 : 0];
    const key = upper ? `${value}`.toUpperCase() : `${value}`;
    return (map.get(key) ?? missing)[isSecondary ? 1 : 0];
  };
  const listed = () => getSortedKeys().filter((k) => map.has(k));
  if (!lazy) return Object.assign(f, { keys: listed() });
  // Lazy values: a dictionary reads them for its `keys` alone, so only when they are read.
  let keys: string[] | undefined;
  return Object.defineProperty(f, 'keys', { get: () => (keys ??= listed()), enumerable: true }) as ValuePaletteFunction;
}

function discrete(p: DiscreteValuePaletteProps, values: Values, cifMissing: boolean): ValuePaletteFunction {
  const bins = p.colors;
  if (bins.length === 0) return () => undefined;
  // As the continuous one: the scale inline, the bins' edges in typed arrays and the value read from `held`, so `at`
  // allocates nothing on a frame (a draw's `palette.at(o)` per position or column).
  const from = Float64Array.from(bins, (x) => x.fromValue);
  const to = Float64Array.from(bins, (x) => x.toValue);
  const colors1 = bins.map((x) => SplitColorValue.toTuple(x.color)[0]);
  const colors2 = bins.map((x) => SplitColorValue.toTuple(x.color)[1]);
  const range = numericRange(p, values, cifMissing);
  const held = new Float64Array(1);
  const core = (isSecondary: boolean): Color | undefined => {
    const xAbs = held[0];
    if (Number.isNaN(xAbs)) return undefined;
    const lo = range[0];
    const hi = range[1];
    const x = Number.isNaN(lo)
      ? xAbs
      : lo === hi
        ? xAbs < lo
          ? -0.5
          : xAbs === lo
            ? 0.5
            : 1.5
        : (xAbs - lo) / (hi - lo);
    for (let i = from.length - 1; i >= 0; i--) {
      if (from[i] <= x && x <= to[i]) return isSecondary ? colors2[i] : colors1[i];
    }
    return undefined;
  };
  const at = (o: { readonly value: number }, isSecondary?: boolean) => {
    held[0] = o.value;
    return core(isSecondary === true);
  };
  return Object.assign(
    (value: ValuePaletteValue, isSecondary?: boolean) => {
      held[0] = toNumber(value, cifMissing);
      return core(isSecondary === true);
    },
    { at },
  );
}

function continuous(p: ContinuousValuePaletteProps, values: Values, cifMissing: boolean): ValuePaletteFunction {
  const stops = p.colors.slice().sort((a, b) => a.checkpoint - b.checkpoint);
  if (stops.length === 0) return () => undefined;
  const cps = Float64Array.from(stops, (s) => s.checkpoint);
  const colors1 = Int32Array.from(stops, (s) => SplitColorValue.toTuple(s.color)[0]);
  const colors2 = Int32Array.from(stops, (s) => SplitColorValue.toTuple(s.color)[1]);
  const underflow = SplitColorValue.toTuple(p.underflowColor);
  const overflow = SplitColorValue.toTuple(p.overflowColor);
  const range = numericRange(p, values, cifMissing);
  // Inline throughout: the scale, the first checkpoint at or past x (`SortedArray.findPredecessorIndex`),
  // `Color.interpolate`. The value is read from `held`: `at`'s loads see only its callers' locations (a load
  // polymorphic over a double field boxes it).
  const held = new Float64Array(1);
  const core = (isSecondary: boolean): Color | undefined => {
    const xAbs = held[0];
    if (Number.isNaN(xAbs)) return undefined;
    const lo = range[0];
    const hi = range[1];
    const x = Number.isNaN(lo)
      ? xAbs
      : lo === hi
        ? xAbs < lo
          ? -0.5
          : xAbs === lo
            ? 0.5
            : 1.5
        : (xAbs - lo) / (hi - lo);
    const colorsAt = isSecondary ? colors2 : colors1;
    const n = cps.length;
    let i = 0;
    if (cps[0] < x) {
      i = n;
      for (let k = 1; k < n; k++)
        if (cps[k] >= x) {
          i = k;
          break;
        }
    }
    if (i === 0) return x === cps[0] ? (colorsAt[0] as Color) : underflow[isSecondary ? 1 : 0];
    if (i === n) return overflow[isSecondary ? 1 : 0];
    const t = (x - cps[i - 1]) / (cps[i] - cps[i - 1]);
    const c1 = colorsAt[i - 1];
    const c2 = colorsAt[i];
    const r1 = (c1 >> 16) & 255;
    const g1 = (c1 >> 8) & 255;
    const b1 = c1 & 255;
    const r = r1 + (((c2 >> 16) & 255) - r1) * t;
    const g = g1 + (((c2 >> 8) & 255) - g1) * t;
    const b = b1 + ((c2 & 255) - b1) * t;
    return ((r << 16) | (g << 8) | b) as Color;
  };
  const at = (o: { readonly value: number }, isSecondary?: boolean) => {
    held[0] = o.value;
    return core(isSecondary === true);
  };
  return Object.assign(
    (value: ValuePaletteValue, isSecondary?: boolean) => {
      held[0] = toNumber(value, cifMissing);
      return core(isSecondary === true);
    },
    { at },
  );
}

/**
 * A `normalized` palette's [min, max] over `xMin`, `xMax`, their `null`s from the values (0..1 without any), the value
 * scaled `(x - min) / (max - min)`; [NaN, NaN]: absolute.
 */
function numericRange(
  p: DiscreteValuePaletteProps | ContinuousValuePaletteProps,
  values: Values,
  cifMissing: boolean,
): Float64Array {
  if (p.mode !== 'normalized') return Float64Array.of(Number.NaN, Number.NaN);
  let xMin = p.xMin;
  let xMax = p.xMax;
  if (xMin === null || xMax === null) {
    // one pass over the values, no distinct strings (91 MB per 1.1M values, perf-audit-3 01/P14): `parseFloat` of a
    // number's string is the number, except −0, which its string makes +0
    let min = Number.POSITIVE_INFINITY;
    let max = Number.NEGATIVE_INFINITY;
    let any = false;
    for (const value of values()) {
      const x = toNumber(value, cifMissing);
      if (Number.isNaN(x)) continue;
      const v = x === 0 ? 0 : x;
      if (v < min) min = v;
      if (v > max) max = v;
      any = true;
    }
    if (any) {
      xMin ??= min;
      xMax ??= max;
    } else {
      xMin ??= 0;
      xMax ??= 1;
    }
  }
  return Float64Array.of(xMin, xMax);
}
