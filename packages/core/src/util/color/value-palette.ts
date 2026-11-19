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

import { hashFnv32a } from '../../data/util/hash-functions.js';
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
  /**
   * A value the palette does not list: `missing` (its `missingColor`, the default) or `hash`, the color at the 32-bit
   * FNV-1a hash of the value's UTF-8 bytes modulo `others` (a list palette: its own list), so it takes one color
   * everywhere and another reader can reproduce it.
   */
  unlisted?: 'missing' | 'hash';
  /** The colors an unlisted value hashes into (a dictionary's); a list palette hashes into its list when absent. */
  others?: { color: SplitColorProp }[];
}

/** How a numeric palette reads values: as they are, by their logarithm (base 10), or by their square root. */
export type ValuePaletteScale = 'linear' | 'log' | 'sqrt';

/** Optional fields of the numeric palettes (absent: a linear scale, no center). */
export interface ValuePaletteScaleProps {
  /** `linear` when absent. */
  scale?: ValuePaletteScale;
  /** On a log scale, a value at or below zero: `missing` (when absent) or `clamp`ed to the domain's low end. */
  nonPositive?: 'missing' | 'clamp';
  /**
   * A diverging center (a value, read through the scale): a normalized palette puts it at 0.5, the domain's low end at
   * 0 and its high end at 1, each half scaled on its own.
   */
  mid?: number | null;
}

export interface DiscreteValuePaletteProps extends ValuePaletteScaleProps {
  colors: { color: OptionalSplitColorProp; fromValue: number; toValue: number }[];
  mode: 'normalized' | 'absolute';
  xMin: number | null;
  xMax: number | null;
}

export interface ContinuousValuePaletteProps extends ValuePaletteScaleProps {
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
  /** A value not listed: `missing` (`missing_color`, the default) or `hash` (FNV-1a over UTF-8 into `others`). */
  readonly unlisted?: 'missing' | 'hash';
  /** A list name or colors that an unlisted value hashes into; a list palette's own list when absent. */
  readonly others?: string | readonly ValuePaletteColor[];
}

/** The optional scale fields of the numeric palettes' spec. */
export interface ValuePaletteScaleSpec {
  readonly scale?: ValuePaletteScale;
  readonly non_positive?: 'missing' | 'clamp';
  readonly mid?: number | null;
}

export interface DiscreteValuePaletteSpec extends ValuePaletteScaleSpec {
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

export interface ContinuousValuePaletteSpec extends ValuePaletteScaleSpec {
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
  readonly categorical: Required<Omit<CategoricalValuePaletteSpec, 'unlisted' | 'others'>> &
    Pick<CategoricalValuePaletteSpec, 'unlisted' | 'others'>;
  readonly discrete: Required<Omit<DiscreteValuePaletteSpec, keyof ValuePaletteScaleSpec>> & ValuePaletteScaleSpec;
  readonly continuous: Required<Omit<ContinuousValuePaletteSpec, keyof ValuePaletteScaleSpec>> & ValuePaletteScaleSpec;
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

const _lists: Partial<Record<ValuePaletteListName, readonly Color[]>> = {};

/**
 * The named color lists a palette's `colors` may name (MVS's `ColorListNameT`), as colors. Each list is converted when
 * first read (a getter per name), so loading the module converts none.
 */
export const ValuePaletteLists = {} as Readonly<Record<ValuePaletteListName, readonly Color[]>>;
for (const name of Object.keys(ListNames) as ValuePaletteListName[]) {
  Object.defineProperty(ValuePaletteLists, name, {
    enumerable: true,
    get: () => (_lists[name] ??= ColorLists[ListNames[name]].list.map(Color.fromColorListEntry)),
  });
}

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
        ...(full.unlisted ? { unlisted: full.unlisted } : {}),
        ...(full.others !== undefined ? { others: otherColors(full.others, options) } : {}),
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
        ...scaleProps(full),
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
        ...scaleProps(full),
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

function scaleProps(full: ValuePaletteScaleSpec): ValuePaletteScaleProps {
  return {
    ...(full.scale && full.scale !== 'linear' ? { scale: full.scale } : {}),
    ...(full.non_positive ? { nonPositive: full.non_positive } : {}),
    ...(typeof full.mid === 'number' ? { mid: full.mid } : {}),
  };
}

function otherColors(
  others: string | readonly ValuePaletteColor[],
  options: ValuePalettePropsOptions,
): { color: SplitColorProp }[] {
  let list: readonly ValuePaletteColor[] | undefined;
  if (typeof others === 'string') {
    list = namedList(others);
    if (!list) options.onUnknownName?.(others);
  } else list = others;
  return (list ?? []).map((c) => ({ color: splitColor(c) }));
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
  const unlisted = p.unlisted === 'hash' ? hashedColors(p) : undefined;
  const f = (value: ValuePaletteValue, isSecondary?: boolean) => {
    // A missing value takes `missingColor` too: a kernel fallback draws "no category".
    if (isMissing(value, cifMissing)) return missing[isSecondary ? 1 : 0];
    const key = upper ? `${value}`.toUpperCase() : `${value}`;
    return (map.get(key) ?? unlisted?.(key) ?? missing)[isSecondary ? 1 : 0];
  };
  const listed = () => getSortedKeys().filter((k) => map.has(k));
  if (!lazy) return Object.assign(f, { keys: listed() });
  // Lazy values: a dictionary reads them for its `keys` alone, so only when they are read.
  let keys: string[] | undefined;
  return Object.defineProperty(f, 'keys', { get: () => (keys ??= listed()), enumerable: true }) as ValuePaletteFunction;
}

const _utf8 = new TextEncoder();

/** The 32-bit FNV-1a hash of a string's UTF-8 bytes: how an unlisted value picks its color. */
export function hashValueUtf8(value: string): number {
  return hashFnv32a(_utf8.encode(value));
}

/** An unlisted value's color: `others` (else the list) at its hash; memoized per value, at most 4096 values. */
function hashedColors(p: CategoricalValuePaletteProps) {
  const source = p.others ?? (p.colors.name === 'list' ? p.colors.params : []);
  const list = source.map((x) => SplitColorValue.toTuple(x.color));
  if (list.length === 0) return undefined;
  const memo = new Map<string, readonly [Color | undefined, Color | undefined]>();
  return (key: string) => {
    let c = memo.get(key);
    if (c) return c;
    c = list[hashValueUtf8(key) % list.length];
    if (memo.size < 4096) memo.set(key, c);
    return c;
  };
}

/** The scale's transform kind: 0 linear, 1 log, 2 sqrt (a number, so the hot path branches without a call). */
function scaleKind(p: ValuePaletteScaleProps): number {
  return p.scale === 'log' ? 1 : p.scale === 'sqrt' ? 2 : 0;
}

/** A value through the scale (log10, sqrt; NaN outside its domain). At setup only, never per value. */
function scaled(kind: number, x: number): number {
  if (kind === 1) return x > 0 ? Math.log10(x) : Number.NaN;
  if (kind === 2) return x >= 0 ? Math.sqrt(x) : Number.NaN;
  return x;
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
  const kind = scaleKind(p);
  const clamp = p.nonPositive === 'clamp';
  const mid = typeof p.mid === 'number' ? scaled(kind, p.mid) : Number.NaN;
  const hasMid = !Number.isNaN(mid);
  if (kind !== 0 && Number.isNaN(range[0])) {
    // absolute edges are values: read through the scale as the values are
    for (let i = 0; i < from.length; i++) {
      if (Number.isFinite(from[i])) from[i] = scaled(kind, from[i]);
      if (Number.isFinite(to[i])) to[i] = scaled(kind, to[i]);
    }
  }
  const core = (isSecondary: boolean): Color | undefined => {
    let xAbs = held[0];
    if (kind === 1) xAbs = xAbs > 0 ? Math.log10(xAbs) : clamp ? range[2] : Number.NaN;
    else if (kind === 2) xAbs = xAbs >= 0 ? Math.sqrt(xAbs) : Number.NaN;
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
        : hasMid
          ? xAbs < mid
            ? mid === lo
              ? -0.5
              : (0.5 * (xAbs - lo)) / (mid - lo)
            : hi === mid
              ? xAbs === mid
                ? 0.5
                : 1.5
              : 0.5 + (0.5 * (xAbs - mid)) / (hi - mid)
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
  const kind = scaleKind(p);
  const clamp = p.nonPositive === 'clamp';
  const mid = typeof p.mid === 'number' ? scaled(kind, p.mid) : Number.NaN;
  const hasMid = !Number.isNaN(mid);
  // absolute checkpoints are values: read through the scale as the values are
  if (kind !== 0 && Number.isNaN(range[0])) for (let i = 0; i < cps.length; i++) cps[i] = scaled(kind, cps[i]);
  const core = (isSecondary: boolean): Color | undefined => {
    let xAbs = held[0];
    if (kind === 1) xAbs = xAbs > 0 ? Math.log10(xAbs) : clamp ? range[2] : Number.NaN;
    else if (kind === 2) xAbs = xAbs >= 0 ? Math.sqrt(xAbs) : Number.NaN;
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
        : hasMid
          ? xAbs < mid
            ? mid === lo
              ? -0.5
              : (0.5 * (xAbs - lo)) / (mid - lo)
            : hi === mid
              ? xAbs === mid
                ? 0.5
                : 1.5
              : 0.5 + (0.5 * (xAbs - mid)) / (hi - mid)
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
 * scaled `(x - min) / (max - min)`; [NaN, NaN]: absolute. Through the palette's scale (log, sqrt): the ends and the
 * values are read through it, a value outside its domain skipped. The third number is what a clamped non-positive
 * value reads as: the low end (absolute: −Infinity, under the first checkpoint).
 */
function numericRange(
  p: DiscreteValuePaletteProps | ContinuousValuePaletteProps,
  values: Values,
  cifMissing: boolean,
): Float64Array {
  if (p.mode !== 'normalized') return Float64Array.of(Number.NaN, Number.NaN, Number.NEGATIVE_INFINITY);
  const kind = scaleKind(p);
  let xMin = p.xMin === null ? null : scaled(kind, p.xMin);
  let xMax = p.xMax === null ? null : scaled(kind, p.xMax);
  if (Number.isNaN(xMin)) xMin = null;
  if (Number.isNaN(xMax)) xMax = null;
  if (xMin === null || xMax === null) {
    // one pass over the values, no distinct strings (91 MB per 1.1M values, perf-audit-3 01/P14): `parseFloat` of a
    // number's string is the number, except −0, which its string makes +0
    let min = Number.POSITIVE_INFINITY;
    let max = Number.NEGATIVE_INFINITY;
    let any = false;
    for (const value of values()) {
      const x = scaled(kind, toNumber(value, cifMissing));
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
  return Float64Array.of(xMin, xMax, xMin);
}
