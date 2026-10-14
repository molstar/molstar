/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * MVS's palettes: the same kinds, params, defaults and evaluation, so a value gets one color in SVS and in MVS
 * (seqspec:format#themes).
 *
 * Copied from Mol* (MIT, (c) 2023-2026 mol* contributors; Mol*'s src/extensions/mvs; author Adam Midlik):
 * the palette types and defaults of `tree/mvs/param-types.ts`; the conversion of `load-helpers.ts`
 * (`palettePropsFromMVSPalette` and the functions it calls); the palette functions of
 * `components/annotation-color-theme.ts` (`makePaletteFunction*`, `makeNumericPaletteScale`,
 * `makeContinuousPaletteCheckpoints`); `SortedArray.findPredecessorIndex` from `mol-data/int` (David Sehnal).
 * Mol*'s imports removed, the two steps (props, then function) joined. Differences, all from the spec: missing
 * values (`null`, `undefined`) and the strings `.` and `?` never take a list color (MVS skips them only in CIF
 * annotations) and are colored by a categorical palette's `missing_color` alone (MVS leaves them uncolored); colors are
 * the first half of a split color; `repeat_color_list` defaults to `true` (a list starts over rather than leaving the
 * values past its end uncolored; MVS's default is `false`).
 */

import { Color } from "./color.ts";
import { ColorDictionaries, type ColorDictionaryName } from "./dictionaries.ts";
import { type ColorListName, ColorLists } from "./lists.ts";

/** A document `Color` string (`#rgb`, `#rrggbb`, a Mol* color name, `a/b`). */
export type ColorString = string;

export interface CategoricalPalette {
  readonly kind: "categorical";
  /** A list name, a dictionary name, colors, or `{value: Color}`. */
  readonly colors?:
    | ColorListName
    | ColorDictionaryName
    | readonly ColorString[]
    | Readonly<Record<string, ColorString>>;
  readonly repeat_color_list?: boolean;
  readonly sort?: "none" | "lexical" | "numeric";
  readonly sort_direction?: "ascending" | "descending";
  readonly case_insensitive?: boolean;
  readonly missing_color?: ColorString | null;
}

export interface DiscretePalette {
  readonly kind: "discrete";
  /** A list name, colors, `[Color, checkpoint]` pairs, or `[Color | null, from | null, to | null]` ranges. */
  readonly colors?:
    | ColorListName
    | readonly ColorString[]
    | readonly (readonly [ColorString, number])[]
    | readonly (readonly [ColorString | null, number | null, number | null])[];
  readonly reverse?: boolean;
  readonly mode?: "normalized" | "absolute";
  readonly value_domain?: readonly [number | null, number | null];
}

export interface ContinuousPalette {
  readonly kind: "continuous";
  /** A list name, colors, or `[Color, checkpoint]` pairs. */
  readonly colors?: ColorListName | readonly ColorString[] | readonly (readonly [ColorString, number])[];
  readonly reverse?: boolean;
  readonly mode?: "normalized" | "absolute";
  readonly value_domain?: readonly [number | null, number | null];
  readonly underflow_color?: "auto" | ColorString | null;
  readonly overflow_color?: "auto" | ColorString | null;
}

export type Palette = CategoricalPalette | DiscretePalette | ContinuousPalette;

/** MVS's defaults (`CategoricalPaletteDefaults`, `DiscretePaletteDefaults`, `ContinuousPaletteDefaults`), but a list repeats. */
export const PaletteDefaults = {
  categorical: {
    kind: "categorical",
    colors: "Category10",
    repeat_color_list: true,
    sort: "none",
    sort_direction: "ascending",
    case_insensitive: false,
    missing_color: null,
  } satisfies Required<CategoricalPalette>,
  discrete: {
    kind: "discrete",
    colors: "YlGn",
    reverse: false,
    mode: "normalized",
    value_domain: [null, null],
  } satisfies Required<DiscretePalette>,
  continuous: {
    kind: "continuous",
    colors: "YlGn",
    reverse: false,
    mode: "normalized",
    value_domain: [null, null],
    underflow_color: null,
    overflow_color: null,
  } satisfies Required<ContinuousPalette>,
} as const;

/** A field's value as a palette reads it: compared as a string, or parsed as a number. */
export type PaletteValue = string | number | boolean | null | undefined;
/**
 * The color of one value; `undefined`: not colored by this layer. A categorical palette's function carries `keys`: the
 * values it colors, in its order (its `sort`, case-folded under `case_insensitive`), what a legend lists.
 */
export type PaletteFunction = ((value: PaletteValue) => Color | undefined) & {
  readonly keys?: readonly string[];
  /**
   * A continuous or discrete palette at `o.value`, the same as the call: per element on a frame no double crosses a call (V8
   * boxes a double passed to one it does not inline), so a draw's `palette(mean)` allocates nothing.
   */
  readonly at?: (o: { readonly value: number }) => Color | undefined;
  /** A normalized palette's `value_domain` as `layerPalette` filled it: what a legend's ramp spans. */
  readonly domain?: readonly [number | null, number | null];
  /**
   * `layerPalette` over streamed data with no vocabulary of the field: its values are colored by their hash
   * (`paletteHash`), so `keys` cannot list them.
   */
  readonly unlisted?: boolean;
};

/** Missing values, `.` and `?`: no list color and no domain; a categorical palette's `missing_color` alone. */
function isMissing(value: PaletteValue): boolean {
  return value === null || value === undefined || value === "." || value === "?";
}

/**
 * The palette's function over `values`: every value of the field over the layer's data, annotations in `data`
 * order, rows in order (first occurrence and the normalized domain come from them).
 */
export function makePaletteFunction(palette: Palette, values: Iterable<PaletteValue>): PaletteFunction {
  if (palette.kind === "categorical")
    return categorical({ ...PaletteDefaults.categorical, ...defined(palette) }, values);
  if (palette.kind === "discrete") return discrete({ ...PaletteDefaults.discrete, ...defined(palette) }, values);
  return continuous({ ...PaletteDefaults.continuous, ...defined(palette) }, values);
}

/** MVS's `objMerge`: `undefined` is a missing key, `null` a value. */
function defined<T extends object>(o: T): Partial<T> {
  const out: Partial<T> = {};
  for (const key in o) if (o[key] !== undefined) out[key] = o[key];
  return out;
}

/** Distinct values in order of first occurrence, as strings (`getDistinctValuesInField`). */
function distinct(values: Iterable<PaletteValue>, upper: boolean): string[] {
  const seen = new Set<string>();
  for (const value of values) {
    if (isMissing(value)) continue;
    const s = upper ? `${value}`.toUpperCase() : `${value}`;
    seen.add(s);
  }
  return [...seen];
}

function named(name: string): readonly Color[] | Readonly<Record<string, Color>> | undefined {
  if (Object.hasOwn(ColorLists, name)) return ColorLists[name as ColorListName];
  if (Object.hasOwn(ColorDictionaries, name)) return ColorDictionaries[name as ColorDictionaryName];
  return undefined;
}

function categorical(p: Required<CategoricalPalette>, values: Iterable<PaletteValue>): PaletteFunction {
  const map = new Map<string, Color>();
  const source = typeof p.colors === "string" ? (named(p.colors) ?? []) : p.colors;
  const decode = (c: Color | ColorString) => (typeof c === "number" ? c : Color.decodePrimary(c));
  const keys = distinct(values, p.case_insensitive);
  if (p.sort === "lexical") keys.sort();
  else if (p.sort === "numeric") keys.sort((a, b) => Number.parseFloat(a) - Number.parseFloat(b));
  if (p.sort_direction === "descending") keys.reverse();
  if (!Array.isArray(source)) {
    for (const [key, color] of Object.entries(source as Record<string, Color | ColorString>)) {
      map.set(p.case_insensitive ? key.toUpperCase() : key, decode(color));
    }
  } else {
    const list = (source as readonly (Color | ColorString)[]).map(decode);
    let next = 0;
    for (const key of keys) {
      if (next < list.length) map.set(key, list[next]);
      next++;
      if (next >= list.length) {
        if (p.repeat_color_list) next = 0;
        else break;
      }
    }
  }
  const missing = p.missing_color === null ? undefined : Color.decodePrimary(p.missing_color);
  const f = (value: PaletteValue) => {
    // A missing value takes `missing_color` too: a kernel fallback draws "no category".
    if (isMissing(value)) return missing;
    const key = p.case_insensitive ? `${value}`.toUpperCase() : `${value}`;
    return map.get(key) ?? missing;
  };
  return Object.assign(f, { keys: keys.filter((k) => map.has(k)) });
}

/**
 * FNV-1a (32-bit) over the UTF-8 bytes of `s` (seqspec:format#palettes): where a streamed value outside a categorical
 * list palette's keys takes its color, the same in every region, viewer and reload. Allocates nothing.
 */
export function paletteHash(s: string): number {
  let h = 0x811c9dc5;
  const byte = (b: number) => {
    h = Math.imul(h ^ b, 0x01000193);
  };
  for (let i = 0; i < s.length; i++) {
    let c = s.charCodeAt(i);
    // A surrogate pair is one code point, four bytes; a lone surrogate is U+FFFD, as TextEncoder writes it.
    if (c >= 0xd800 && c < 0xdc00 && i + 1 < s.length && (s.charCodeAt(i + 1) & 0xfc00) === 0xdc00) {
      c = 0x10000 + ((c - 0xd800) << 10) + (s.charCodeAt(++i) - 0xdc00);
    } else if (c >= 0xd800 && c < 0xe000) c = 0xfffd;
    if (c < 0x80) byte(c);
    else if (c < 0x800) {
      byte(0xc0 | (c >> 6));
      byte(0x80 | (c & 63));
    } else if (c < 0x10000) {
      byte(0xe0 | (c >> 12));
      byte(0x80 | ((c >> 6) & 63));
      byte(0x80 | (c & 63));
    } else {
      byte(0xf0 | (c >> 18));
      byte(0x80 | ((c >> 12) & 63));
      byte(0x80 | ((c >> 6) & 63));
      byte(0x80 | (c & 63));
    }
  }
  return h >>> 0;
}

/** A categorical palette's list colors (a named list or a given list); null for a dictionary or `{value: Color}`. */
export function categoricalList(palette: CategoricalPalette): readonly Color[] | null {
  const colors = palette.colors ?? PaletteDefaults.categorical.colors;
  const source = typeof colors === "string" ? named(colors) : colors;
  if (!Array.isArray(source)) return null;
  return (source as readonly (Color | ColorString)[]).map((c) => (typeof c === "number" ? c : Color.decodePrimary(c)));
}

/** A value as a number (`parseFloat` of its string, as MVS). */
function toNumber(value: PaletteValue): number {
  if (typeof value === "number") return value;
  return isMissing(value) ? Number.NaN : Number.parseFloat(`${value}`);
}

function listColors(colors: string | readonly unknown[]): readonly unknown[] {
  return typeof colors === "string"
    ? Object.hasOwn(ColorLists, colors)
      ? ColorLists[colors as ColorListName]
      : []
    : colors;
}

function discrete(p: Required<DiscretePalette>, values: Iterable<PaletteValue>): PaletteFunction {
  let colors = listColors(p.colors);
  let bins: { color: Color | undefined; from: number; to: number }[] = [];
  const decode = (c: unknown) =>
    c === null ? undefined : typeof c === "number" ? c : Color.decodePrimary(c as string);
  if (colors.every((t) => !Array.isArray(t))) {
    if (p.reverse) colors = colors.slice().reverse();
    const section = 1 / colors.length;
    bins = colors.map((c, i) => ({ color: decode(c), from: i * section, to: (i + 1) * section }));
  } else if (colors.every((t) => Array.isArray(t) && t.length === 2)) {
    const pairs = colors as readonly (readonly [ColorString, number])[];
    bins = pairs.map(([c, from], i) => ({ color: decode(c), from, to: pairs[i + 1]?.[1] ?? Number.POSITIVE_INFINITY }));
  } else if (colors.every((t) => Array.isArray(t) && t.length === 3)) {
    const ranges = colors as readonly (readonly [ColorString | null, number | null, number | null])[];
    bins = ranges.map(([c, from, to]) => ({
      color: decode(c),
      from: from ?? Number.NEGATIVE_INFINITY,
      to: to ?? Number.POSITIVE_INFINITY,
    }));
  }
  if (bins.length === 0) return () => undefined;
  // As the continuous one: the scale inline, the bins' edges in typed arrays and the value read from `held`, so `at`
  // allocates nothing on a frame (a draw's `palette(mean)` per position or column).
  const from = Float64Array.from(bins, (x) => x.from);
  const to = Float64Array.from(bins, (x) => x.to);
  const colorOf = bins.map((x) => x.color);
  const range = numericRange(p, values);
  const held = new Float64Array(1);
  const core = (): Color | undefined => {
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
    for (let i = from.length - 1; i >= 0; i--) if (from[i] <= x && x <= to[i]) return colorOf[i];
    return undefined;
  };
  const at = (o: { readonly value: number }) => {
    held[0] = o.value;
    return core();
  };
  return Object.assign(
    (value: PaletteValue) => {
      held[0] = toNumber(value);
      return core();
    },
    { at },
  );
}

function continuous(p: Required<ContinuousPalette>, values: Iterable<PaletteValue>): PaletteFunction {
  let colors = listColors(p.colors);
  let stops: { color: Color; checkpoint: number }[];
  const decode = (c: unknown) => (typeof c === "number" ? c : Color.decodePrimary(c as string));
  if (colors.every((t) => !Array.isArray(t))) {
    if (p.reverse) colors = colors.slice().reverse();
    const n = colors.length - 1;
    stops = colors.map((c, i) => ({ color: decode(c), checkpoint: i / n }));
  } else {
    // `reverse` has no effect on explicit checkpoints, as in MVS.
    stops = (colors as readonly (readonly [ColorString, number])[]).map(([c, checkpoint]) => ({
      color: decode(c),
      checkpoint,
    }));
  }
  if (stops.length === 0) return () => undefined;
  const flow = (option: "auto" | ColorString | null, pick: (a: number, b: number) => boolean) => {
    if (option === null) return undefined;
    if (option !== "auto") return Color.decodePrimary(option);
    return stops.reduce((a, b) => (pick(a.checkpoint, b.checkpoint) ? a : b)).color;
  };
  const underflow = flow(p.underflow_color, (a, b) => a < b);
  const overflow = flow(p.overflow_color, (a, b) => a > b);
  stops.sort((a, b) => a.checkpoint - b.checkpoint);
  const cps = Float64Array.from(stops, (s) => s.checkpoint);
  const colorsAt = Int32Array.from(stops, (s) => s.color);
  const range = numericRange(p, values);
  // Inline throughout: the scale, the first checkpoint at or past x (`findPredecessorIndex`), `Color.interpolate`.
  // The value is read from `held`: `at`'s loads see only its callers' locations (a load polymorphic over a double
  // field boxes it).
  const held = new Float64Array(1);
  const core = (): Color | undefined => {
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
    if (i === 0) return x === cps[0] ? colorsAt[0] : underflow;
    if (i === n) return overflow;
    const t = (x - cps[i - 1]) / (cps[i] - cps[i - 1]);
    const c1 = colorsAt[i - 1];
    const c2 = colorsAt[i];
    const r1 = (c1 >> 16) & 255;
    const g1 = (c1 >> 8) & 255;
    const b1 = c1 & 255;
    const r = r1 + (((c2 >> 16) & 255) - r1) * t;
    const g = g1 + (((c2 >> 8) & 255) - g1) * t;
    const b = b1 + ((c2 & 255) - b1) * t;
    return (r << 16) | (g << 8) | b;
  };
  const at = (o: { readonly value: number }) => {
    held[0] = o.value;
    return core();
  };
  return Object.assign(
    (value: PaletteValue) => {
      held[0] = toNumber(value);
      return core();
    },
    { at },
  );
}

/**
 * A `normalized` palette's [min, max] over `value_domain`, its `null`s from the values (0..1 without any), the value
 * scaled `(x - min) / (max - min)`; [NaN, NaN]: absolute.
 */
function numericRange(
  p: Required<DiscretePalette> | Required<ContinuousPalette>,
  values: Iterable<PaletteValue>,
): Float64Array {
  if (p.mode !== "normalized") return Float64Array.of(Number.NaN, Number.NaN);
  let [xMin, xMax] = p.value_domain;
  if (xMin === null || xMax === null) {
    const numbers = distinct(values, false)
      .map((v) => Number.parseFloat(v))
      .filter((x) => !Number.isNaN(x));
    if (numbers.length > 0) {
      xMin ??= numbers.reduce((a, b) => (a < b ? a : b));
      xMax ??= numbers.reduce((a, b) => (a > b ? a : b));
    } else {
      xMin ??= 0;
      xMax ??= 1;
    }
  }
  return Float64Array.of(xMin, xMax);
}
