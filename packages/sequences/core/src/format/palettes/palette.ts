/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * MVS's palettes: the same kinds, params, defaults and evaluation, so a value gets one color in SVS and in MVS
 * (seqspec:format#themes). One implementation, Mol*'s `@molstar/core/util/color/value-palette` (MVS's runtime evaluates
 * the same), over the spec's shape: the document's palette to its props with SVS's defaults, then its function.
 * Differences from MVS, all from the spec: missing values (`null`, `undefined`) and the strings `.` and `?` never take a
 * list color (MVS skips them only in CIF annotations) and are colored by a categorical palette's `missing_color` alone
 * (`cifMissing`); colors are the first half of a split color (SVS viewers draw that one); `repeat_color_list` defaults to
 * `true` (a list starts over rather than leaving the values past its end uncolored; MVS's default is `false`).
 */

import { makeValuePaletteFunction, valuePalettePropsFromSpec } from '@molstar/core/util/color/value-palette';
import { Color } from './color.js';
import { ColorDictionaries, type ColorDictionaryName } from './dictionaries.js';
import { type ColorListName, ColorLists } from './lists.js';

/** A document `Color` string (`#rgb`, `#rrggbb`, a Mol* color name, `a/b`). */
export type ColorString = string;

export interface CategoricalPalette {
  readonly kind: 'categorical';
  /** A list name, a dictionary name, colors, or `{value: Color}`. */
  readonly colors?:
    | ColorListName
    | ColorDictionaryName
    | readonly ColorString[]
    | Readonly<Record<string, ColorString>>;
  readonly repeat_color_list?: boolean;
  readonly sort?: 'none' | 'lexical' | 'numeric';
  readonly sort_direction?: 'ascending' | 'descending';
  readonly case_insensitive?: boolean;
  readonly missing_color?: ColorString | null;
}

export interface DiscretePalette {
  readonly kind: 'discrete';
  /** A list name, colors, `[Color, checkpoint]` pairs, or `[Color | null, from | null, to | null]` ranges. */
  readonly colors?:
    | ColorListName
    | readonly ColorString[]
    | readonly (readonly [ColorString, number])[]
    | readonly (readonly [ColorString | null, number | null, number | null])[];
  readonly reverse?: boolean;
  readonly mode?: 'normalized' | 'absolute';
  readonly value_domain?: readonly [number | null, number | null];
}

export interface ContinuousPalette {
  readonly kind: 'continuous';
  /** A list name, colors, or `[Color, checkpoint]` pairs. */
  readonly colors?: ColorListName | readonly ColorString[] | readonly (readonly [ColorString, number])[];
  readonly reverse?: boolean;
  readonly mode?: 'normalized' | 'absolute';
  readonly value_domain?: readonly [number | null, number | null];
  readonly underflow_color?: 'auto' | ColorString | null;
  readonly overflow_color?: 'auto' | ColorString | null;
}

export type Palette = CategoricalPalette | DiscretePalette | ContinuousPalette;

/** MVS's defaults (`CategoricalPaletteDefaults`, `DiscretePaletteDefaults`, `ContinuousPaletteDefaults`), but a list repeats. */
export const PaletteDefaults = {
  categorical: {
    kind: 'categorical',
    colors: 'Category10',
    repeat_color_list: true,
    sort: 'none',
    sort_direction: 'ascending',
    case_insensitive: false,
    missing_color: null,
  } satisfies Required<CategoricalPalette>,
  discrete: {
    kind: 'discrete',
    colors: 'YlGn',
    reverse: false,
    mode: 'normalized',
    value_domain: [null, null],
  } satisfies Required<DiscretePalette>,
  continuous: {
    kind: 'continuous',
    colors: 'YlGn',
    reverse: false,
    mode: 'normalized',
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

/**
 * The palette's function over `values`: every value of the field over the layer's data, annotations in `data`
 * order, rows in order (first occurrence and the normalized domain come from them).
 */
export function makePaletteFunction(palette: Palette, values: Iterable<PaletteValue>): PaletteFunction {
  const props = valuePalettePropsFromSpec(palette, PaletteDefaults, { dictionaries: ColorDictionaries });
  return makeValuePaletteFunction(props, values, { cifMissing: true });
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
  const source =
    typeof colors === 'string'
      ? Object.hasOwn(ColorLists, colors)
        ? ColorLists[colors as ColorListName]
        : null
      : colors;
  if (!Array.isArray(source)) return null;
  return (source as readonly (Color | ColorString)[]).map((c) => (typeof c === 'number' ? c : Color.decodePrimary(c)));
}
