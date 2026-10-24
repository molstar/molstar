/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from './color.js';
import { decodeColor } from './utils.js';

/** Utility functions for handling "split color" strings with primary and secondary color, e.g. 'red/white', '#ff0000/#ffffff'.
 * Regular color strings are treated as the same color repeated (e.g. 'red' == 'red/red'). */
export const SplitColor = {
  /** Return number-encoded colors before and after slash, or the same value twice if there is no slash.
   * e.g. 'red/white' -> [0xff0000, 0xffffff], 'red' -> [0xff0000, 0xff0000],
   * Throw error on invalid input. */
  decodeStrict(colorString: string): [Color, Color] {
    const out = SplitColor.decode(colorString);
    if (out === undefined) throw new Error(`Invalid color "${colorString}".`);
    return out;
  },
  /** Return number-encoded colors before and after slash, or the same value twice if there is no slash.
   * e.g. 'red/white' -> [0xff0000, 0xffffff], 'red' -> [0xff0000, 0xff0000],
   * Return `undefined` on invalid input. */
  decode(colorString: string): [Color, Color] | undefined {
    const out: [Color | undefined, Color | undefined] = [undefined, undefined];
    SplitColor.decodeTo(colorString, out);
    if (out[0] !== undefined && out[1] !== undefined) {
      return out as [Color, Color];
    } else {
      return undefined;
    }
  },
  /** Parse `colorString` into number-encoded colors before and after slash, or the same value twice if there is no slash.
   * Write result into `out`.
   * e.g. 'red/white' -> [0xff0000, 0xffffff], 'red' -> [0xff0000, 0xff0000],
   * Write [`undefined`, `undefined` on invalid input. */
  decodeTo(colorString: string, out: [Color | undefined, Color | undefined]): void {
    const slashIdx = colorString.indexOf('/');
    if (slashIdx >= 0) {
      out[0] = decodeColor(colorString.substring(0, slashIdx));
      out[1] = decodeColor(colorString.substring(slashIdx + 1, undefined));
    } else {
      out[0] = out[1] = decodeColor(colorString);
    }
  },
  /** Convert to hex code string, '#ff0000/#ffffff', '#ff0000' */
  toHexStyle(primaryColor: Color, secondaryColor?: Color): string {
    if (secondaryColor === undefined || secondaryColor === primaryColor) {
      return Color.toHexStyle(primaryColor);
    } else {
      return `${Color.toHexStyle(primaryColor)}/${Color.toHexStyle(secondaryColor)}`;
    }
  },
};

/** A regular color or a split color, as the values of `SplitColorProp.PD()` (`value-palette-params.ts`). */
export type SplitColorProp =
  | { name: 'oneColor'; params: { color: Color } }
  | { name: 'splitColor'; params: { color1: Color; color2: Color } };
/** A regular color, a split color or no color, as the values of `SplitColorProp.PDOptional()`. */
export type OptionalSplitColorProp = SplitColorProp | NoneSplitColorProp;
export type NoneSplitColorProp = { name: 'none'; params: {} };
type ColorTuple<T extends OptionalSplitColorProp> = T extends SplitColorProp ? [Color, Color] : [undefined, undefined];

const FALLBACK_COLOR = Color(0x000000);

const _tmpColors: [Color | undefined, Color | undefined] = [undefined, undefined];

/** The values of a split color param, without its param definition (`SplitColorProp` in `value-palette-params.ts` adds it). */
export const SplitColorValue = {
  /** Invalid colors are black; `null` is no color. */
  fromString<T extends string | null>(colorString: T): T extends string ? SplitColorProp : NoneSplitColorProp {
    if (colorString == null) {
      return SplitColorValue.none() satisfies NoneSplitColorProp as any;
    }
    SplitColor.decodeTo(colorString, _tmpColors);
    if (_tmpColors[0] === _tmpColors[1]) {
      return { name: 'oneColor', params: { color: _tmpColors[0] ?? FALLBACK_COLOR } } satisfies SplitColorProp as any;
    } else {
      return {
        name: 'splitColor',
        params: { color1: _tmpColors[0] ?? FALLBACK_COLOR, color2: _tmpColors[1] ?? FALLBACK_COLOR },
      } satisfies SplitColorProp as any;
    }
  },
  none(): NoneSplitColorProp {
    return { name: 'none', params: {} };
  },
  toString(value: OptionalSplitColorProp): string {
    if (value.name === 'none') return 'None';
    if (value.name === 'oneColor') return Color.toHexStyle(value.params.color);
    else return `${Color.toHexStyle(value.params.color1)}/${Color.toHexStyle(value.params.color2)}`;
  },
  toTuple<T extends OptionalSplitColorProp>(value: T): ColorTuple<T> {
    if (value.name === 'none') return [undefined, undefined] as ColorTuple<T>;
    if (value.name === 'oneColor') return [value.params.color, value.params.color] as ColorTuple<T>;
    else return [value.params.color1, value.params.color2] as ColorTuple<T>;
  },
};
