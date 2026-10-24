/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Colors as numbers (0xRRGGBB), parsed from the document's `Color` strings: Mol*'s `decodeColor` and
 * `Color.interpolate` (`@molstar/core/util/color`). `isValid` is MVS's `ColorT` check (`isColorT`, author Adam Midlik).
 * A split color `a/b` decodes to its first half: SVS viewers draw that one (seqspec:format#positions); black on invalid
 * input, as MVS's `SplitColorProp.fromString`.
 */

import { Color as MolstarColor } from '@molstar/core/util/color/color';
import { decodeColor } from '@molstar/core/util/color/utils';
import { ColorNames } from './names.js';

/** 0xRRGGBB. */
export type Color = number;

/** MVS's `ColorT` hex form: `#rgb` or `#rrggbb`. */
const mvsHexRegex = /^#([0-9A-F]{3}){1,2}$/i;

export const Color = {
  toHex(color: Color): string {
    return `#${`000000${color.toString(16)}`.slice(-6)}`;
  },

  /** Mol*'s `decodeColor`: `#rgb`, `#rgba`, `#rrggbb`, `#rrggbbaa` (alpha dropped), a color name, `rgb(r, g, b)`. */
  decode(colorString: string | null | undefined): Color | undefined {
    return decodeColor(colorString);
  },

  /** A document `Color`'s first half (a split color's); undefined when invalid or not a string. */
  decodeFirst(colorString: unknown): Color | undefined {
    if (typeof colorString !== 'string') return undefined;
    const slash = colorString.indexOf('/');
    return Color.decode(slash >= 0 ? colorString.substring(0, slash) : colorString);
  },

  /** A document `Color` (split colors by their first half); black when invalid, as MVS's `SplitColorProp`, or not a
   * string. */
  decodePrimary(colorString: unknown): Color {
    return Color.decodeFirst(colorString) ?? 0x000000;
  },

  /** MVS's `ColorT` check: `#rgb`, `#rrggbb`, a Mol* color name (case-sensitive, as MVS), or `a/b` of two of them. */
  isValid(value: unknown): value is string {
    if (typeof value !== 'string') return false;
    const parts = value.split('/');
    if (parts.length > 2) return false;
    return parts.every((p) => mvsHexRegex.test(p) || ColorNames.has(p));
  },

  /** Mol*'s `Color.interpolate`. */
  interpolate(c1: Color, c2: Color, t: number): Color {
    return MolstarColor.interpolate(c1 as MolstarColor, c2 as MolstarColor, t);
  },
};

/** A list of colors written as hex digits (`"fff5f0 fee0d2 …"`), the compact form of the copied tables. */
export function colorList(src: string): readonly Color[] {
  return src
    .trim()
    .split(/\s+/)
    .map((hex) => Number.parseInt(hex, 16));
}

/** A dictionary written as groups, a color then its keys (`"#ff0080 helix H #ffc800 strand B E"`). */
export function colorDict(src: string): Readonly<Record<string, Color>> {
  const out: Record<string, Color> = {};
  let color = 0;
  for (const token of src.trim().split(/\s+/)) {
    if (token.startsWith('#')) color = Color.decodePrimary(token);
    else out[token] = color;
  }
  return out;
}
