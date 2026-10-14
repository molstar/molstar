/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Colors as numbers (0xRRGGBB), parsed from the document's `Color` strings.
 *
 * Copied from Mol* (MIT, (c) 2018-2026 mol* contributors; Mol*'s src): `mol-util/color/color.ts`
 * (`toHexStyle`, `fromHexStyle`, `interpolate`; authors Alexander Rose, David Sehnal), `mol-util/color/utils.ts`
 * (`decodeColor`; David Sehnal, Zach Charlop-Powers), `extensions/mvs/tree/mvs/param-types.ts` (`isColorT`) and
 * `extensions/mvs/components/split-uniform-color-theme.ts` (`SplitColorProp.fromString`, black on invalid input;
 * Adam Midlik). Mol*'s imports removed. A split color `a/b` decodes to its first half: SVS viewers draw that one
 * (seqspec:format#positions).
 */

import { ColorNames } from "./names.ts";

/** 0xRRGGBB. */
export type Color = number;

const hexColorRegex = /^#([0-9A-F]{3}|[0-9A-F]{4}|[0-9A-F]{6}|[0-9A-F]{8})$/i;
const rgbColorRegex = /^rgb\(\s*(\d{1,3})\s*,\s*(\d{1,3})\s*,\s*(\d{1,3})\s*\)$/i;
/** MVS's `ColorT` hex form: `#rgb` or `#rrggbb`. */
const mvsHexRegex = /^#([0-9A-F]{3}){1,2}$/i;

export const Color = {
  toHex(color: Color): string {
    return `#${`000000${color.toString(16)}`.slice(-6)}`;
  },

  /** Mol*'s `decodeColor`: `#rgb`, `#rgba`, `#rrggbb`, `#rrggbbaa` (alpha dropped), a color name, `rgb(r, g, b)`. */
  decode(colorString: string | null | undefined): Color | undefined {
    if (colorString === undefined || colorString === null) return undefined;
    let s = colorString;
    if (hexColorRegex.test(s)) {
      if (s.length === 4 || s.length === 5) s = `#${s[1]}${s[1]}${s[2]}${s[2]}${s[3]}${s[3]}`;
      else if (s.length === 9) s = s.substring(0, 7);
      const result = Number.parseInt(s.replace("#", "0x"), 16);
      if (!Number.isNaN(result)) return result;
    }
    const named = ColorNames.get(s.toLowerCase());
    if (named !== undefined) return named;
    const rgb = rgbColorRegex.exec(s);
    if (rgb) {
      const r = Number.parseInt(rgb[1], 10);
      const g = Number.parseInt(rgb[2], 10);
      const b = Number.parseInt(rgb[3], 10);
      if (r <= 255 && g <= 255 && b <= 255) return (r << 16) | (g << 8) | b;
    }
    return undefined;
  },

  /** A document `Color`'s first half (a split color's); undefined when invalid or not a string. */
  decodeFirst(colorString: unknown): Color | undefined {
    if (typeof colorString !== "string") return undefined;
    const slash = colorString.indexOf("/");
    return Color.decode(slash >= 0 ? colorString.substring(0, slash) : colorString);
  },

  /** A document `Color` (split colors by their first half); black when invalid, as MVS's `SplitColorProp`, or not a
   * string. */
  decodePrimary(colorString: unknown): Color {
    return Color.decodeFirst(colorString) ?? 0x000000;
  },

  /** MVS's `ColorT` check: `#rgb`, `#rrggbb`, a Mol* color name (case-sensitive, as MVS), or `a/b` of two of them. */
  isValid(value: unknown): value is string {
    if (typeof value !== "string") return false;
    const parts = value.split("/");
    if (parts.length > 2) return false;
    return parts.every((p) => mvsHexRegex.test(p) || ColorNames.has(p));
  },

  /** Mol*'s `Color.interpolate`. */
  interpolate(c1: Color, c2: Color, t: number): Color {
    const r1 = (c1 >> 16) & 255;
    const g1 = (c1 >> 8) & 255;
    const b1 = c1 & 255;
    const r = r1 + (((c2 >> 16) & 255) - r1) * t;
    const g = g1 + (((c2 >> 8) & 255) - g1) * t;
    const b = b1 + ((c2 & 255) - b1) * t;
    return (r << 16) | (g << 8) | b;
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
    if (token.startsWith("#")) color = Color.decodePrimary(token);
    else out[token] = color;
  }
  return out;
}
