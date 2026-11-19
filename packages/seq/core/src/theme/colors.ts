/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color/color';

// A spec's colors (_spec/README.md "Themes", Colors): a hex literal, a token (`token:<name>`), or `light-dark(a, b)`.
// Colors here are 32-bit RGBA (`0xRRGGBBAA`), as the draw lists hold them; Mol*'s are RGB (`0xRRGGBB`).

/** What resolving a color reads: the table's theme (`light`, `dark`) and its colors by name. */
export interface TokenTable {
  readonly theme: string;
  readonly colors: Readonly<Record<string, number>>;
  /** Equal tables have equal versions. */
  readonly version: number;
}

/** `#rgb`, `#rgba`, `#rrggbb`, `#rrggbbaa` as `0xRRGGBBAA`; undefined for anything else. */
export function parseColor(s: string | undefined): number | undefined {
  if (typeof s !== 'string' || !/^#([0-9a-f]{3,4}|[0-9a-f]{6}|[0-9a-f]{8})$/i.test(s)) return undefined;
  let h = s.slice(1);
  if (h.length <= 4) h = [...h].map((c) => c + c).join('');
  if (h.length === 6) h += 'ff';
  return Number.parseInt(h, 16) >>> 0;
}

const LightDark = /^light-dark\(\s*([^,()\s]+)\s*,\s*([^,()\s]+)\s*\)$/;

/**
 * A spec's color: a literal (`#rrggbb`), a token (`token:accent` reads `seq-accent`, then `msui-accent`, then `accent`),
 * or `light-dark(a, b)` (`a` in a light table, `b` in a dark one; each a literal or a token). Undefined when it is none
 * of these, or a token the table lacks.
 */
export function resolveColor(value: unknown, tokens: TokenTable): number | undefined {
  if (typeof value !== 'string') return undefined;
  if (value.startsWith('token:')) {
    const name = value.slice(6);
    const colors = tokens.colors;
    return colors[`seq-${name}`] ?? colors[`msui-${name}`] ?? colors[name];
  }
  if (value.startsWith('light-dark(')) {
    const m = LightDark.exec(value);
    return m ? resolveColor(tokens.theme === 'dark' ? m[2] : m[1], tokens) : undefined;
  }
  return parseColor(value);
}

/** Whether a string is spelled as a color (a literal, a token name, `light-dark` of two): the check's rule. */
export function isColorSpelling(value: string): boolean {
  if (value.startsWith('token:')) return /^token:[a-z][a-z0-9-]*$/.test(value);
  const m = LightDark.exec(value);
  if (m) return isColorSpelling(m[1]) && isColorSpelling(m[2]) && !m[1].startsWith('light-dark');
  return parseColor(value) !== undefined;
}

/** `#rrggbbaa`. */
export function colorHex(c: number): string {
  return `#${(c >>> 0).toString(16).padStart(8, '0')}`;
}

/** Mol*'s RGB color as RGBA, opaque. */
export function fromRgb(c: Color | number): number {
  return ((c << 8) | 0xff) >>> 0;
}

/** An RGBA color's RGB, as Mol*'s colors are. */
export function toRgb(c: number): Color {
  return Color(c >>> 8);
}

/** WCAG contrast of two RGBA colors (their alpha ignored): Mol*'s `Color.contrast`. */
export function contrastRgba(a: number, b: number): number {
  return Color.contrast(toRgb(a), toRgb(b));
}

/** Pure black and white, RGBA. */
export const Black = 0x000000ff;
export const White = 0xffffffff;

/**
 * Ink on a fill (_spec/README.md "Themes", Colors: ink follows its fill, 4.5:1): the table's dark or light ink,
 * whichever contrasts more; where neither reaches 4.5:1 (a mid-tone fill between inks that are not pure), pure black
 * or white, whichever is higher.
 */
export function inkFor(fill: number, dark: number, light: number): number {
  const cd = contrastRgba(fill, dark);
  const cl = contrastRgba(fill, light);
  const best = cd >= cl ? dark : light;
  if (Math.max(cd, cl) >= 4.5) return best;
  return contrastRgba(fill, Black) >= contrastRgba(fill, White) ? Black : White;
}
