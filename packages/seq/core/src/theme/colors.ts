/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '@molstar/core/util/color/color';

// A spec's colors (_spec/README.md "Themes", Colors): a hex literal, a token (`token:<name>`), or `light-dark(a, b)`.
// Colors here are 32-bit RGBA (`0xRRGGBBAA`), as the draw lists hold them; Mol*'s are RGB (`0xRRGGBB`).

export { isColorSpelling, parseColor, resolveColor, type TokenTable } from './spelling.js';

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
  // a direct-mapped memo: the contrast arithmetic boxes doubles, and a layout asks for the same few fills every frame
  const slot =
    (Math.imul(fill >>> 8, 0x9e3779b1) ^ Math.imul((dark >>> 8) ^ Math.imul(light >>> 8, 31), 0x85ebca6b)) >>> 22;
  if (InkFills[slot] === fill && InkDark[slot] === dark && InkLight[slot] === light && InkUsed[slot] === 1)
    return InkOut[slot];
  const ink = computeInk(fill, dark, light);
  InkFills[slot] = fill;
  InkDark[slot] = dark;
  InkLight[slot] = light;
  InkOut[slot] = ink;
  InkUsed[slot] = 1;
  return ink;
}

const InkSlots = 1024;
const InkFills = new Uint32Array(InkSlots);
const InkDark = new Uint32Array(InkSlots);
const InkLight = new Uint32Array(InkSlots);
const InkOut = new Uint32Array(InkSlots);
const InkUsed = new Uint8Array(InkSlots);

function computeInk(fill: number, dark: number, light: number): number {
  const cd = contrastRgba(fill, dark);
  const cl = contrastRgba(fill, light);
  const best = cd >= cl ? dark : light;
  if (Math.max(cd, cl) >= 4.5) return best;
  return contrastRgba(fill, Black) >= contrastRgba(fill, White) ? Black : White;
}
