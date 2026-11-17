/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ink of a canvas (stage §2): what is drawn over a canvas (a toolbar, a label, a band) takes its theme from the
// canvas's background color, not from the UI's theme, so it stays readable on whatever the app draws. Generic: a CSS
// color in, a theme out. DOM-free.

/** The theme of the chrome over a canvas: `dark` chrome (light text) on a dark background, `light` on a light one. */
export type InkTheme = 'light' | 'dark';

/** What a control drawn over a canvas needs: the canvas's color (a CSS color) and the chrome's theme on it. */
export interface CanvasInk {
  readonly theme: InkTheme;
  readonly color: string;
}

/** The relative luminance (WCAG 2, 0..1) of a `0xRRGGBB` color. */
export function luminance(color: number): number {
  return (
    0.2126 * channel((color >> 16) & 0xff) + 0.7152 * channel((color >> 8) & 0xff) + 0.0722 * channel(color & 0xff)
  );
}

function channel(c: number): number {
  const v = c / 255;
  return v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4;
}

/** Where black and white text have the same contrast against a background (WCAG: (L + 0.05)² = 1.05 × 0.05). */
export const InkThreshold = 0.179;

export function themeForLuminance(l: number): InkTheme {
  return l < InkThreshold ? 'dark' : 'light';
}

export function themeForColor(color: number): InkTheme {
  return themeForLuminance(luminance(color));
}

/** A few named colors an app sets a background to; others are not parsed (the ink then counts them as white). */
const Named: Record<string, number> = { white: 0xffffff, black: 0x000000 };

/**
 * Parses `#rgb`, `#rrggbb`, `#rrggbbaa`, `rgb()` / `rgba()` (comma or space separated) and `white` / `black` into
 * `0xRRGGBB` (alpha ignored); undefined otherwise.
 */
export function parseColor(text: string): number | undefined {
  const t = text.trim().toLowerCase();
  if (t in Named) return Named[t];
  const hex = /^#?([0-9a-f]{3}|[0-9a-f]{6}|[0-9a-f]{8})$/.exec(t)?.[1];
  if (hex) {
    const six = hex.length === 3 ? [...hex].map((c) => c + c).join('') : hex.slice(0, 6);
    return Number.parseInt(six, 16);
  }
  const rgb = /^rgba?\(\s*(\d+(?:\.\d+)?)[\s,]+(\d+(?:\.\d+)?)[\s,]+(\d+(?:\.\d+)?)/.exec(t);
  if (rgb) {
    const c = (s: string) => Math.max(0, Math.min(255, Math.round(Number(s))));
    return (c(rgb[1]) << 16) | (c(rgb[2]) << 8) | c(rgb[3]);
  }
  return undefined;
}

/** `0xRRGGBB` as `#rrggbb`. */
export function colorHex(color: number): string {
  return `#${(color & 0xffffff).toString(16).padStart(6, '0')}`;
}

/** The ink on a background color (a CSS color); a color that does not parse counts as white (light chrome). */
export function inkForColor(color: string): CanvasInk {
  const parsed = parseColor(color);
  return parsed === undefined ? { theme: 'light', color: '#ffffff' } : { theme: themeForColor(parsed), color };
}
