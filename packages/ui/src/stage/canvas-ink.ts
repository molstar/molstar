/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ColorNames } from '@molstar/core/util/color/names';
import { Hsl } from '@molstar/core/util/color/spaces/hsl';
import { Rgb } from '@molstar/core/util/color/spaces/rgb';

// The ink of a canvas (stage §2): what is drawn over a canvas (a toolbar, a label, a band) takes its theme from the
// canvas's background color, not from the UI's theme, so it stays readable on whatever the app draws. Generic: a CSS
// color in, a theme out. DOM-free (a color it cannot read is the view's to resolve, `CanvasInk.resolved`).

/** The theme of the chrome over a canvas: `dark` chrome (light text) on a dark background, `light` on a light one. */
export type InkTheme = 'light' | 'dark';

/** What a control drawn over a canvas needs: the canvas's color (a CSS color) and the chrome's theme on it. */
export interface CanvasInk {
  readonly theme: InkTheme;
  /** The app's color as given (the stage's background). */
  readonly color: string;
  /**
   * Whether `theme` was read from the color. False for a color this module does not parse (`oklch()`, `var()`) or one
   * not opaque (alpha below 1: what is behind it shows through): `theme` is then `light`, and a view resolves it
   * through the CSS engine (stage §2).
   */
  readonly resolved: boolean;
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

const byte = (v: number) => Math.max(0, Math.min(255, Math.round(v)));
const NumberRe = /^[+-]?(\d+\.?\d*|\.\d+)(e[+-]?\d+)?$/;

/** A channel: a number (0..255) or a percentage. */
function rgbChannel(s: string): number | undefined {
  const percent = s.endsWith('%');
  const n = percent ? s.slice(0, -1) : s;
  if (!NumberRe.test(n)) return undefined;
  return byte(percent ? (Number(n) * 255) / 100 : Number(n));
}

/** An alpha: a number (0..1) or a percentage; absent, 1. */
function alphaOf(s: string | undefined): number | undefined {
  if (s === undefined) return 1;
  const percent = s.endsWith('%');
  const n = percent ? s.slice(0, -1) : s;
  if (!NumberRe.test(n)) return undefined;
  return Math.max(0, Math.min(1, percent ? Number(n) / 100 : Number(n)));
}

/** A hue in degrees (`deg`, `turn`, `rad`, `grad` or a bare number). */
function hueOf(s: string): number | undefined {
  const m = /^(.*?)(deg|turn|rad|grad)?$/.exec(s);
  if (!m || !NumberRe.test(m[1])) return undefined;
  const v = Number(m[1]);
  const deg = m[2] === 'turn' ? v * 360 : m[2] === 'rad' ? (v * 180) / Math.PI : m[2] === 'grad' ? v * 0.9 : v;
  return ((deg % 360) + 360) % 360;
}

function fraction(s: string): number | undefined {
  const n = s.endsWith('%') ? s.slice(0, -1) : s;
  if (!NumberRe.test(n)) return undefined;
  return Math.max(0, Math.min(1, Number(n) / 100));
}

const hsl = Hsl();
const rgb = Rgb();

/**
 * A CSS color as `[0xRRGGBB, alpha]`: `#rgb`, `#rgba`, `#rrggbb`, `#rrggbbaa`; `rgb()` / `rgba()` and `hsl()` / `hsla()`
 * (comma or space separated, channels as numbers or percentages, an alpha after `/` or as a fourth); a CSS color name;
 * `transparent`. Undefined otherwise (`oklch()`, `color()`, `var()`: the view's to resolve).
 */
export function parseColorAlpha(text: string): readonly [color: number, alpha: number] | undefined {
  const t = text.trim().toLowerCase();
  if (t === 'transparent') return [0, 0];
  const hex = /^#([0-9a-f]{3,4}|[0-9a-f]{6}|[0-9a-f]{8})$/.exec(t)?.[1];
  if (hex) {
    const full = hex.length <= 4 ? [...hex].map((c) => c + c).join('') : hex;
    const alpha = full.length === 8 ? Number.parseInt(full.slice(6), 16) / 255 : 1;
    return [Number.parseInt(full.slice(0, 6), 16), alpha];
  }
  const fn = /^(rgba?|hsla?)\(\s*([^()]*?)\s*\)$/.exec(t);
  if (fn) {
    const [body, alphaPart] = fn[2].split(/\s*\/\s*/);
    const parts = body.split(/\s*,\s*|\s+/).filter(Boolean);
    if (alphaPart === undefined ? parts.length !== 3 && parts.length !== 4 : parts.length !== 3) return undefined;
    const alpha = alphaOf(alphaPart ?? parts[3]);
    if (alpha === undefined) return undefined;
    if (fn[1].startsWith('rgb')) {
      const r = rgbChannel(parts[0]);
      const g = rgbChannel(parts[1]);
      const b = rgbChannel(parts[2]);
      return r === undefined || g === undefined || b === undefined ? undefined : [(r << 16) | (g << 8) | b, alpha];
    }
    const h = hueOf(parts[0]);
    const s = fraction(parts[1]);
    const l = fraction(parts[2]);
    if (h === undefined || s === undefined || l === undefined) return undefined;
    hsl[0] = h;
    hsl[1] = s;
    hsl[2] = l;
    // Rounded per channel (`Hsl.toColor` truncates: 127.5 is 0x7f).
    Hsl.toRgb(rgb, hsl);
    return [(byte(rgb[0] * 255) << 16) | (byte(rgb[1] * 255) << 8) | byte(rgb[2] * 255), alpha];
  }
  return Object.hasOwn(ColorNames, t) ? [ColorNames[t as keyof typeof ColorNames] as number, 1] : undefined;
}

/** `0xRRGGBB` as `#rrggbb`. */
export function colorHex(color: number): string {
  return `#${(color & 0xffffff).toString(16).padStart(6, '0')}`;
}

/**
 * The ink on a background color (a CSS color, kept as given). A color that does not parse, or one not opaque (what is
 * behind it shows through: a nearly transparent black is the page's color, not black), is not resolved: light chrome
 * until a view composites it over what is behind it through the CSS engine (stage §2).
 */
export function inkForColor(color: string): CanvasInk {
  const parsed = parseColorAlpha(color);
  return parsed === undefined || parsed[1] < 1
    ? { theme: 'light', color, resolved: false }
    : { theme: themeForColor(parsed[0]), color, resolved: true };
}
