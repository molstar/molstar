/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from '../color.js';

export { Oklab };

/**
 * Oklab (Björn Ottosson, 2020): a perceptual color space where a step reads as the same step anywhere, so palettes
 * interpolate in it, distances in it compare colors (ΔE), and its lightness orders a ramp.
 *
 * - L [0..1]: lightness
 * - a, b [about -0.4..0.4]: green to red, blue to yellow
 */
interface Oklab extends Array<number> {
  [d: number]: number;
  '@type': 'oklab';
  length: 3;
}

function Oklab() {
  return Oklab.zero();
}

namespace Oklab {
  export function zero(): Oklab {
    return [0, 0, 0] as unknown as Oklab;
  }

  export function create(l: number, a: number, b: number): Oklab {
    const out = zero();
    out[0] = l;
    out[1] = a;
    out[2] = b;
    return out;
  }

  function toLinear(c: number): number {
    const s = c / 255;
    return s <= 0.04045 ? s / 12.92 : ((s + 0.055) / 1.055) ** 2.4;
  }

  function fromLinear(x: number): number {
    return x <= 0.0031308 ? 12.92 * x : 1.055 * x ** (1 / 2.4) - 0.055;
  }

  /** From linear sRGB (each 0..1). */
  export function fromLinearRgb(out: Oklab, r: number, g: number, b: number): Oklab {
    const l = Math.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b);
    const m = Math.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b);
    const s = Math.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b);
    out[0] = 0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s;
    out[1] = 1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s;
    out[2] = 0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s;
    return out;
  }

  /** To linear sRGB, unclamped (out of gamut where a channel leaves 0..1). */
  export function toLinearRgb(out: [number, number, number], lab: Oklab): [number, number, number] {
    const l = (lab[0] + 0.3963377774 * lab[1] + 0.2158037573 * lab[2]) ** 3;
    const m = (lab[0] - 0.1055613458 * lab[1] - 0.0638541728 * lab[2]) ** 3;
    const s = (lab[0] - 0.0894841775 * lab[1] - 1.291485548 * lab[2]) ** 3;
    out[0] = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s;
    out[1] = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s;
    out[2] = -0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s;
    return out;
  }

  export function fromColor(out: Oklab, color: Color): Oklab {
    return fromLinearRgb(out, toLinear((color >> 16) & 255), toLinear((color >> 8) & 255), toLinear(color & 255));
  }

  const _rgb: [number, number, number] = [0, 0, 0];

  /** The color, each channel clamped into the gamut. */
  export function toColor(lab: Oklab): Color {
    toLinearRgb(_rgb, lab);
    const ch = (x: number) => Math.round(Math.max(0, Math.min(1, fromLinear(Math.max(0, x)))) * 255);
    return Color.fromRgb(ch(_rgb[0]), ch(_rgb[1]), ch(_rgb[2]));
  }

  /** Whether the color is inside the sRGB gamut (a small tolerance for rounding). */
  export function inGamut(lab: Oklab): boolean {
    toLinearRgb(_rgb, lab);
    const e = 1e-4;
    return _rgb.every((x) => x >= -e && x <= 1 + e);
  }

  /** Euclidean distance: Oklab's ΔE. */
  export function distance(a: Oklab, b: Oklab): number {
    const x = a[0] - b[0];
    const y = a[1] - b[1];
    const z = a[2] - b[2];
    return Math.sqrt(x * x + y * y + z * z);
  }

  export function interpolate(out: Oklab, a: Oklab, b: Oklab, t: number): Oklab {
    out[0] = a[0] + (b[0] - a[0]) * t;
    out[1] = a[1] + (b[1] - a[1]) * t;
    out[2] = a[2] + (b[2] - a[2]) * t;
    return out;
  }

  /** The chroma scaled down (hue and lightness kept) until the color is in gamut: a bisection, 16 steps. */
  export function clampChroma(out: Oklab, lab: Oklab): Oklab {
    out[0] = Math.max(0, Math.min(1, lab[0]));
    out[1] = lab[1];
    out[2] = lab[2];
    if (inGamut(out)) return out;
    let lo = 0;
    let hi = 1;
    for (let i = 0; i < 16; i++) {
      const k = (lo + hi) / 2;
      out[1] = lab[1] * k;
      out[2] = lab[2] * k;
      if (inGamut(out)) lo = k;
      else hi = k;
    }
    out[1] = lab[1] * lo;
    out[2] = lab[2] * lo;
    return out;
  }

  const _a = zero();
  const _b = zero();

  /** ΔE between two colors. */
  export function colorDistance(a: Color, b: Color): number {
    return distance(fromColor(_a, a), fromColor(_b, b));
  }

  /** A color's Oklab lightness, 0 to 1. */
  export function lightness(c: Color): number {
    return fromColor(_a, c)[0];
  }

  /** Interpolation of two colors in Oklab (`t` in 0 to 1), the result in gamut. */
  export function interpolateColors(a: Color, b: Color, t: number): Color {
    fromColor(_a, a);
    fromColor(_b, b);
    return toColor(clampChroma(_a, interpolate(_a, _a, _b, t)));
  }

  /**
   * The color with its lightness mapped from `from` into `to` (linearly, clamped), hue kept and chroma clamped to the
   * gamut: a dark theme's variant of a light theme's color (`from` the light band, `to` the dark one).
   */
  export function mapLightness(c: Color, from: readonly [number, number], to: readonly [number, number]): Color {
    fromColor(_a, c);
    const t = from[1] === from[0] ? 0.5 : (_a[0] - from[0]) / (from[1] - from[0]);
    _a[0] = to[0] + Math.max(0, Math.min(1, t)) * (to[1] - to[0]);
    return toColor(clampChroma(_a, _a));
  }

  /** The color with its lightness set (hue kept, chroma clamped). */
  export function withLightness(c: Color, l: number): Color {
    fromColor(_a, c);
    _a[0] = l;
    return toColor(clampChroma(_a, _a));
  }
}
