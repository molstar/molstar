/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SchemeEntry } from '@molstar/seqx-core/extension';

// `rondo:hues` (_spec/README.md, Colour): n hues evenly round the Oklch wheel at the categorical tokens' lightness and
// chroma (L 0.68, C 0.13), the dark variant mapped into the dark band; **not colour-blind safe** (24 hues: adjacent
// ones 15° apart), so its substitute is `categories`. Named by the showcase's spec, never the layer's default. The
// conversions are written here: the extension imports the core's extension entry alone.

const L = 0.68;
const C = 0.13;
/** Lightness from the light band into the dark one (the core's derived variants: [0.4, 1] → [0.3, 0.62]). */
const darkL = (l: number) => 0.3 + ((l - 0.4) / 0.6) * 0.32;

/** sRGB channel from linear, clamped. */
function gamma(x: number): number {
  const c = Math.max(0, Math.min(1, x));
  return c <= 0.0031308 ? 12.92 * c : 1.055 * c ** (1 / 2.4) - 0.055;
}

/** An Oklch color as `0xRRGGBB` (Björn Ottosson's Oklab, clamped to the sRGB gamut). */
export function oklch(l: number, c: number, hDeg: number): number {
  const h = (hDeg * Math.PI) / 180;
  const a = c * Math.cos(h);
  const b = c * Math.sin(h);
  const l_ = (l + 0.3963377774 * a + 0.2158037573 * b) ** 3;
  const m_ = (l - 0.1055613458 * a - 0.0638541728 * b) ** 3;
  const s_ = (l - 0.0894841775 * a - 1.291485548 * b) ** 3;
  const r = gamma(4.0767416621 * l_ - 3.3077115913 * m_ + 0.2309699292 * s_);
  const g = gamma(-1.2684380046 * l_ + 2.6097574011 * m_ - 0.3413193965 * s_);
  const bl = gamma(-0.0041960863 * l_ - 0.7034186147 * m_ + 1.707614701 * s_);
  return (Math.round(r * 255) << 16) | (Math.round(g * 255) << 8) | Math.round(bl * 255);
}

/** `n` hues evenly round the wheel from 30° (light theme's), as `#rrggbb`. */
export function hueList(n: number, dark = false): string[] {
  const out: string[] = [];
  for (let i = 0; i < n; i++) {
    const c = oklch(dark ? darkL(L) : L, C, 30 + (360 * i) / Math.max(1, n));
    out.push(`#${c.toString(16).padStart(6, '0')}`);
  }
  return out;
}

/** The hues as draw colors (`0xRRGGBBAA`, opaque). */
export function hueColors(n: number, dark = false): Uint32Array {
  const out = new Uint32Array(n);
  for (let i = 0; i < n; i++) out[i] = ((oklch(dark ? darkL(L) : L, C, 30 + (360 * i) / n) << 8) | 0xff) >>> 0;
  return out;
}

/** The hues a genome of 24 chromosomes takes by default (`n` its param). */
export const DefaultHueCount = 24;

/** `rondo:hues`: a categorical palette of `n` hues over a field (`system` when absent), in its domain's order. */
export const RondoHues: SchemeEntry = {
  name: 'rondo:hues',
  label: 'Hues (24)',
  subject: 'categories',
  cvdSafe: false,
  substitute: 'categories',
  dark: 'derive',
  expand: (p) => ({
    kind: 'field',
    params: {
      name: typeof p.field === 'string' ? p.field : 'system',
      palette: {
        kind: 'categorical',
        params: {
          colors: hueList(typeof p.n === 'number' && p.n > 0 ? Math.floor(p.n) : DefaultHueCount),
          repeat: false,
        },
      },
    },
  }),
};
