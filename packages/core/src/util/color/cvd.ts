/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Color } from './color.js';
import { Oklab } from './spaces/oklab.js';

/**
 * Color-vision deficiency: how a color looks under protanopia, deuteranopia or tritanopia (Machado, Oliveira and
 * Fernandes 2009, "A physiologically-based model for simulation of color vision deficiency", severity 1.0, applied in
 * linear sRGB), and a check that a palette's colors stay apart under each (ΔE in Oklab).
 */
export type ColorVision = 'normal' | 'protanopia' | 'deuteranopia' | 'tritanopia';

export const ColorVisions: readonly ColorVision[] = ['normal', 'protanopia', 'deuteranopia', 'tritanopia'];

const Machado: Readonly<Record<Exclude<ColorVision, 'normal'>, readonly number[]>> = {
  protanopia: [0.152286, 1.052583, -0.204868, 0.114503, 0.786281, 0.099216, -0.003882, -0.048116, 1.051998],
  deuteranopia: [0.367322, 0.860646, -0.227968, 0.280085, 0.672501, 0.047413, -0.01182, 0.04294, 0.968881],
  tritanopia: [1.255528, -0.076749, -0.178779, -0.078411, 0.930809, 0.147602, 0.004733, 0.691367, 0.3039],
};

function toLinear(c: number): number {
  const s = c / 255;
  return s <= 0.04045 ? s / 12.92 : ((s + 0.055) / 1.055) ** 2.4;
}

function fromLinear(x: number): number {
  const v = Math.max(0, Math.min(1, x));
  return Math.round((v <= 0.0031308 ? 12.92 * v : 1.055 * v ** (1 / 2.4) - 0.055) * 255);
}

/** The color as seen under a deficiency (itself under `normal`). */
export function simulateColorVision(color: Color, vision: ColorVision): Color {
  if (vision === 'normal') return color;
  const m = Machado[vision];
  const r = toLinear((color >> 16) & 255);
  const g = toLinear((color >> 8) & 255);
  const b = toLinear(color & 255);
  return Color.fromRgb(
    fromLinear(m[0] * r + m[1] * g + m[2] * b),
    fromLinear(m[3] * r + m[4] * g + m[5] * b),
    fromLinear(m[6] * r + m[7] * g + m[8] * b),
  );
}

export interface ColorVisionCheck {
  /** Every pair at least `threshold` apart under every vision checked. */
  readonly ok: boolean;
  /** The closest pair: its vision, indices and ΔE (Infinity for fewer than two colors). */
  readonly worst: { readonly vision: ColorVision; readonly i: number; readonly j: number; readonly distance: number };
}

/**
 * Whether a palette's colors stay apart: every two at least `threshold` (ΔE in Oklab, default 0.06) under normal vision
 * and the three simulated deficiencies (or the visions given).
 */
export function checkColorVision(
  colors: readonly Color[],
  threshold = 0.06,
  visions: readonly ColorVision[] = ColorVisions,
): ColorVisionCheck {
  let worst = { vision: 'normal' as ColorVision, i: -1, j: -1, distance: Number.POSITIVE_INFINITY };
  const a = Oklab.zero();
  const b = Oklab.zero();
  for (const vision of visions) {
    const seen = colors.map((c) => simulateColorVision(c, vision));
    for (let i = 0; i < seen.length; i++) {
      Oklab.fromColor(a, seen[i]);
      for (let j = i + 1; j < seen.length; j++) {
        const d = Oklab.distance(a, Oklab.fromColor(b, seen[j]));
        if (d < worst.distance) worst = { vision, i, j, distance: d };
      }
    }
  }
  return { ok: worst.distance >= threshold, worst };
}
