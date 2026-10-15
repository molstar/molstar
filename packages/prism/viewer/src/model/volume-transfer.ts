/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';
import type { VolumeLevelHistogram, VolumeLevelStats, VolumeLevelUnit } from './volume-levels.js';

// A direct volume's transfer function (parity.md §2, v5's control points editor, improved): opacity and color over the
// map's values on one axis. The opacity is the representation's control points; the color is its theme's: the
// volume-value theme's stops (their offsets are along the theme's domain, placed on this axis as the renderer does), or
// one color. Both edit in place, a drag one undo step.

/** A control point: `x` from 0 to 1 over the map's [min, max], `alpha` the opacity per step (small: 0.05 is a lot). */
export interface TransferPoint {
  readonly x: number;
  readonly alpha: number;
}

/** A color stop on the same axis as the points. */
export interface TransferStop {
  readonly x: number;
  readonly color: number;
}

export type TransferColors =
  /** The volume-value theme: its stops, and where its domain lies on the axis (the stops move within it). */
  | {
      readonly kind: 'value';
      readonly stops: readonly TransferStop[];
      readonly interpolate: boolean;
      readonly domain: readonly [number, number];
    }
  | { readonly kind: 'uniform'; readonly color: number }
  /** Another theme: shown as its name; its colors are not on the axis. */
  | { readonly kind: 'theme'; readonly label: string };

export interface TransferChange {
  readonly points?: readonly TransferPoint[];
  readonly stops?: readonly TransferStop[];
  readonly uniform?: number;
}

export interface TransferPreset {
  readonly label: string;
  readonly points: readonly TransferPoint[];
}

/** Colors for the peaks: categorical (in order) or sequential (spread over the peaks, low to high). */
export interface TransferPalette {
  readonly id: string;
  readonly label: string;
  readonly colors: readonly number[];
  readonly sequential?: boolean;
}

/** The colors a palette gives `n` peaks, low to high. */
export function paletteColors(palette: TransferPalette, n: number): number[] {
  const c = palette.colors;
  if (!palette.sequential || n <= 1 || c.length < 2) {
    return Array.from({ length: n }, (_, i) => c[(palette.sequential ? c.length - 1 : i) % c.length]);
  }
  return Array.from({ length: n }, (_, i) => c[Math.round((i * (c.length - 1)) / (n - 1))]);
}

export interface VolumeTransferModel {
  /** Which representation, when a volume has more than one. */
  readonly title?: string;
  readonly stats: VolumeLevelStats;
  readonly histogram: ReadonlyState<VolumeLevelHistogram | undefined>;
  /** In order of `x`, with the changes not applied yet already in. */
  readonly points: ReadonlyState<readonly TransferPoint[]>;
  readonly colors: ReadonlyState<TransferColors>;
  /** The opacity axis: logarithmic from this up to 1 (Mol*'s `yAxis`), as opacities per step are small. */
  readonly logMin: number;
  /** Curves made from the map's statistics (Mol*'s presets). */
  presets(): readonly TransferPreset[];
  readonly unit: ReadonlyState<VolumeLevelUnit>;
  setUnit(unit: VolumeLevelUnit): void;
  /** Changes in quick succession are one undo step (a drag), as the levels' (`VolumeLevelsModel.change`). */
  change(change: TransferChange): void;
  settle(): void;
  /** Palettes to color the peaks with (`colorsAtPeaks`). */
  readonly palettes: readonly TransferPalette[];
  /**
   * Color stops at the curve's peaks, at each one's ends (a peak is one color, the colors blend between peaks): a
   * palette's colors (by id), or with none the theme's own colors in order (not its whites and greys), then distinct
   * ones; on another theme, the volume-value theme with them.
   */
  colorsAtPeaks(palette?: string): void;
  /**
   * A preset's curve, and the color stops moved onto its peaks with the colors they have (on the volume-value theme), so
   * the colors stay with the parts of the map the curve shows: one step.
   */
  applyPreset(points: readonly TransferPoint[]): void;
}

/** A peak of the curve: a run of points above zero, from its first to its last. */
export interface TransferPeak {
  readonly from: number;
  readonly to: number;
}

/** The curve's peaks: each run of points above zero. */
export function transferPeaks(points: readonly TransferPoint[]): TransferPeak[] {
  const peaks: TransferPeak[] = [];
  let run: TransferPoint[] = [];
  const flush = () => {
    if (run.length) peaks.push({ from: run[0].x, to: run[run.length - 1].x });
    run = [];
  };
  for (const p of points) {
    if (p.alpha > 0) run.push(p);
    else flush();
  }
  flush();
  return peaks;
}

/** The color painted at `x` (as the renderer does: the theme's domain on the axis, clamped at its ends). */
export function transferColorAt(colors: TransferColors, x: number): number | undefined {
  if (colors.kind === 'uniform') return colors.color;
  if (colors.kind !== 'value' || !colors.stops.length) return undefined;
  const [d0, d1] = colors.domain;
  const t = d1 > d0 ? Math.min(1, Math.max(0, (x - d0) / (d1 - d0))) : 0;
  const stops = colors.stops;
  const offset = (s: TransferStop) => (d1 > d0 ? (s.x - d0) / (d1 - d0) : 0);
  if (t <= offset(stops[0])) return stops[0].color;
  for (let i = 1; i < stops.length; i++) {
    const a = offset(stops[i - 1]);
    const b = offset(stops[i]);
    if (t > b) continue;
    if (!colors.interpolate) return t < (a + b) / 2 ? stops[i - 1].color : stops[i].color;
    return mix(stops[i - 1].color, stops[i].color, b > a ? (t - a) / (b - a) : 1);
  }
  return stops[stops.length - 1].color;
}

function mix(c1: number, c2: number, t: number): number {
  const r = Math.round(((c1 >> 16) & 255) * (1 - t) + ((c2 >> 16) & 255) * t);
  const g = Math.round(((c1 >> 8) & 255) * (1 - t) + ((c2 >> 8) & 255) * t);
  const b = Math.round((c1 & 255) * (1 - t) + (c2 & 255) * t);
  return (r << 16) | (g << 8) | b;
}
