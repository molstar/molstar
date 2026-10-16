/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui/reactive/readonly-state';

// A map's levels (parity.md §2, v5's density controls): what people change most on a volume, its iso values and colors,
// without the representation's whole form. One control for a loaded map (its isosurfaces) and for streamed density (its
// channels): the map's values as a histogram, a marker per level to drag, each level's color and opacity.

/** The map's values: what σ is measured against, and how many voxels fall in each bin of [min, max]. */
export interface VolumeLevelStats {
  readonly min: number;
  readonly max: number;
  readonly mean: number;
  readonly sigma: number;
}

/** Voxels per bin over [min, max], equal bins. */
export interface VolumeLevelHistogram {
  readonly min: number;
  readonly max: number;
  readonly counts: ArrayLike<number>;
}

/** One level: an isosurface (or a streamed channel) at an iso value, in a color. */
export interface VolumeLevel {
  /** The representation's ref, or the channel's name. */
  readonly id: string;
  readonly label: string;
  /** The iso value in map units (σ is `(value - mean) / sigma`). */
  readonly value: number;
  /** Its uniform color; undefined: colored by a theme (a pick makes it uniform). */
  readonly color: number | undefined;
  /** 0 (clear) to 1. */
  readonly opacity: number;
  readonly hidden: boolean;
}

/** A change to a level: any of its settings. */
export interface VolumeLevelChange {
  /** In map units. */
  readonly value?: number;
  readonly color?: number;
  readonly opacity?: number;
}

export interface VolumeLevelsModel {
  /** Which map, when a row has more than one (streamed X-ray density: 2Fo-Fc, Fo-Fc). */
  readonly title?: string;
  /** The map's statistics; streamed density: the whole map's, from the server. */
  readonly stats: VolumeLevelStats;
  /** The map's values; undefined: not known (streamed density before it arrives). */
  readonly histogram: ReadonlyState<VolumeLevelHistogram | undefined>;
  /** Its levels, with the changes not applied yet already in (what is being dragged shows where it is). */
  readonly levels: ReadonlyState<readonly VolumeLevel[]>;
  /** How values are shown and typed: σ above the mean, or map units. Shared by every map's control. */
  readonly unit: ReadonlyState<VolumeLevelUnit>;
  setUnit(unit: VolumeLevelUnit): void;
  /**
   * Changes a level. Changes in quick succession (a drag, keys held, the color field dragged across) are one undo
   * step: applied one at a time as the viewer can, the latest winning, until `settle` or a moment's quiet.
   */
  change(id: string, change: VolumeLevelChange): void;
  /** The changes so far are one step; the next starts another (the pointer let go). */
  settle(): void;
}

export type VolumeLevelUnit = 'sigma' | 'absolute';
