/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { configItem } from '@molstar/seqx-core/extension';

// The extension's settings (_spec/README.md, Configuration): how stems and ribbons are built and bounded, keys under
// its namespace (`rondo.*`). The look shared with bridges is core's `ribbons.*` (`RenderConfig`).

export interface RondoConfig {
  /** A stem's target width, CSS px of the layout at the ring's outer circumference (never device px). */
  readonly stemPx: number;
  /** The neck's depth inward from the mouth, a share of the hole's radius. */
  readonly stemDepth: number;
  /** The widest stem in view's neck, CSS px: one px a count for the view. */
  readonly maxWidthPx: number;
  /** The items gate and the drawn bound in one. */
  readonly maxRibbons: number;
  /** A stub's length, a share of the hole's radius. */
  readonly stubShare: number;
  /** The least distance between two labelled stub tips, CSS px (their boxes never overlap either). */
  readonly stubLabelPx: number;
  /** The hit grid's cell over the hole, CSS px. */
  readonly hitCellPx: number;
  /** The wash over the rest of the hole when a stem or a block is hovered (0: none). */
  readonly fade: number;
  /** Segments a cubic is flattened into (drawn and hit alike). */
  readonly flatten: number;
}

export const DefaultRondoConfig: RondoConfig = {
  stemPx: 16,
  stemDepth: 0.12,
  maxWidthPx: 24,
  maxRibbons: 1024,
  stubShare: 0.2,
  stubLabelPx: 24,
  hitCellPx: 16,
  fade: 0.6,
  flatten: 12,
};

/** The settings as items (`rondo.stemPx`, …): what a spec's `config` names. */
export const RondoConfigItems = {
  stemPx: configItem<number>('rondo.stemPx'),
  stemDepth: configItem<number>('rondo.stemDepth'),
  maxWidthPx: configItem<number>('rondo.maxWidthPx'),
  maxRibbons: configItem<number>('rondo.maxRibbons'),
  stubShare: configItem<number>('rondo.stubShare'),
  stubLabelPx: configItem<number>('rondo.stubLabelPx'),
  hitCellPx: configItem<number>('rondo.hitCellPx'),
  fade: configItem<number>('rondo.fade'),
  flatten: configItem<number>('rondo.flatten'),
} as const;

/** The settings with some overridden (a host's `rondoWith`). */
export function rondoConfig(over: Partial<RondoConfig> = {}): RondoConfig {
  return { ...DefaultRondoConfig, ...over };
}
