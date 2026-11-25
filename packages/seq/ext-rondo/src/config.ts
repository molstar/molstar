/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { configItem, settingOf } from '@molstar/seqx-core/extension';

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
  /** A thin ribbon's or stub's hit: within this of its centre line, CSS px (either side). */
  readonly hitSlopPx: number;
  /** The wash over the rest of the hole when a stem or a block is hovered (0: none). */
  readonly fade: number;
  /** Segments a cubic is flattened into (drawn and hit alike). */
  readonly flatten: number;
  /**
   * The stems' opacity on a ring of one chromosome or a pair, a share of the ribbons' (the information is in the
   * ribbons and stubs there; the trumpets would dominate).
   */
  readonly focusStems: number;
  /** A stub of fewer links than this folds with the stem's other small stubs into one "other chromosomes" stub. */
  readonly stubMin: number;
}

export const DefaultRondoConfig: RondoConfig = {
  stemPx: 16,
  stemDepth: 0.12,
  maxWidthPx: 24,
  maxRibbons: 1024,
  stubShare: 0.2,
  stubLabelPx: 24,
  hitCellPx: 16,
  hitSlopPx: 4,
  fade: 0.6,
  flatten: 12,
  focusStems: 0.35,
  stubMin: 2,
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
  hitSlopPx: configItem<number>('rondo.hitSlopPx'),
  fade: configItem<number>('rondo.fade'),
  flatten: configItem<number>('rondo.flatten'),
  focusStems: configItem<number>('rondo.focusStems'),
  stubMin: configItem<number>('rondo.stubMin'),
} as const;

/** The settings with some overridden (a host's `rondoWith`). */
export function rondoConfig(over: Partial<RondoConfig> = {}): RondoConfig {
  return { ...DefaultRondoConfig, ...over };
}

/** The settings as composed, by the descriptor's numbers they override (kept: one object a pair). */
const Composed = new WeakMap<RondoConfig, WeakMap<object, RondoConfig>>();

/**
 * The settings a layout reads (_spec/README.md, Configuration): the composition's `config` (`rondo.*` items, read
 * from the render config a layer type is given, `settingOf`) over the descriptor's (`rondoWith`'s, else the defaults).
 * The same object while neither changes: a warm layout makes none.
 */
export function composedConfig(base: RondoConfig, render: { readonly settings?: Readonly<Record<string, unknown>> }) {
  const settings = render.settings;
  if (!settings) return base;
  let byBase = Composed.get(base);
  if (!byBase) Composed.set(base, (byBase = new WeakMap()));
  let out = byBase.get(settings);
  if (out) return out;
  const over: Record<string, number> = {};
  for (const [name, item] of Object.entries(RondoConfigItems)) {
    const v = settingOf(render, item);
    if (typeof v === 'number' && Number.isFinite(v)) over[name] = v;
  }
  out = { ...base, ...over };
  byBase.set(settings, out);
  return out;
}
