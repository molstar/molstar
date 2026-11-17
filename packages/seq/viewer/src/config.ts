/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's configuration items of the first slice (_spec/README.md "Configuration"): one default object of the
// shape the models read. Where a behavior could go more than one way, it is an item here with the recommended default.

import { configItem, type SeqConfigItem } from '@molstar/seqx-core/context';

/** What the models and behaviors read; read-only after composition. */
export interface ViewerConfig {
  readonly foundations: { readonly preset: 'comfortable' | 'compact' };
  readonly density: {
    /** The least pixels a position needs for its letter to be drawn. */
    readonly letterMinPx: number;
    /** The least pixels a position needs for a dot. */
    readonly dotMinPx: number;
  };
  readonly bindings: { readonly preset: 'default' | 'genome-browser' };
  readonly store: {
    /** The loaded bytes the data store keeps (64 MB). */
    readonly cacheBytes: number;
    /** Requests in flight at most (6). */
    readonly inFlight: number;
    /** Bins per tile (1,024). */
    readonly tileBins: number;
  };
  readonly hover: {
    /** Project the hover onto every linked section (true); false: only sections of the hovered system. */
    readonly projectLinked: boolean;
    /** Close the hover before projecting: a nucleotide lights its whole codon and its residue (true, EC-30, EC-144). */
    readonly closeCodons: boolean;
  };
  readonly viewport: {
    /** The most pixels a position may take when zoomed in: the zoom's upper bound. */
    readonly maxPxPerPosition: number;
    /** Sections on equal axes share one viewport (true, EC-113); false: a viewport per section. */
    readonly groupEqualAxes: boolean;
    /**
     * When a section's width changes: keep its zoom about the center (`keep-zoom`, the default: a section hidden for
     * lack of width comes back as it was, EC-112) or keep its range (`keep-range`). A fitted viewport stays fitted.
     */
    readonly onResize: 'keep-zoom' | 'keep-range';
  };
  readonly follow: {
    /**
     * A move asked of a follower (followers never lead, EC-114): `ignore` it (the default) or `detach` the follower,
     * which then has its own viewport until it is attached again (`attach`, or a fit).
     */
    readonly followerMoves: 'ignore' | 'detach';
    /**
     * A follower's range over a cyclic origin: kept as the wrap, two pieces (`wrap`, the default), or the hull of
     * the image (`hull`).
     */
    readonly cyclicArc: 'wrap' | 'hull';
  };
  readonly navigate: {
    /** The zoom factor of one wheel notch (100 px of delta) with the zoom modifier, and of a key zoom. */
    readonly wheelZoomStep: number;
    /** The part of the viewport an arrow key pans. */
    readonly keyPanFraction: number;
    /** Wheel deltas within one frame are applied once, at the frame (true); false: at each event. */
    readonly coalesceWheel: boolean;
  };
  readonly theme: {
    /**
     * Where a mounted viewer takes its token table from: the computed style of its root (`root`, the default: the
     * `--msui-*` and `--seq-*` properties, redrawn when `data-msui-theme` changes) or the context's built-in table of
     * the `theme` item (`config`). A headless viewer always uses the config.
     */
    readonly source: 'root' | 'config';
  };
  readonly canvas: {
    /**
     * The canvas's height: the layout's content (`content`, the default: tracks below the page are not laid out, so
     * the content is measured against `maxHeight`) or a fixed number of CSS px.
     */
    readonly height: 'content' | number;
    /** The tallest the canvas grows with `height: 'content'`, CSS px. */
    readonly maxHeight: number;
    /** The device pixel ratio: the device's (`device`, the default) or a fixed one (an export, a test). */
    readonly dpr: 'device' | number;
    /** How long after the last move of the view it counts as settled (open domains computed, lanes packed), ms. */
    readonly settleMs: number;
  };
  readonly pointer: {
    /**
     * A resting pointer is hit-tested again once a frame while the view moves under it (EC-143), so the hover follows
     * the data (true, the default); false: the hover changes with the pointer only.
     */
    readonly rehit: boolean;
    /**
     * A trackpad's pinch arrives as a wheel with Ctrl held and small deltas: it is the `pinch` gesture (true, the
     * default); false: it is a wheel with Ctrl.
     */
    readonly pinchFromCtrlWheel: boolean;
    /** A pinch's zoom factor is `exp(-deltaY · pinchSpeed)`. */
    readonly pinchSpeed: number;
    /**
     * A Ctrl+wheel with a larger delta than this (px) is a mouse wheel, not a pinch (not on macOS, where Ctrl+wheel is
     * always a pinch: Ctrl is no zoom modifier there).
     */
    readonly pinchMaxDelta: number;
    /** Pixels a line of a wheel in line mode (`deltaMode` 1) scrolls: 16. A page (mode 2) is the viewport's own axis. */
    readonly lineHeight: number;
  };
}

export const DefaultViewerConfig: ViewerConfig = {
  foundations: { preset: 'comfortable' },
  density: { letterMinPx: 7, dotMinPx: 2 },
  bindings: { preset: 'default' },
  store: { cacheBytes: 64 * 1024 * 1024, inFlight: 6, tileBins: 1024 },
  hover: { projectLinked: true, closeCodons: true },
  viewport: { maxPxPerPosition: 32, groupEqualAxes: true, onResize: 'keep-zoom' },
  follow: { followerMoves: 'ignore', cyclicArc: 'wrap' },
  navigate: { wheelZoomStep: 1.25, keyPanFraction: 0.1, coalesceWheel: true },
  theme: { source: 'root' },
  canvas: { height: 'content', maxHeight: 2048, dpr: 'device', settleMs: 150 },
  pointer: { rehit: true, pinchFromCtrlWheel: true, pinchSpeed: 0.01, pinchMaxDelta: 50, lineHeight: 16 },
};

type DeepPartial<T> = { readonly [K in keyof T]?: T[K] extends object ? DeepPartial<T[K]> : T[K] };

/** The defaults with the given items over them, one level of groups deep. */
export function viewerConfig(over: DeepPartial<ViewerConfig> = {}): ViewerConfig {
  const out: Record<string, unknown> = {};
  for (const [group, items] of Object.entries(DefaultViewerConfig))
    out[group] = { ...items, ...(over as Record<string, object | undefined>)[group] };
  return out as unknown as ViewerConfig;
}

/** A config item's key: a group and an item (`hover.projectLinked`). */
export type ViewerConfigKey = {
  [G in keyof ViewerConfig & string]: `${G}.${keyof ViewerConfig[G] & string}`;
}[keyof ViewerConfig & string];

/** The value type of a config item's key. */
export type ViewerConfigValue<K extends ViewerConfigKey> = K extends `${infer G}.${infer I}`
  ? G extends keyof ViewerConfig
    ? I extends keyof ViewerConfig[G]
      ? ViewerConfig[G][I]
      : never
    : never
  : never;

/** The viewer's config items, one per setting: `viewerItem('hover.projectLinked')`, a pair `[item, value]` in a spec. */
export function viewerItem<K extends ViewerConfigKey>(key: K): SeqConfigItem<ViewerConfigValue<K>> {
  return configItem(key);
}

/** The config the viewer's items in `pairs` set over the defaults (items that are not the viewer's are ignored). */
export function viewerConfigOf(pairs: readonly (readonly [SeqConfigItem<any>, unknown])[] = []): ViewerConfig {
  const over: Record<string, Record<string, unknown>> = {};
  for (const [item, value] of pairs) {
    const dot = item.key.indexOf('.');
    if (dot < 0) continue;
    const group = item.key.slice(0, dot);
    if (!Object.hasOwn(DefaultViewerConfig, group)) continue;
    (over[group] ??= {})[item.key.slice(dot + 1)] = value;
  }
  return viewerConfig(over as DeepPartial<ViewerConfig>);
}
