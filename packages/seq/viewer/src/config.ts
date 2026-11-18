/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The viewer's configuration items of the first slice (_spec/README.md "Configuration"): one default object of the
// shape the models read. Where a behavior could go more than one way, it is an item here with the recommended default.

import { configItem, type SeqConfigItem } from '@molstar/seqx-core/context';
import { DefaultMaxBlockValues } from '@molstar/seqx-core/data';
import { DefaultSpecBounds, type SpecBounds } from '@molstar/seqx-core/lang';

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
  readonly selection: {
    /**
     * The dynamic nodes the selection names (its elements' owners, its systems) are held in the store while selected
     * (true): a selected transcript stays when its region scrolls away; false: only while shown.
     */
    readonly holdNodes: boolean;
    /**
     * What an extend (Shift+click, Shift+drag, Shift+arrows) changes: the range the anchor is in, the other ranges and
     * elements kept (`range`, the default), or the whole selection, which becomes that range (`selection`).
     */
    readonly extendScope: 'range' | 'selection';
  };
  readonly highlight: {
    /** A highlight channel's color when `colors` names none: a hex literal or `token:<name>` (`token:cat-5`). */
    readonly color: string;
    /** Colors by channel (`link-hover`: `token:hover`, the hover's ink, tinted). */
    readonly colors: Readonly<Record<string, string>>;
  };
  readonly cursor: {
    /**
     * The keys that move the keyboard cursor: a position back and on, and to the section above and below; Shift with
     * `prev` or `next` extends the selection from it (`ArrowLeft`, `ArrowRight`, `ArrowUp`, `ArrowDown`).
     */
    readonly keys: { readonly prev: string; readonly next: string; readonly up: string; readonly down: string };
    /** A cursor moved out of view pans its section to show it (true); false: it moves unseen. */
    readonly reveal: boolean;
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
     * A move asked of a follower (a pan, a zoom, a pinch): `lead`, the default, moves its leader through the follow's
     * inverse, so the follower moves as asked and every section in step; `ignore` it; or `detach` the follower, which
     * then has its own viewport until it is attached again (`attach`, or a fit).
     */
    readonly followerMoves: 'lead' | 'ignore' | 'detach';
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
  readonly address: {
    /** The most nodes or elements an address may name (10,000): more is an issue, never a truncation. */
    readonly limit: number;
    /** How deep below the root resolving an address may walk (64). */
    readonly maxDepth: number;
    /** How many refs a chain may pass through (8). */
    readonly refDepth: number;
  };
  readonly links: {
    /** The prefix of this viewer's keys in a link (`seq-`: `!seq-focus=/refs/dbd`); other keys are another viewer's. */
    readonly prefix: string;
    /** The longest link read, in characters (8,192). */
    readonly maxLength: number;
    /** The most actions one link runs (16). */
    readonly maxActions: number;
    /** The highlight channel a link without `channel` highlights on (`link`). */
    readonly highlightChannel: string;
    /** The channel a hover over a link (and the `address.hover` command) uses, its own (EC-182): `link-hover`. */
    readonly hoverChannel: string;
  };
  readonly revisions: {
    /** The recent roots the spec keeps by revision, for a host's rebase of a stale patch (64). */
    readonly keepRecent: number;
  };
  readonly edit: {
    /**
     * A stale patch (made against an older revision) is refused (false, the default: EC-81), or rebased through
     * `merge3` from its base root (true: the host opts in).
     */
    readonly rebaseStale: boolean;
  };
  readonly preview: {
    /**
     * A preview that no longer applies after the root changed is dropped, with a warning in the spec's issues
     * (`issue`, the default) or quietly (`silent`).
     */
    readonly dropped: 'issue' | 'silent';
  };
  /**
   * The bounds on what a spec or an answer sizes, each refused with a coded issue before the work starts: the spec's
   * (`SpecBounds`: size 64 MB, nesting 256 and queries 64, a column 2^31 − 1, issues kept 100 a node and 10 places a
   * code, held by the viewer's patches) and `maxBlockValues`, the most values (rows times columns) a table of an
   * answer may hold (2^24). The nested sources' depth is the resolve config's (`maxNestedDepth`, 4).
   */
  readonly bounds: SpecBounds & { readonly maxBlockValues: number };
  readonly copyLink: {
    /** The forms an address is copied in, in order: a ref naming exactly it, then ids and forms (`['ref', 'ids']`). */
    readonly prefer: readonly ('ref' | 'ids')[];
    /** The numbering labels are copied in (`ordinal`); one that does not read back falls back to it. */
    readonly numbering: string;
    /** The action a copied link carries (`focus`). */
    readonly action: 'focus' | 'reveal' | 'select' | 'highlight';
  };
}

export const DefaultViewerConfig: ViewerConfig = {
  foundations: { preset: 'comfortable' },
  density: { letterMinPx: 7, dotMinPx: 2 },
  bindings: { preset: 'default' },
  store: { cacheBytes: 64 * 1024 * 1024, inFlight: 6, tileBins: 1024 },
  hover: { projectLinked: true, closeCodons: true },
  selection: { holdNodes: true, extendScope: 'range' },
  highlight: { color: 'token:cat-5', colors: { 'link-hover': 'token:hover' } },
  cursor: { keys: { prev: 'ArrowLeft', next: 'ArrowRight', up: 'ArrowUp', down: 'ArrowDown' }, reveal: true },
  viewport: { maxPxPerPosition: 32, groupEqualAxes: true, onResize: 'keep-zoom' },
  follow: { followerMoves: 'lead', cyclicArc: 'wrap' },
  navigate: { wheelZoomStep: 1.25, keyPanFraction: 0.1, coalesceWheel: true },
  theme: { source: 'root' },
  canvas: { height: 'content', maxHeight: 2048, dpr: 'device', settleMs: 150 },
  pointer: { rehit: true, pinchFromCtrlWheel: true, pinchSpeed: 0.01, pinchMaxDelta: 50, lineHeight: 16 },
  address: { limit: 10_000, maxDepth: 64, refDepth: 8 },
  links: { prefix: 'seq-', maxLength: 8192, maxActions: 16, highlightChannel: 'link', hoverChannel: 'link-hover' },
  revisions: { keepRecent: 64 },
  edit: { rebaseStale: false },
  preview: { dropped: 'issue' },
  bounds: { ...DefaultSpecBounds, maxBlockValues: DefaultMaxBlockValues },
  copyLink: { prefer: ['ref', 'ids'], numbering: 'ordinal', action: 'focus' },
};

type DeepPartial<T> = {
  readonly [K in keyof T]?: T[K] extends readonly unknown[] ? T[K] : T[K] extends object ? DeepPartial<T[K]> : T[K];
};

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
