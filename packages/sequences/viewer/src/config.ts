/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Config items (seqspec:viewer#config): typed settings with defaults, as Mol*'s PluginConfigItem, read by the managers
// at use.

import type { SeqViewerContext } from "./context.ts";
import type { RenderClamp } from "./layout.ts";

export interface SeqViewerConfigItem<T = unknown> {
  readonly key: string;
  readonly defaultValue: T;
  /** What it sets, for the reference. */
  readonly description: string;
}
export function SeqViewerConfigItem<T>(key: string, defaultValue: T, description = ""): SeqViewerConfigItem<T> {
  return { key, defaultValue, description };
}

const item = SeqViewerConfigItem;

/**
 * A track menu's built-ins (seqspec:interaction#the-default-ui): `rows` expands or collapses a split track's rows, only
 * on one; `colors` opens a layer's color schemes (seqspec:interaction#color-schemes-at-runtime), only where more than
 * one fits and `UI.ColorSchemes` is on.
 */
export type TrackActionId = "filter" | "hide" | "remove" | "focus" | "rows" | "colors" | "retry";
export interface CustomTrackAction {
  readonly id: string;
  readonly label: string;
  readonly danger?: boolean;
  /**
   * Absent: on every track. `split`, here and below: the value of the split row whose header the menu was opened on,
   * null on the track's own line (a host acts on that row alone, or on the whole track).
   */
  when?(ctx: SeqViewerContext, section: string, track: string, split: string | null): boolean;
  /** A check mark beside it (the choice it makes is the current one); absent: none. */
  checked?(ctx: SeqViewerContext, section: string, track: string, split: string | null): boolean;
  /** Its label for this track or split row; absent: `label`. */
  labelFor?(ctx: SeqViewerContext, section: string, track: string, split: string | null): string;
  run(ctx: SeqViewerContext, section: string, track: string, split: string | null): void;
}
export type TrackActionList = readonly (TrackActionId | CustomTrackAction)[] | "none";
/** How a UI shows the legend (`UI.Legend`, a UI spec's `legend`). */
export type LegendMode = "collapsed" | "expanded" | "hidden";

export interface PresetDefaults {
  readonly blocks: string;
  readonly columns: string;
  readonly text: string;
}

export const SeqViewerConfig = {
  Viewport: {
    SettleMs: item(
      "viewport.settleMs",
      150,
      "ms after the last interaction before the viewport settles and commits its view",
    ),
    AnimationMs: item("viewport.animationMs", 250, "an animated focus or zoom, ms"),
    OnAxisChange: item<"keepCenter" | "fit">(
      "viewport.onAxisChange",
      "keepCenter",
      "when a section's axis changes: keep the view's center, or fit the whole axis",
    ),
    ReducedMotion: item<boolean | null>(
      "viewport.reducedMotion",
      null,
      "true: no animation; null: the engine reads `prefers-reduced-motion`",
    ),
  },
  Input: {
    DragSlopPx: item("input.dragSlopPx", 3, "CSS px a press moves before it is a drag"),
    DoubleClickMs: item("input.doubleClickMs", 300, "the longest gap between the clicks of a double click, ms"),
  },
  Tasks: { BusyDelayMs: item("tasks.busyDelayMs", 750, "ms a task runs before the UI shows it as busy") },
  Issues: {
    PerNode: item(
      "issues.perNode",
      100,
      "issues a graph node keeps, the most severe first; the rest are counted in one `seqv.issues.truncated`",
    ),
  },
  Stream: {
    SampleBudget: item<number | null>(
      "stream.sampleBudget",
      null,
      "samples per view at most, where the source lists a level that coarse; null: one per device pixel of the data area",
    ),
    MaxInFlight: item("stream.maxInFlight", 6, "stream requests in flight at once"),
    CancelGraceMs: item(
      "stream.cancelGraceMs",
      250,
      "ms a request no view needs any more stays in flight before it is aborted (sooner when its slot is wanted); 0: at once",
    ),
    CacheBytes: item("stream.cacheBytes", 128 * 2 ** 20, "bytes of streamed regions held before eviction"),
    Margin: item("stream.margin", 1, "view widths requested on each side of the view"),
    LettersPerPx: item(
      "stream.lettersPerPx",
      16,
      "a letters stream is asked only where a device px holds at most this many positions",
    ),
    RegionBytes: item<number | null>(
      "stream.regionBytes",
      1 << 20,
      "a data request's bound on what its provider reads (`max_bytes`); past it the provider answers with counts",
    ),
    ShimmerFps: item("stream.shimmerFps", 30, "frames a second of a pending region's shimmer; 0: a still band"),
    Prefetch: item(
      "stream.prefetch",
      true,
      "request ahead: a view width past the margin on the side last panned toward, a zoom's next level, and the regions of the snapshots next to the one shown",
    ),
  },
  Text: { CacheEntries: item("text.cacheEntries", 20_000, "measured strings the text measurer keeps") },
  Presets: {
    Default: item<PresetDefaults>(
      "presets.default",
      { blocks: "tracks", columns: "msa", text: "wrapped" },
      "the preset by section kind (blocks, alignment columns, text layout) when neither the UI nor the section's `style` names one",
    ),
    Alias: item<Readonly<Record<string, string>>>("presets.alias", {}, "a section `style` key → a preset name"),
  },
  Render: {
    Clamp: item<RenderClamp | null>(
      "render.clamp",
      null,
      "bounds on the numeric render props, over presets and styles; null: none",
    ),
  },
  History: { Limit: item("history.limit", 200, "undo steps kept") },
  /** Read by the default UIs when they mount. */
  UI: {
    Chrome: item<"full" | "clean">(
      "ui.chrome",
      "clean",
      "the default layout's composition: `clean`, the canvas (hover, selection, the tooltip and hover card, the navigator), the player, the description, search, busy and the states; `full` adds the tabs, section headers, track menus, legend, export and issues",
    ),
    Legend: item<LegendMode>(
      "ui.legend",
      "collapsed",
      "the legend `collapsed` (its toggle only), `expanded` or `hidden`; a UI's own option, the person's toggle and `view.legend` go first",
    ),
    LegendItems: item("ui.legendItems", 24, 'items a legend group shows before its "+n more"'),
    LegendGroups: item("ui.legendGroups", 8, 'groups the legend shows before its "+n more"'),
    TrackActions: item<TrackActionList>(
      "ui.trackActions",
      ["rows", "filter", "focus", "colors", "hide", "retry"],
      "the track menu's items in order; `none`: no menu (`remove` is built in but not listed by default)",
    ),
    ColorSchemes: item("ui.colorSchemes", true, "the color scheme pickers (a legend group's title, a track's Colors)"),
    ColorSchemeWrites: item(
      "ui.colorSchemeWrites",
      false,
      "true: a picked scheme is written into the document; false: UI state only",
    ),
    ItemActions: item<readonly string[]>(
      "ui.itemActions",
      [],
      "the registered item actions a hover card offers, in order; empty: no card",
    ),
    HoverCardDelayMs: item(
      "ui.hoverCardDelay",
      400,
      "ms the pointer rests on an item with actions before its card opens",
    ),
    OwnRowHover: item<"item" | "residue">(
      "ui.ownRowHover",
      "residue",
      "a bar on a split row on a unit of its own: `item` hovers the whole item, `residue` the residue under the pointer",
    ),
  },
} as const;

export interface ConfigManager {
  get<T>(item: SeqViewerConfigItem<T>): T;
  /** Read by its manager at its next use. */
  set<T>(item: SeqViewerConfigItem<T>, value: T): void;
}

/**
 * Values by an item's key, not its object: a spec made with another copy of the viewer (its own `SeqViewerConfig`
 * objects, as registration accepts either copy's spec) sets the same items.
 */
export function createConfigManager(initial: readonly (readonly [SeqViewerConfigItem, unknown])[] = []): ConfigManager {
  const values = new Map<string, unknown>(initial.map(([i, v]) => [i.key, v]));
  return {
    get: <T>(i: SeqViewerConfigItem<T>) => (values.has(i.key) ? (values.get(i.key) as T) : i.defaultValue),
    set: (i, value) => void values.set(i.key, value),
  };
}
