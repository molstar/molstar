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
import {
  DefaultRowsZoom,
  DefaultZoomStrip,
  type PositionsRailMode,
  type RailPlace,
  type RowsRailMode,
} from '@molstar/seqx-core/render';
import type { ZoomOverBlocks } from '@molstar/seqx-core/render/window';
import { DefaultStoryConfig, type StoryConfig } from '@molstar/seqx-core/story';

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
    /**
     * The hover card (the `hoverCard` behavior): what the pointer rests on, in words, beside it (true); false: none (the
     * live region still says it).
     */
    readonly card: boolean;
    /**
     * The items the hover card lists at most (the features at a place, a group's members): past it the first ones and a
     * row `+k more` (8); the card's data keeps the rest (`HoverCard.more`) and the live region says every one.
     */
    readonly cardItems: number;
    /**
     * The lines a card row's value wraps to within the card's width, then an ellipsis (3): a tooltip's text, a long INFO
     * value. The default adapter's; a host's card reads it from the show's options.
     */
    readonly tooltipLines: number;
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
    /**
     * What a drag or an extend across a collapsed break selects: the positions shown between its ends (`shown`, the
     * default: what is selected is what is seen), or `span`, the hidden stretches too (a gene's whole span to export).
     */
    readonly acrossBreaks: 'shown' | 'span';
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
  };
  readonly reveal: {
    /** A cursor moved out of view brings its section to it, the layout's least move (true); false: it moves unseen. */
    readonly cursor: boolean;
    /**
     * A selection made elsewhere (another section, a structure) that projects onto a section with none of it in view
     * scrolls its first part in, as a move by hand (no history step): in text sections (`text`, the default: a text
     * hides most of itself, a window shows the projection where it is), in `all` sections, or `none`.
     */
    readonly linked: 'text' | 'all' | 'none';
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
  /**
   * The rows' zoom (_spec/README.md "Configuration", Rows; core's _spec Rows, The rows rail): the least and the most px a
   * row (`minPx` 1/64, `maxPx` 40: the rows' zoom keys and the rows rail's grips stop there), and whether the rows zoom at
   * all (`enabled`, true; false: the rows rail is a pan-only bar, no grips, the zoom keys inert, the wheel and the rail's
   * body still scroll). A source's own floor (a lazy alignment serving no row levels: a device pixel a row, below which
   * its band is blank) applies whatever `minPx` says: the effective floor is the larger of the two.
   */
  readonly rows: {
    readonly zoom: { readonly minPx: number; readonly maxPx: number; readonly enabled: boolean };
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
  /**
   * Navigators (_spec/README.md "Configuration", Navigators): the default navigator of every linear section (`rail`, a
   * ruler rail; `none`), what a drag from outside a window does (`sweep` a new window and zoom to it; `pan`: centre
   * there, then pan), when a drag pushes a step back (`far`: farther than one window width; `never`, `always`), what a
   * double-click does (`toggle-fit`: fit, or fitted, zoom so the window is `zoomPx` navigator px wide; `fit`, `none`),
   * which way a navigator runs (`display`: the bound section's display order, a drag right moves the view right;
   * `genomic`), and when a navigator gets a zoom strip (`zoomStrip`: its reachable span narrower than `minPx`, 48; the
   * strip `heightPx` more, 16; over blocks `overBlocks`: `always`, the chromosome the control row at any window width,
   * `narrow` only for a narrow window, `never`), and when a rows band draws its rows rail (`rowsRail`: `always`, beside a declared
   * overview too; `unlessOverview`, replaced by a navigator over its rows; `none`), when a linear section draws its
   * positions rail (`positionsRail`: `always`, beside a declared navigator too; `unlessNavigator`, replaced by one;
   * `none`) and where (`railPlace`: `ruler` below the ruler, `bottom` under the last track, `navigator` by a declared
   * horizontal navigator, else `ruler`; `pinned` at the frame's foot, in reach at any scroll; one place or a list, a
   * rail each on the one window).
   */
  readonly navigator: {
    readonly form: 'rail' | 'none';
    readonly rowsRail: RowsRailMode;
    readonly positionsRail: PositionsRailMode;
    readonly railPlace: RailPlace | readonly RailPlace[];
    readonly dragOutside: 'sweep' | 'pan';
    readonly dragStep: 'far' | 'never' | 'always';
    readonly doubleClick: 'toggle-fit' | 'fit' | 'none';
    readonly zoomPx: number;
    /** A square's muted line of its gestures under it (true); false: none (a figure's, a host's own words). */
    readonly squareHint: boolean;
    readonly direction: 'display' | 'genomic';
    readonly zoomStrip: { readonly minPx: number; readonly heightPx: number; readonly overBlocks: ZoomOverBlocks };
  };
  readonly view: {
    /** The steps back the view's history of jumps keeps (50): the oldest goes first. */
    readonly historyLimit: number;
    /** A double-click focusing what the view already shows focused goes back (true); false: it does nothing. */
    readonly backOnFocused: boolean;
    /** How long a jump (focus, fit, reveal, back, forward, a step zoom) animates, ms (240); 0: a cut. */
    readonly jumpMs: number;
    /**
     * Whether motion plays, for jumps and stories alike (one item): `platform` (the default) animates unless the
     * platform asks for reduced motion (then a cut), `cut` always a cut, `animate` always animated.
     */
    readonly motion: 'platform' | 'cut' | 'animate';
    /**
     * The keys of the view's commands, by command id: a chord is a key with `Shift+`, `Alt+`, `Mod+` (Cmd on macOS, Ctrl
     * elsewhere) before it; a printable key matches whatever Shift it needs (`+`). The first is the command's shortcut.
     * An override merges by command id (`[]`: none). Read when the viewer is made.
     */
    readonly keys: Readonly<Record<string, readonly string[]>>;
  };
  readonly navigate: {
    /** The zoom factor of one wheel notch (100 px of delta) with the zoom modifier, and of a key zoom. */
    readonly wheelZoomStep: number;
    /** The zoom factor a double-click's zoom out (Shift+double-click) takes about the pointer (2). */
    readonly doubleClickZoomStep: number;
    /** The part of the viewport an arrow key pans. */
    readonly keyPanFraction: number;
    /** Wheel deltas within one frame are applied once, at the frame (true); false: at each event. */
    readonly coalesceWheel: boolean;
    /**
     * On a row of blocks (a genome's chromosomes): one viewport over the blocks end to end, a pan carrying across a gap
     * into the next (`cross`, the default), or a pan, a zoom and a fit kept in the block the window is in (`stop`: one
     * chromosome at a time; a jump to another block moves there).
     */
    readonly acrossBlocks: 'cross' | 'stop';
  };
  readonly theme: {
    /**
     * Where a mounted viewer takes its token table from: the computed style of its root (`root`, the default: the
     * `--msui-*` and `--seq-*` properties, redrawn when `data-msui-theme` changes) or the context's built-in table of
     * the `theme` item (`config`). A headless viewer always uses the config.
     */
    readonly source: 'root' | 'config';
    /**
     * The reader's color-blind substitution, at start (off): each scheme that is not safe drawn as its declared
     * substitute, a categorical palette as Okabe–Ito; the legend says so. A reader's setting, never the spec's
     * (`theme.cvd-safe` toggles it).
     */
    readonly cvdSafe: boolean;
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
  /**
   * Stories (the core's `StoryConfig`): a step's transition (`duration` 240 ms, `easing` ease in and out) unless the
   * step says, its `linger` (2,000 ms), the view's path (`auto`), the color space (`oklab`), the cut (0.5), the start of a step (`spec`), the encoder's
   * form (`patch`) and its chain (8).
   */
  readonly story: StoryConfig;
}

export const DefaultViewerConfig: ViewerConfig = {
  foundations: { preset: 'comfortable' },
  density: { letterMinPx: 7, dotMinPx: 2 },
  bindings: { preset: 'default' },
  store: { cacheBytes: 64 * 1024 * 1024, inFlight: 6, tileBins: 1024 },
  hover: { projectLinked: true, closeCodons: true, card: true, cardItems: 8, tooltipLines: 3 },
  selection: { holdNodes: true, extendScope: 'range', acrossBreaks: 'shown' },
  highlight: { color: 'token:cat-5', colors: { 'link-hover': 'token:hover' } },
  cursor: { keys: { prev: 'ArrowLeft', next: 'ArrowRight', up: 'ArrowUp', down: 'ArrowDown' } },
  reveal: { cursor: true, linked: 'text' },
  viewport: { maxPxPerPosition: 32, groupEqualAxes: true, onResize: 'keep-zoom' },
  rows: { zoom: DefaultRowsZoom },
  follow: { followerMoves: 'lead', cyclicArc: 'wrap' },
  navigator: {
    form: 'rail',
    rowsRail: 'always',
    positionsRail: 'always',
    railPlace: 'ruler',
    dragOutside: 'sweep',
    dragStep: 'far',
    doubleClick: 'toggle-fit',
    zoomPx: 48,
    squareHint: true,
    direction: 'display',
    zoomStrip: DefaultZoomStrip,
  },
  view: {
    historyLimit: 50,
    backOnFocused: true,
    jumpMs: 240,
    motion: 'platform',
    keys: {
      'view.back': ['Alt+ArrowLeft', 'Backspace'],
      'view.forward': ['Alt+ArrowRight'],
      'view.fit': ['0'],
      'view.zoom-in': ['=', '+'],
      'view.zoom-out': ['-'],
      'view.pan-left': ['Shift+Alt+ArrowLeft'],
      'view.pan-right': ['Shift+Alt+ArrowRight'],
      'view.pan-up': ['Shift+Alt+ArrowUp'],
      'view.pan-down': ['Shift+Alt+ArrowDown'],
      'selection.copy-letters': ['Mod+c'],
      'selection.bookmark': ['Mod+d'],
      'section.reverse': ['Shift+R'],
      // the next and the previous block of a row (a chromosome): Alt+Shift would type a character on macOS (EC-155)
      'block.next': ['Alt+PageDown'],
      'block.previous': ['Alt+PageUp'],
      'section.collapse': ['Shift+C'],
      // each matches only its own `e.key`: typed alone on the common layouts, never a browser's zoom key (EC-155)
      'rows.zoom-in': [']'],
      'rows.zoom-out': ['['],
      // a page of rows: the window moves by the rows it shows less one, the cursor shown on them with it (the rows'
      // own `rows.page-*` take no key: the cursor's scroll the rows when it is not shown)
      'cursor.page-down': ['PageDown'],
      'cursor.page-up': ['PageUp'],
      // the scroller under the pointer to its start or end (a rows or lanes band, a section past its cap)
      'cursor.scroll-start': ['Home'],
      'cursor.scroll-end': ['End'],
      // the first and the last row (the page's keys when there is no cursor on rows)
      'cursor.first-row': ['Mod+ArrowUp'],
      'cursor.last-row': ['Mod+ArrowDown'],
      // the keyboard's way to the column titles' click: the next column, ascending then descending
      'rows.order-next': ['Shift+O'],
      // a narrow frame's collapsed tree opened, or closed to its toggle again (the toggle's click)
      'rows.tree': ['Shift+T'],
      // the keys panel: what the view answers to (a printable key: whatever Shift types it)
      'help.keys': ['?'],
      // a track header's menu at the keyboard cursor (the `headers` behavior's; Enter on a header too)
      'headers.menu': ['ContextMenu', 'Shift+F10'],
    },
  },
  navigate: {
    wheelZoomStep: 1.25,
    doubleClickZoomStep: 2,
    keyPanFraction: 0.1,
    coalesceWheel: true,
    acrossBlocks: 'cross',
  },
  theme: { source: 'root', cvdSafe: false },
  canvas: { height: 'content', maxHeight: 2048, dpr: 'device', settleMs: 150 },
  pointer: { rehit: true, pinchFromCtrlWheel: true, pinchSpeed: 0.01, pinchMaxDelta: 50, lineHeight: 16 },
  address: { limit: 10_000, maxDepth: 64, refDepth: 8 },
  links: { prefix: 'seq-', maxLength: 8192, maxActions: 16, highlightChannel: 'link', hoverChannel: 'link-hover' },
  revisions: { keepRecent: 64 },
  edit: { rebaseStale: false },
  preview: { dropped: 'issue' },
  bounds: { ...DefaultSpecBounds, maxBlockValues: DefaultMaxBlockValues },
  copyLink: { prefer: ['ref', 'ids'], numbering: 'ordinal', action: 'focus' },
  story: DefaultStoryConfig,
};

type DeepPartial<T> = {
  readonly [K in keyof T]?: T[K] extends readonly unknown[] ? T[K] : T[K] extends object ? DeepPartial<T[K]> : T[K];
};

const ZoomOffDefault = { ...DefaultRowsZoom, enabled: false };

/**
 * `rows.zoom` as the keys, the grips and the layout read it: finite and positive bounds, `minPx <= maxPx`, else the
 * defaults (a bad value never reaches a layout).
 */
export function validRowsZoom(z: ViewerConfig['rows']['zoom']): ViewerConfig['rows']['zoom'] {
  const ok = Number.isFinite(z.minPx) && Number.isFinite(z.maxPx) && z.minPx > 0 && z.maxPx > 0 && z.minPx <= z.maxPx;
  return ok ? z : z.enabled === false ? ZoomOffDefault : DefaultRowsZoom;
}

/** The defaults with the given items over them, one level of groups deep. */
export function viewerConfig(over: DeepPartial<ViewerConfig> = {}): ViewerConfig {
  const out: Record<string, unknown> = {};
  for (const [group, items] of Object.entries(DefaultViewerConfig))
    out[group] = { ...items, ...(over as Record<string, object | undefined>)[group] };
  // the view's keys merge by command id: an override of one command's keys keeps the others' (`[]` removes them)
  const keys = over.view?.keys;
  if (keys) (out.view as { keys: unknown }).keys = { ...DefaultViewerConfig.view.keys, ...keys };
  // the zoom strip's items merge: one set keeps the other's default
  const zoomStrip = over.navigator?.zoomStrip;
  if (zoomStrip)
    (out.navigator as { zoomStrip: unknown }).zoomStrip = { ...DefaultViewerConfig.navigator.zoomStrip, ...zoomStrip };
  // the rows' zoom items merge as the zoom strip's do
  const rowsZoom = over.rows?.zoom;
  if (rowsZoom) (out.rows as { zoom: unknown }).zoom = validRowsZoom({ ...DefaultViewerConfig.rows.zoom, ...rowsZoom });
  // `view.motion: 'auto'`, the earlier spelling of `platform` (before one item decided jumps and stories), reads as it
  const view = out.view as { motion: string };
  if (view.motion === 'auto') view.motion = 'platform';
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
