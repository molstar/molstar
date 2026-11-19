/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { HalfOpenRanges } from '@molstar/core/data/int/half-open-ranges';
import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { type Coverage, nodeKey, type Want } from '@molstar/seqx-core/data';
import type { Location, LocationPart } from '@molstar/seqx-core/query';
import { type OverlayNode, projectRanges, type StaticModel } from '@molstar/seqx-core/resolve';
import {
  axisBreaks,
  breakHolding,
  expandTemplate,
  followRange,
  followWindow,
  offsetOf,
  piecesTotal,
  templateNodes,
  type ResolvedLayer,
  type ResolvedSection,
  type ResolvedTrack,
  type ResolvedView,
  type TrackTemplate,
} from '@molstar/seqx-core/view';
import { AxisTransform, Cut, type Viewport } from './axis.js';
import { DefaultRenderConfig, fontSize, type RenderConfig } from './config.js';
import {
  type DrawLists,
  drawPass,
  NoColor,
  type PassName,
  Passes,
  type PassWriter,
  passDraws,
  Round,
} from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type {
  Box,
  Emphasis,
  EmphasisHighlight,
  Hit,
  LabelSink,
  LayerType,
  LayoutContext,
  Overlay,
  StoreReader,
} from './layer-type.js';
import { walkStripes } from './highlight-stripes.js';
import { boundedLevel } from './layer-data.js';
import { EmptyOverlay, sectionId } from './layer-type.js';
import type { LettersLay } from './letters.js';
import { breakHit } from './marks.js';
import type { ScaleParams } from '@molstar/seqx-core/theme/index';
import { DefaultThemeContext, type ThemeContext } from './channels.js';
import { elide, type TextMetrics } from './text.js';
import { CategoryCount, type SeqTokens } from './tokens.js';
import {
  BuiltInSizings,
  ContentSizing,
  heightsFit,
  LabelLane,
  type LayerContent,
  lanesHeight,
  type SizedTrack,
  type TrackSizing,
  type TrackSizingHook,
} from './sizing.js';
import { PlacePoint, StripPlacement, type Strips, WarpKind } from './layout/strips.js';
import { TranslateWarp } from './layout/warp.js';
import { chipText } from './layout/blocks.js';
import { TextLayout } from './layout/text.js';
import { TracksLayout } from './layout/tracks.js';
import { GridPlacement } from './layout/grid.js';
import type { PlaceContext, Placement, SectionLayoutType } from './layout/types.js';
import { drawRail, hitRail, placeRail, placeRowsRail, Rail, SlimRailPx } from './rail.js';
import { RowsLay } from './rows.js';
import { rectOver } from './grid.js';
import { columnsRoom, fitRowColumns, RowHeight, type RowColumn, rowColumnsOf } from './rows-model.js';
import type { ImageCache } from './images.js';
import { wellClipArgs, wellCorners } from './wells.js';

// The tracks layout (_spec/README.md "Render", Tracks layout): sections stacked, each with its axis transform from its
// system and range (a follower's from its leader's view), tracks with heights in px under one label column, each
// track's layers in order over its box. Tracks below the page are not laid out: the draw lists are bounded by the
// screen. The layout keeps each layer's last layout, so buffers and lanes carry over frame to frame.

/** No layer types, no layouts: what a layout is given when its input names none (it draws nothing and says so). */
const NoTypes: Readonly<Record<string, LayerType<any, any>>> = {};
const NoLayouts: Readonly<Record<string, SectionLayoutType<any>>> = {};

export interface LayerLayout {
  readonly key: string;
  readonly layer: ResolvedLayer;
  readonly type: LayerType<any, any>;
  readonly lay: unknown;
}

export interface TrackLayout {
  readonly track: ResolvedTrack;
  /** The data box: its band, inside its well by its padding (what its layers lay out, hit and emphasize in). */
  readonly box: Box;
  /** Its well: the data box grown by its padding (_spec/README.md "Render", Tracks layout), filled in the background. */
  readonly well: Box;
  /**
   * The label as drawn (cut to the label column with an ellipsis, never a hard clip).
   * @internal
   */
  readonly labelText: string;
  /** The label's measured width and its baseline (drawing never measures). */
  readonly labelWidth: number;
  readonly labelY: number;
  readonly layers: readonly LayerLayout[];
  /**
   * The "n filtered" chip: the rows a layer's filter left out (the first layer with some), its text and its box at the
   * label column's right end; a click on it clears the filter (`target: 'filter'`).
   */
  readonly filter?: {
    readonly layer: string;
    readonly count: number;
    readonly text: string;
    readonly x: number;
    readonly width: number;
    /** Its line below the label's (px), on a track tall enough for two lines: the label then whole; else 0, beside it. */
    readonly dy: number;
  };
}

export interface SectionLayout {
  readonly section: ResolvedSection;
  /** The section's id (its address's last step): what a hit names it by. */
  readonly id: string;
  /** Its layout (_spec/README.md "Render", Section layouts) and where it placed the section. */
  readonly layout: SectionLayoutType<any>;
  readonly placement: Placement;
  /** Undefined while the section's length is unknown (pending); else its transform over its window (a follower reads it). */
  readonly axis?: AxisTransform;
  /** Its tracks on the page and the gaps between them, the data's x: the column a `column` emphasis spans. */
  readonly band?: Box;
  readonly y: number;
  readonly height: number;
  /** Its tracks as bands: on the page in tracks; in text, those of its first line. */
  readonly tracks: readonly TrackLayout[];
  /** The breaks in view, `[x, hiddenStart, hiddenEnd, …]` (x unsnapped): a hairline each across the band. */
  readonly breaks: readonly number[];
  /** Its ruler rail, when it draws one (_spec/README.md "Render", Navigators). */
  readonly rail?: Rail;
  /** Its rows layers' rails (upright, at the right edge of each band), when they draw them. */
  readonly rowsRails?: readonly Rail[];
  /** The direction mark in the first track's label row (`axis.directionMark`): its text, x and measured width. */
  readonly mark?: { readonly text: string; readonly x: number; readonly width: number };
  /** The layers' lays a strip and band (`strip * bands + band`), in the strip's frame. */
  readonly lays: readonly (readonly LayerLayout[])[];
  /** Each strip's breaks in view (its frame), as `breaks`. */
  readonly stripBreaks: readonly (readonly number[])[];
  /** Its `section` layers' lays: laid out once a section, drawn in page coordinates. */
  readonly sectionLays: readonly LayerLayout[];
  /** Each band's box a strip (`strip * bands + band`, its frame) and each strip's band column (none when empty). */
  readonly boxes: readonly Box[];
  readonly columns: readonly (Box | undefined)[];
  /** What its chrome reads of the view. */
  readonly pageHeight: number;
  readonly tokens: SeqTokens;
  readonly config: RenderConfig;
}

/** How far from a break's notch a point hits it, CSS px (fixed, as other hit boxes are). */
export const BreakSlop = 3;
/** A break's notch above the band: its width and height, CSS px. */
const Notch = { width: 5, height: 3 };

/** A section's direction mark, or none: a reversed axis says so (`reversed`, or `along` what it reads along). */
function directionMark(section: ResolvedSection, config: RenderConfig): string | undefined {
  if (config.axis.directionMark === 'none' || section.axis.dir !== -1) return undefined;
  const along = section.axis.along;
  return along === undefined ? '← reversed' : `← along ${along.slice(along.lastIndexOf('/') + 1)}`;
}

export interface ViewLayout {
  readonly width: number;
  readonly height: number;
  readonly dpr: number;
  readonly tokens: SeqTokens;
  readonly config: RenderConfig;
  readonly sections: readonly SectionLayout[];
  /** What the layers need from the store now. */
  readonly wants: readonly Want[];
  /** The issues of this layout (a color it cannot read, a layer it cannot draw), each once. */
  readonly issues: readonly Issue[];
  /** The dynamic nodes this layout shows, with what they name (with `nodes.holdShown`): its owner holds them. */
  readonly shown: readonly OverlayNode[];
}

export interface LayoutInput {
  readonly model: StaticModel;
  readonly view: ResolvedView;
  readonly size: { readonly width: number; readonly height: number };
  /** The window of display units of every section that does not follow another; all of each by default. */
  readonly viewport?: Viewport;
  /** A section's own window (by its address), before `viewport`: a caller placing sections one by one. */
  readonly viewports?: Readonly<Record<string, Viewport>>;
  readonly dpr?: number;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  readonly config?: RenderConfig;
  readonly data: StoreReader;
  readonly overlay?: Overlay;
  readonly settled?: boolean;
  /** The layer types drawn, by name (the composition's registry): none when absent (every layer is said, not drawn). */
  readonly types?: Readonly<Record<string, LayerType<any, any>>>;
  /** The schemes and measures themes read (the built-ins by default) and the reader's color-blind substitution. */
  readonly theme?: Partial<Pick<ThemeContext, 'schemes' | 'measures' | 'cvdSafe'>>;
  /** The section layouts by name (the composition's registry): tracks for every section when absent. */
  readonly layouts?: Readonly<Record<string, SectionLayoutType<any>>>;
  /** Lays out every line of a text section, its box its content (an export: the figure). */
  readonly figure?: boolean;
  /** Draws the windows a figure leaves out by default (a ring's in an export). */
  readonly windows?: boolean;
  /** The clades a rows layer folds by its track's address, where its band gives none (an export's `folds`). */
  readonly folds?: Readonly<Record<string, readonly string[]>>;
  /**
   * The track sizings a section's `heights` may name besides the built-ins (`content`, `uniform`, `fixed`): the
   * composition's `track-sizing` registry entries, by name (_spec/README.md "Render", Track heights).
   */
  readonly sizings?: Readonly<Record<string, TrackSizing>>;
  /** The host's `trackSizing` port: the last word on each section's tracks' heights. */
  readonly trackSizing?: TrackSizingHook;
  /**
   * A navigator's picture: its tracks small, no label lane (its labels are the view's).
   * @internal
   */
  readonly picture?: boolean;
  /** The shared image cache (_spec/README.md "Render", The shared reduction): the renderer's. */
  readonly images?: ImageCache;
  /**
   * The sections that draw a ruler rail, and the rows tracks that draw a rows rail, by address (_spec/README.md
   * "Render", Navigators): part of each section's key by identity (a slim rail takes a row, a rows rail a column at the
   * right edge). Which can move is set on the rails each frame by the renderer, out of the key: idle wakes with no
   * relayout.
   */
  readonly rails?: Readonly<Record<string, unknown>>;
  /**
   * Another layout's held domains (_spec/README.md "Render", Navigators: a navigator's picture colors by its view's,
   * never by a domain of its own over the whole axis): looked up by `map(key)`; this layout then holds none.
   */
  readonly domains?: HeldDomains;
}

/** A layout's held domains as another reads them: its state, and a key here to a key there. */
export interface HeldDomains {
  readonly state: LayoutState;
  map(key: string): string;
}

/** What a layout keeps between frames: each layer's prepared form (per resolved layer) and its last layout. */
export class LayoutState {
  readonly preps = new WeakMap<ResolvedLayer, unknown>();
  readonly lays = new Map<string, { type: LayerType<any, any>; lay: unknown }>();
  /**
   * Open scales' domains (_spec/README.md "Themes", Palettes and scales): the union over the layers naming each, of their
   * values in view, computed when the view settles and held while it moves.
   */
  domains = new Map<string, readonly [number, number]>();
  /** Bumped when `domains` changes: part of a section's reuse. */
  domainsVersion = 0;
  /** Each section's last layout and what it was laid out for (_spec/README.md "Render", The key: per section). */
  readonly sections = new Map<string, SectionCache>();
  /** The row columns' fit, kept while what it read is the same. */
  columnsFit: unknown;
  /** The template instances shown last (their tracks' addresses): kept while a zoom's new level is pending. */
  readonly instancesShown = new Set<string>();
  /**
   * Track labels cut to the label column, by label and width (the metrics and font are the state's renderer's), for
   * the metrics' version (`elidedVersion`: a web font arriving cuts them again).
   */
  readonly elided = new Map<string, Map<number, string>>();
  elidedVersion = -1;
  /** Each section's tracks' heights and contents while their key is (`bandHeights`), and as a leader's reading. */
  readonly sizes = new Map<string, SizedHeights>();
  readonly leadSizes = new Map<string, SizedHeights>();
  /** The bands' label lanes, by band index (one a band, reset before its layers are laid out). */
  readonly labelLanes: LabelLane[] = [];
}

/** A section's tracks' heights with what they were sized for (`layoutView`'s `bandHeights`). */
interface SizedHeights {
  readonly tracks: readonly ResolvedTrack[];
  readonly fit: number;
  readonly config: RenderConfig;
  readonly layout: SectionLayoutType<any>;
  readonly version: number;
  readonly sizings: LayoutInput['sizings'];
  readonly hook: LayoutInput['trackSizing'];
  readonly heights: readonly number[];
  /** Each track's content (one lane block a track), undefined for a track not sized by its content. */
  readonly contents: readonly (LayerContent | undefined)[];
  readonly issue?: Issue;
}

/** Two lists of the same tracks (by identity). */
function sameTracks(a: readonly ResolvedTrack[], b: readonly ResolvedTrack[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/** A layer's key: its address, else its track's and its index. */
export function layerKey(track: ResolvedTrack, layer: ResolvedLayer): string {
  return layer.address ?? `${track.address}:${layer.index}`;
}

/**
 * A section's own pieces: its axis's (the array itself: lookups are kept per pieces array), or the whole system once its
 * source gave the length; undefined while unknown.
 */
function ownPieces(section: ResolvedSection, model: StaticModel, data: StoreReader): readonly number[] | undefined {
  if (section.axis.ranges) return section.axis.ranges;
  const system = model.systems.get(section.axis.system);
  if (!system?.source) return undefined;
  const info = data
    .coverage(system.source)
    ?.systems.find((s) => s.name === system.name || s.aliases?.includes(system.name));
  return info?.length === undefined ? undefined : placedPieces(info.length);
}

/** The whole system's one piece for a length, kept per length (lookups are kept per pieces array). */
const placed = new Map<number, readonly number[]>();
function placedPieces(length: number): readonly number[] {
  let p = placed.get(length);
  if (!p) {
    if (placed.size > 64) placed.clear();
    placed.set(length, (p = [0, length]));
  }
  return p;
}

/**
 * A follower's pieces and window from its leader's view (_spec/README.md "Render", Tracks layout): the images of the
 * leader's window edges as offsets along the follower's own pieces (`followWindow`, as the viewer's follow), so a
 * follower moves with its leader and keeps its direction and breaks. A window over a cyclic follower's origin is drawn
 * over the arc's pieces (`followRange`).
 */
function followAxis(
  section: ResolvedSection,
  own: readonly number[] | undefined,
  leader: AxisTransform,
  window: Viewport,
  config: RenderConfig,
) {
  const axis = { ...section.axis, ...(own ? { ranges: own } : {}) };
  // a cyclic leader's window past its system's end (a linear view across its origin, its pieces twice) wraps
  const turn = leader.cyclic && leader.length !== undefined ? Math.min(leader.total, leader.length) : leader.total;
  if (own && window.end <= turn) {
    const w = followWindow(section, leader, axis, window.start, window.end, config.followArc);
    if (w && !w.wraps) return { pieces: own, viewport: { start: w.start, end: w.end } };
  }
  const pieces = followRange(section, windowPositions(leader, window), { arc: config.followArc });
  return pieces ? { pieces } : undefined;
}

/**
 * The positions a window of a leader's display units shows, `[start, end, …]`: the transform's own window as it is; a
 * cyclic window past the total wraps (two runs, the second from the origin).
 */
function windowPositions(leader: AxisTransform, w: Viewport): number[] {
  if (w.start === leader.u0 && w.end === leader.u1) return leader.visible([]);
  const out: number[] = [];
  const at = (start: number, end: number) =>
    new AxisTransform(leader.system, leader.pieces, 0, 1, 1, { start, end }, { dir: leader.dir }).visible(out);
  const total = leader.total;
  if (leader.cyclic && w.end > total) {
    at(w.start, total);
    at(0, Math.min(w.end - total, w.start));
  } else at(w.start, w.end);
  return out;
}

const Twice = new WeakMap<readonly number[], readonly number[]>();

/**
 * A cyclic axis's pieces under a window past its total (_spec/README.md "Render", Tracks layout: a linear view panned
 * across its origin): the whole system twice, so the window's units go on over the origin as the second turn's (two
 * runs, no break between); else the pieces as they are. Kept per pieces array (a pan makes none).
 */
function turned(section: ResolvedSection, own: readonly number[], viewport: Viewport | undefined): readonly number[] {
  const { cyclic, length } = section.axis;
  if (!cyclic || !viewport || own.length !== 2 || own[0] !== 0 || own[1] !== length || viewport.end <= own[1])
    return own;
  let t = Twice.get(own);
  if (!t) Twice.set(own, (t = [own[0], own[1], own[0], own[1]]));
  return t;
}

/**
 * A following ring's window on its own axis (_spec/README.md "Render", Rings layout: a ring that follows a linear
 * section has the leader's range as its window): the follow's window as it is, else the first piece of the image its
 * offsets on, as long as the image (over the origin: past the total).
 */
function followedWindow(
  section: ResolvedSection,
  own: readonly number[],
  made: { pieces: readonly number[]; viewport?: Viewport },
): Viewport {
  if (made.viewport) return made.viewport;
  const shape = { ...section.axis, ranges: own };
  const p = made.pieces;
  const len = piecesTotal(p);
  // the image's first piece in display order: its low edge forward, its high edge reversed
  const start = offsetOf(shape, section.axis.dir === -1 ? p[1] : p[0]) ?? 0;
  return { start, end: start + len };
}

/** A section's layout: its own `layout`, else its view's. */
export function sectionLayoutName(section: ResolvedSection, view: ResolvedView): string {
  return section.layout ?? view.layout;
}

/** Whether a section's blocks are on more than one system (_spec/README.md "Render", Blocks: it neither follows nor leads). */
export function multiSystem(section: ResolvedSection): boolean {
  const b = section.blocks;
  if (!b) return false;
  for (let i = 1; i < b.length; i++) if (b[i].axis.system !== b[0].axis.system) return true;
  return false;
}

const NoLays: readonly LayerLayout[] = [];
const NoSectionLays: LayerLayout[] = [];

/** A layer's type: the first of its fallbacks drawn here (a loop: no array made a layer a frame). */
function firstType(
  layer: ResolvedLayer,
  types: Readonly<Record<string, LayerType<any, any>>>,
): LayerType<any, any> | undefined {
  const f = layer.fallbacks;
  for (let i = 0; i < f.length; i++) {
    const t = types[f[i]];
    if (t !== undefined) return t;
  }
  return undefined;
}

/**
 * A section's reads of the store (_spec/README.md "Render", The key: per section): the entries and coverages its
 * layout read with their versions, so a section whose reads hold is not laid out again. Its arrays are kept from layout
 * to layout (a count, not a new array: a layout allocates a bounded amount).
 */
class SectionReads implements StoreReader {
  private keys: string[] = [];
  private versions: number[] = [];
  private n = 0;
  private sources: string[] = [];
  private coverages: (Coverage | undefined)[] = [];
  private m = 0;
  private data!: StoreReader;
  begin(data: StoreReader): void {
    this.data = data;
    this.n = 0;
    this.m = 0;
  }
  coverage(source: string) {
    const c = this.data.coverage(source);
    this.sources[this.m] = source;
    this.coverages[this.m++] = c;
    return c;
  }
  keysOf(w: Want) {
    return this.data.keysOf(w);
  }
  tooMany(w: Want) {
    return this.data.tooMany?.(w) ?? false;
  }
  get settings() {
    return this.data.settings;
  }
  get(key: string) {
    const e = this.data.get(key);
    this.keys[this.n] = key;
    this.versions[this.n++] = e?.version ?? -1;
    return e;
  }
  /** Whether they read the same again (read through `data`, so a renderer's reader records them as read). */
  hold(data: StoreReader): boolean {
    for (let i = 0; i < this.n; i++) if ((data.get(this.keys[i])?.version ?? -1) !== this.versions[i]) return false;
    for (let i = 0; i < this.m; i++) if (data.coverage(this.sources[i]) !== this.coverages[i]) return false;
    return true;
  }
}

/** Where a section's layout reports to this frame: the view's wants, issues, held domains, lays used, nodes shown. */
interface FrameSinks {
  readonly wants: Want[];
  readonly report: (i: Issue) => void;
  readonly theme: ThemeContext;
  readonly used: Set<string>;
  readonly showing: Set<string>;
}

/**
 * What a section kept from its last layout, what it was laid out for, and what that layout gave the frame (wants,
 * issues, extents, lays used, nodes shown), replayed when it is reused. One a section, kept across frames: its arrays,
 * its reads, its theme and its sinks are made once (a layout allocates a bounded amount).
 */
class SectionCache {
  layout!: SectionLayout;
  section!: ResolvedSection;
  stamp!: ViewStamp;
  y = 0;
  axis: AxisTransform | undefined;
  gridY: PlaceContext['y'];
  viewport: Viewport | undefined;
  readonly reads = new SectionReads();
  private wants: Want[] = [];
  private nw = 0;
  private issues: Issue[] = [];
  private ni = 0;
  private used: string[] = [];
  private nu = 0;
  private showing: string[] = [];
  private ns = 0;
  private extentKeys: string[] = [];
  private extentValues: number[] = [];
  private ne = 0;
  private sinks!: FrameSinks;

  /** Starts a layout of the section for a frame. */
  begin(data: StoreReader, sinks: FrameSinks): void {
    this.reads.begin(data);
    this.sinks = sinks;
    this.nw = 0;
    this.ni = 0;
    this.nu = 0;
    this.ns = 0;
    this.ne = 0;
  }

  readonly want = (w: Want): void => {
    this.wants[this.nw++] = w;
    this.sinks.wants.push(w);
  };

  readonly report = (i: Issue): void => {
    this.issues[this.ni++] = i;
    this.sinks.report(i);
  };

  addUsed(key: string): void {
    this.used[this.nu++] = key;
    this.sinks.used.add(key);
  }

  addShowing(a: string, b: string): void {
    this.showing[this.ns++] = a;
    this.showing[this.ns++] = b;
    this.sinks.showing.add(a).add(b);
  }

  /** The themes a section's layers read: the frame's, with the extents it saw recorded. */
  readonly theme: ThemeContext = (() => {
    const self = this;
    return {
      get schemes() {
        return self.sinks.theme.schemes;
      },
      get measures() {
        return self.sinks.theme.measures;
      },
      get cvdSafe() {
        return self.sinks.theme.cvdSafe;
      },
      scale: (a) => self.sinks.theme.scale(a),
      domain: (k) => self.sinks.theme.domain(k),
      extent: (k, lo, hi) => {
        self.extentKeys[self.ne] = k;
        self.extentValues[2 * self.ne] = lo;
        self.extentValues[2 * self.ne + 1] = hi;
        self.ne++;
        self.sinks.theme.extent(k, lo, hi);
      },
    } as ThemeContext;
  })();

  /** Gives this frame what the kept layout gave its own. */
  replay(sinks: FrameSinks): void {
    for (let i = 0; i < this.nw; i++) sinks.wants.push(this.wants[i]);
    for (let i = 0; i < this.ni; i++) sinks.report(this.issues[i]);
    for (let i = 0; i < this.nu; i++) sinks.used.add(this.used[i]);
    for (let i = 0; i < this.ns; i++) sinks.showing.add(this.showing[i]);
    for (let i = 0; i < this.ne; i++)
      sinks.theme.extent(this.extentKeys[i], this.extentValues[2 * i], this.extentValues[2 * i + 1]);
  }
}

/** The view's parts a section's layout reads (size, ratio, tokens, text, overlay, config, types, themes, domains). */
interface ViewStamp {
  readonly model: string;
  readonly width: number;
  readonly height: number;
  readonly dpr: number;
  readonly tokens: number;
  readonly text: number;
  readonly overlay: number;
  readonly settled: boolean;
  readonly config: RenderConfig;
  readonly types: unknown;
  readonly schemes: unknown;
  readonly measures: unknown;
  readonly cvdSafe: boolean;
  readonly domains: number;
  readonly figure: boolean;
  readonly windows: boolean;
  readonly rails: unknown;
}

function viewStamp(
  input: LayoutInput,
  config: RenderConfig,
  dpr: number,
  overlay: Overlay,
  state: LayoutState,
): ViewStamp {
  return {
    model: input.model.hash,
    width: input.size.width,
    height: input.size.height,
    dpr,
    tokens: input.tokens.version,
    text: input.text.version,
    overlay: overlay.version,
    settled: input.settled ?? false,
    config,
    types: input.types,
    schemes: input.theme?.schemes,
    measures: input.theme?.measures,
    cvdSafe: input.theme?.cvdSafe ?? false,
    domains: input.domains ? input.domains.state.domainsVersion : state.domainsVersion,
    figure: input.figure ?? false,
    windows: input.windows ?? false,
    rails: input.rails,
  };
}

function sameStamp(a: ViewStamp, b: ViewStamp): boolean {
  return (
    a.model === b.model &&
    a.width === b.width &&
    a.height === b.height &&
    a.dpr === b.dpr &&
    a.tokens === b.tokens &&
    a.text === b.text &&
    a.overlay === b.overlay &&
    a.settled === b.settled &&
    a.config === b.config &&
    a.types === b.types &&
    a.schemes === b.schemes &&
    a.measures === b.measures &&
    a.cvdSafe === b.cvdSafe &&
    a.domains === b.domains &&
    a.figure === b.figure &&
    a.windows === b.windows &&
    a.rails === b.rails
  );
}

/** Two transforms over the same window of the same pieces, at the same place. */
function sameAxis(a: AxisTransform | undefined, b: AxisTransform | undefined): boolean {
  if (a === b) return true;
  if (!a || !b) return false;
  return (
    a.u0 === b.u0 &&
    a.u1 === b.u1 &&
    a.x0 === b.x0 &&
    a.width === b.width &&
    a.dir === b.dir &&
    samePieces(a.pieces, b.pieces)
  );
}

/** Two pieces arrays alike by value (a follower's are made again each frame). */
function samePieces(a: readonly number[], b: readonly number[]): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}

/** Whether a section's kept layout is this frame's: its view parts, place, window, viewport by its layout, reads. */
/** Two grids' y place alike: the same pieces and window. */
function sameGridY(a: PlaceContext['y'], b: PlaceContext['y']): boolean {
  if (a === b) return true;
  if (!a || !b || !samePieces(a.pieces, b.pieces) || a.axis.dir !== b.axis.dir) return false;
  const va = a.viewport;
  const vb = b.viewport;
  return va === vb || (!!va && !!vb && va.start === vb.start && va.end === vb.end);
}

function reusable(
  kept: SectionCache,
  section: ResolvedSection,
  layout: SectionLayoutType<any>,
  stamp: ViewStamp,
  y: number,
  axis: AxisTransform | undefined,
  viewport: Viewport | undefined,
  data: StoreReader,
): boolean {
  return (
    kept.section === section &&
    kept.layout.layout === layout &&
    kept.y === y &&
    axis !== undefined &&
    sameStamp(kept.stamp, stamp) &&
    sameAxis(kept.axis, axis) &&
    layout.sameKey(kept.viewport, viewport) &&
    kept.reads.hold(data)
  );
}

/** The union of a strip's bands in its frame (`ctx.span`: what a layer with `across` outlines). */
function stripSpan(strips: Strips, i: number): Box {
  const a = strips.axes[i];
  let y0 = Infinity;
  let y1 = -Infinity;
  for (let j = 0; j < strips.bands; j++) {
    const k = i * strips.bands + j;
    if (strips.bandH[k] <= 0) continue;
    if (strips.bandY[k] < y0) y0 = strips.bandY[k];
    if (strips.bandY[k] + strips.bandH[k] > y1) y1 = strips.bandY[k] + strips.bandH[k];
  }
  return y1 > y0 ? { x: a.x0, y: y0, width: a.width, height: y1 - y0 } : { x: a.x0, y: 0, width: a.width, height: 0 };
}

/**
 * The bands of the tracks a layer's `across` names on a strip, the first to the last and those between (`ctx.span`: a
 * frame's outline), in its frame; the whole strip's when it names none shown.
 */
function acrossSpan(strips: Strips, i: number, shown: readonly ResolvedTrack[], across: readonly unknown[]): Box {
  let a = -1;
  let b = -1;
  for (let j = 0; j < shown.length; j++) {
    const id = shown[j].address.slice(shown[j].address.lastIndexOf('/') + 1);
    if (id !== across[0] && id !== across[across.length - 1]) continue;
    if (a < 0) a = j;
    b = j;
  }
  if (a < 0) return stripSpan(strips, i);
  const x = strips.axes[i];
  let y0 = Infinity;
  let y1 = -Infinity;
  for (let j = a; j <= b; j++) {
    const k = i * strips.bands + j;
    if (strips.bandH[k] <= 0) continue;
    y0 = Math.min(y0, strips.bandY[k]);
    y1 = Math.max(y1, strips.bandY[k] + strips.bandH[k]);
  }
  return y1 > y0 ? { x: x.x0, y: y0, width: x.width, height: y1 - y0 } : stripSpan(strips, i);
}

/** A strip's band column in its frame: its bands above the page's bottom and the gaps between them; none when empty. */
function stripColumn(strips: Strips, i: number, pageHeight: number): Box | undefined {
  const a = strips.axes[i];
  let y0 = Infinity;
  let y1 = -Infinity;
  for (let j = 0; j < strips.bands; j++) {
    const k = i * strips.bands + j;
    if (strips.bandH[k] <= 0 || strips.bandY[k] + strips.dy[i] >= pageHeight) continue;
    if (strips.bandY[k] < y0) y0 = strips.bandY[k];
    if (strips.bandY[k] + strips.bandH[k] > y1) y1 = strips.bandY[k] + strips.bandH[k];
  }
  return y1 > y0 ? { x: a.x0, y: y0, width: a.width, height: y1 - y0 } : undefined;
}

/** A pending section (no axis yet): no strips, its tracks' height kept. */
function pendingPlacement(
  top: number,
  box: Box,
  heights: readonly number[],
  pads: readonly number[],
  config: RenderConfig,
): Placement {
  const p = new StripPlacement();
  let h = 0;
  for (let j = 0; j < heights.length; j++) h += (j ? config.trackGap : 0) + pads[2 * j] + heights[j] + pads[2 * j + 1];
  p.strips.reset(heights.length);
  p.contentHeight = h;
  p.box = { x: config.inset, y: top, width: box.width, height: h };
  return p;
}

/**
 * The rows' header in the view's one label column (_spec/README.md "Render", The rows layer): a rows layer's tree in
 * tree order and its first residue numbers widen the label column by their width (the labels keep theirs), its last
 * numbers take a right gutter. The widest of a view's tracks sections, so every section's data starts at one x. Its
 * row columns (Row columns, in a text figure's too) widen the label column by the left ones and the right gutter by
 * the right ones, by layer.
 */
function rowsHeader(
  view: ResolvedView,
  rails: LayoutInput['rails'],
): {
  readonly tree: number;
  readonly numbers: boolean;
  readonly rail: boolean;
} {
  let tree = 0;
  let numbers = false;
  let rail = false;
  for (const s of view.sections) {
    if (sectionLayoutName(s, view) !== 'tracks') continue;
    for (const t of s.tracks)
      for (const l of t.layers) {
        if (l.type !== 'rows' || !l.placed || l.params.hidden === true) continue;
        if (l.params.numbers !== 'none') numbers = true;
        if (rails?.[t.address] !== undefined) rail = true;
        const w = l.params.tree_width;
        if (l.params.order === 'tree' && l.params.tree !== undefined)
          tree = Math.max(tree, typeof w === 'number' ? w : 120);
      }
  }
  return { tree, numbers, rail };
}

const NoColumns: ReadonlyMap<string, readonly RowColumn[]> = new Map();
const ViewColumns = new WeakMap<ResolvedView, ReadonlyMap<string, readonly RowColumn[]>>();

/** The row columns of a view's rows layers (tracks and text sections), by layer: read once a view. */
function viewColumns(view: ResolvedView): ReadonlyMap<string, readonly RowColumn[]> {
  let out = ViewColumns.get(view);
  if (out) return out;
  let columns: Map<string, readonly RowColumn[]> | undefined;
  for (const s of view.sections) {
    const name = sectionLayoutName(s, view);
    if (name !== 'tracks' && name !== 'text') continue;
    for (const t of s.tracks)
      for (const l of t.layers) {
        if (l.type !== 'rows' || !l.placed || l.params.hidden === true) continue;
        const cols = l.params.columns === undefined ? undefined : rowColumnsOf(l.params.columns);
        if (cols?.length) (columns ??= new Map()).set(layerKey(t, l), cols);
      }
  }
  out = columns ?? NoColumns;
  ViewColumns.set(view, out);
  return out;
}

/** The title row a rows track's layer with columns takes above its rows (_spec/README.md "Render", Row columns). */
function titleRowOf(track: ResolvedTrack, config: RenderConfig): number {
  for (const l of track.layers)
    if (l.type === 'rows' && Array.isArray(l.params.columns) && l.params.columns.length)
      return config.rows.columns.titlePx;
  return 0;
}

/**
 * The row columns' widths a view grants by layer (_spec/README.md "Render", Row columns: EC-107): what the view lacks
 * for the data's least width at every column's own width is taken from the columns, right ones first, each down to a
 * floor and then hidden; and the room the widest layer's left and right columns take.
 */
function columnsFit(
  state: LayoutState,
  columns: ReadonlyMap<string, readonly RowColumn[]>,
  width: number,
  fixed: number,
  config: RenderConfig,
  rail: boolean,
): ColumnsFit | undefined {
  if (!columns.size) return undefined;
  const kept = state.columnsFit as ColumnsFit | undefined;
  if (kept && kept.columns === columns && kept.width === width && kept.fixed === fixed && kept.config === config)
    return kept;
  const { floorPx, gapPx, dataMinPx } = config.rows.columns;
  const widths = new Map<string, number[]>();
  let left = 0;
  let right = 0;
  for (const [key, cols] of columns) {
    const own = cols.map((c) => c.width);
    const need = fixed + dataMinPx + columnsRoom(cols, own, 'left', gapPx) + columnsRoom(cols, own, 'right', gapPx);
    const w = fitRowColumns(cols, need - width, floorPx, gapPx, []);
    widths.set(key, w);
    left = Math.max(left, columnsRoom(cols, w, 'left', gapPx));
    right = Math.max(right, columnsRoom(cols, w, 'right', gapPx));
  }
  // a rows rail keeps the gutter's gap from the last right column
  if (rail && right) right += gapPx;
  const rightX = width - config.inset - (rail ? SlimRailPx : 0) - right;
  const grants = new Map<string, { widths: readonly number[]; rightX: number }>();
  for (const [key, w] of widths) grants.set(key, { widths: w, rightX });
  const fit: ColumnsFit = { columns, width, fixed, config, left, right, grants };
  state.columnsFit = fit;
  return fit;
}

/** The row columns' fit a view's layout keeps (what it was made for, the rooms, each layer's grant). */
interface ColumnsFit {
  readonly columns: ReadonlyMap<string, readonly RowColumn[]>;
  readonly width: number;
  readonly fixed: number;
  readonly config: RenderConfig;
  readonly left: number;
  readonly right: number;
  readonly grants: ReadonlyMap<string, { readonly widths: readonly number[]; readonly rightX: number }>;
}

/** A config with its label column widened by `extra` px, one a base and width. */
const Widened = new WeakMap<RenderConfig, Map<number, RenderConfig>>();
function widened(base: RenderConfig, extra: number): RenderConfig {
  if (!extra) return base;
  let m = Widened.get(base);
  if (!m) Widened.set(base, (m = new Map()));
  let c = m.get(extra);
  if (!c) m.set(extra, (c = { ...base, labelWidth: base.labelWidth + extra }));
  return c;
}

/**
 * What a layer with no color is colored by (_spec/README.md "Themes", A layer with no color): its key, the data it
 * shows (a ring and its detail showing one node: one color) or, for a template's layer, the template's (its instances:
 * one color); undefined for a layer the spec colors or a type that has its own default.
 */
function colorKey(track: ResolvedTrack, layer: ResolvedLayer): string | undefined {
  if (layer.type !== 'ranges' || layer.theme.color !== undefined) return undefined;
  return track.template ? `${track.address}:${layer.index}` : layer.data?.address;
}

/**
 * The colourless layers' places in the categorical palette, by their keys' first appearance in the view (a layer's
 * order): made once a view.
 */
const DefaultColors = new WeakMap<ResolvedView, Map<string, number>>();
function defaultColors(view: ResolvedView): Map<string, number> {
  let m = DefaultColors.get(view);
  if (m) return m;
  m = new Map();
  for (const s of view.sections)
    for (const t of s.tracks)
      for (const l of t.layers) {
        const k = colorKey(t, l);
        if (k !== undefined && !m.has(k)) m.set(k, m.size % CategoryCount);
      }
  DefaultColors.set(view, m);
  return m;
}

/**
 * The label column's width a view lays out with (_spec/README.md "Render", Configuration: `labelColumn`): the widest
 * track label, between the column's least and most; `labelWidth` when the
 * column is not measured. What the viewer's data box starts after, before any layout (its subpath, off the barrel).
 */
export function labelColumnWidth(view: ResolvedView, config: RenderConfig, text: TextMetrics): number {
  const c = config.labelColumn;
  if (!c.measure) return config.labelWidth;
  const kept = Columns.get(view);
  if (kept && kept.config === config && kept.text === text && kept.version === text.version) return kept.width;
  const width = measureColumn(view, config, text);
  Columns.set(view, { config, text, version: text.version, width });
  return width;
}

/** The label column's width by view, with the config and metrics it was measured for. */
const Columns = new WeakMap<
  ResolvedView,
  { readonly config: RenderConfig; readonly text: TextMetrics; readonly version: number; readonly width: number }
>();

function measureColumn(view: ResolvedView, config: RenderConfig, text: TextMetrics): number {
  const c = config.labelColumn;
  let w = 0;
  // the labels alone: a direction mark (a section reversed) fades at the column's end, so a reverse never moves the
  // data box (its windows keep their px a position)
  for (const s of view.sections)
    for (const t of s.tracks) {
      // an alignment's rows label themselves in the column: their names, not the tracks', are what it holds (the
      // column keeps its fixed width, its row labels their own ellipsis)
      for (const l of t.layers) if (l.type === 'rows') return config.labelWidth;
      const m = text.measure(config.labelFont, t.label);
      if (m > w) w = m;
    }
  return Math.min(c.max, Math.max(c.min, Math.ceil(w)));
}

/** A config with its label column `width` px wide, one a base and width. */
const Measured = new WeakMap<RenderConfig, Map<number, RenderConfig>>();
function measured(base: RenderConfig, width: number): RenderConfig {
  if (width === base.labelWidth) return base;
  let m = Measured.get(base);
  if (!m) Measured.set(base, (m = new Map()));
  let c = m.get(width);
  if (!c) m.set(width, (c = { ...base, labelWidth: width }));
  return c;
}

/** A track label cut to `max` px with an ellipsis, kept by label and width (a layout makes no string a track). */
function elidedLabel(state: LayoutState, text: TextMetrics, font: string, label: string, max: number): string {
  const kept = state.elided;
  if (state.elidedVersion !== text.version) {
    kept.clear();
    state.elidedVersion = text.version;
  }
  let m = kept.get(label);
  if (!m) {
    if (kept.size > 4096) kept.clear();
    kept.set(label, (m = new Map()));
  }
  let t = m.get(max);
  if (t === undefined) m.set(max, (t = elide(text, font, label, max)));
  return t;
}

/** The right gutter's width when a rows layer draws residue numbers: five digits and a gap. */
function numbersGutter(text: TextMetrics, config: RenderConfig): number {
  return text.measure(config.itemFont, '00000') + config.labelGap;
}

export function layoutView(input: LayoutInput, state = new LayoutState(), again = false): ViewLayout {
  const header = rowsHeader(input.view, input.rails);
  const given = input.config ?? DefaultRenderConfig;
  const base = measured(given, labelColumnWidth(input.view, given, input.text));
  // the first numbers take their own room in the label column, the labels keep theirs
  const first = header.numbers ? Math.ceil(input.text.measure(base.itemFont, '00000')) + 4 : 0;
  // the last numbers' gutter, and beside it a rows rail's column (fixed: a rail idle or live keeps the data's width)
  const gutter = (header.numbers ? numbersGutter(input.text, base) : 0) + (header.rail ? SlimRailPx : 0);
  // the row columns: the left ones widen the label column, the right ones the right gutter (and a gap before a rows
  // rail), yielding to a narrow view; kept while the view, the width and the room they read are
  const fixed = 2 * base.inset + base.labelWidth + base.labelGap + header.tree + first + gutter;
  const fit = columnsFit(state, viewColumns(input.view), input.size.width, fixed, base, header.rail);
  const config = widened(base, header.tree + first + (fit?.left ?? 0));
  const right = gutter + (fit?.right ?? 0);
  const grants = fit?.grants;
  // a text section's lines leave the right columns their room
  const textWidth = input.size.width - (fit?.right ?? 0);
  const extents = new Map<string, [number, number]>();
  const theme: ThemeContext = {
    schemes: input.theme?.schemes ?? DefaultThemeContext.schemes,
    ...(input.theme?.measures ? { measures: input.theme.measures } : {}),
    cvdSafe: input.theme?.cvdSafe ?? false,
    scale: (a) => input.view.scales?.[a] as ScaleParams | undefined,
    domain: input.domains
      ? (key) => input.domains!.state.domains.get(input.domains!.map(key))
      : (key) => state.domains.get(key),
    extent: (key, lo, hi) => {
      const e = extents.get(key);
      if (!e) extents.set(key, [lo, hi]);
      else {
        if (lo < e[0]) e[0] = lo;
        if (hi > e[1]) e[1] = hi;
      }
    },
  };
  const dpr = input.dpr ?? 1;
  const types = input.types ?? NoTypes;
  const overlay = input.overlay ?? EmptyOverlay;
  const { width, height } = input.size;
  const wants: Want[] = [];
  const issues: Issue[] = [];
  const seen = new Set<string>();
  const report = (i: Issue) => {
    const k = `${i.code}|${i.address ?? ''}`;
    if (seen.has(k)) return;
    seen.add(k);
    issues.push(i);
  };
  const dataX = config.inset + config.labelWidth + config.labelGap;
  const dataW = Math.max(1, width - dataX - config.inset - right);
  // axes first, leaders before their followers (a follower may come first in the view)
  const bySection = new Map(input.view.sections.map((s) => [s.address, s]));
  const layouts = input.layouts ?? NoLayouts;
  const sizings = BuiltInSizings;
  /** A section's layout: its name's in the layouts given, else tracks (a layout not given is said, once). */
  const layoutOf = (section: ResolvedSection) => {
    const name = sectionLayoutName(section, input.view);
    const layout = layouts[name];
    if (layout) return layout;
    if (name !== TracksLayout.name)
      report(
        issueOf(
          RenderIssues,
          'seq.render.layout',
          `The layout ${name} is not drawn here: the section is laid out as tracks`,
          {
            address: section.address,
            hint: RenderIssues.codes['seq.render.layout'].hint,
          },
        ),
      );
    return TracksLayout;
  };
  /**
   * What a follower reads of its leader (_spec/README.md "Render", Section layouts): the leader's viewport as its layout
   * reads it (`followWindow`: a ring's window, a text's lines in view), else its transform's window (tracks).
   */
  const leaderWindow = (section: ResolvedSection, axis: AxisTransform): Viewport => {
    const layout = layoutOf(section);
    if (!layout.followWindow) return { start: axis.u0, end: axis.u1 };
    const viewport = input.viewports?.[section.address] ?? input.viewport;
    const heights = bandHeights(section, section.tracks, layout, true, report).heights;
    const pads = section.tracks.flatMap((t) => {
      const p = t.padding ?? config.trackPadding;
      return [p.top, p.bottom];
    });
    return layout.followWindow(section, axis, viewport, {
      width: layout === TextLayout ? textWidth : width,
      config,
      text: input.text,
      heights,
      pads,
    });
  };
  /**
   * A section's tracks' heights (_spec/README.md "Render", Track heights): the section's `heights` sizing over its
   * tracks' content heights (a layout with heights of its own, text, keeps them), a track's own height over it, the
   * host's port over all. Kept in the state by the section (and `lead`, its leader's reading) while its key is: its
   * shown tracks, the px a position of its whole extent, the config, the layout, the metrics' version, the registry and
   * the port; so a pan runs neither the sizing nor the port. Each track's content (one lane block a track) is kept with
   * it.
   */
  const bandHeights = (
    section: ResolvedSection,
    tracks: readonly ResolvedTrack[],
    layout: SectionLayoutType<any>,
    lead: boolean,
    say: (i: Issue) => void,
  ): SizedHeights => {
    // the px a position with the section's own whole extent across the data box, whatever its window shows (a pan over
    // a cyclic origin, a follower's range, a zoom): a label fitting there fits at every zoom in
    const own = layout.heights ? undefined : ownPieces(section, input.model, input.data);
    let total = 0;
    if (own) for (let k = 0; k + 1 < own.length; k += 2) total += Math.abs(own[k + 1] - own[k]);
    const fit = total > 0 ? dataW / total : 0;
    const memo = lead ? state.leadSizes : state.sizes;
    const kept = memo.get(section.address);
    if (
      kept &&
      kept.fit === fit &&
      kept.config === config &&
      kept.layout === layout &&
      kept.version === input.text.version &&
      kept.sizings === input.sizings &&
      kept.hook === input.trackSizing &&
      sameTracks(kept.tracks, tracks)
    ) {
      if (kept.issue) say(kept.issue);
      return kept;
    }
    const made = sizeTracks(section, tracks, layout, fit);
    memo.set(section.address, made);
    if (made.issue) say(made.issue);
    return made;
  };
  /** The sizing pass itself (`bandHeights` keeps what it gives). */
  const sizeTracks = (
    section: ResolvedSection,
    tracks: readonly ResolvedTrack[],
    layout: SectionLayoutType<any>,
    fit: number,
  ): SizedHeights => {
    const byType = layout.heights?.(config) ?? config.heights;
    const fixed = tracks.map((t) => trackHeight(t, config, byType));
    const base = {
      tracks: [...tracks],
      fit,
      config,
      layout,
      version: input.text.version,
      sizings: input.sizings,
      hook: input.trackSizing,
    };
    if (layout.heights || !tracks.length) return { ...base, heights: fixed, contents: [] };
    const cctx = { config, text: input.text, fit, labelLane: layout === TracksLayout && !input.picture };
    const sized: SizedTrack[] = [];
    const contents: (LayerContent | undefined)[] = [];
    for (let j = 0; j < tracks.length; j++) {
      const t = tracks[j];
      let lanes: number | undefined;
      let labels = false;
      let other = 0;
      let kind = '';
      for (const l of t.layers) {
        // a hidden layer is no part of the track's kind or height
        if (l.params.hidden === true) continue;
        kind = kind ? `${kind}+${l.type}` : l.type;
        const type = l.placed ? firstType(l, types) : undefined;
        let c: LayerContent | undefined;
        if (type?.content) {
          let prep = state.preps.get(l);
          if (prep === undefined) state.preps.set(l, (prep = type.prepare(l, input.model)));
          c = type.content(prep, cctx);
        }
        if (c) {
          lanes = Math.max(lanes ?? 0, c.lanes);
          labels ||= c.labels;
        } else other = Math.max(other, byType[l.type as keyof RenderConfig['heights']] ?? 0);
      }
      // one lane block a track: the most lanes its layers need, one label lane for all of them
      contents.push(lanes === undefined ? undefined : { lanes, labels });
      sized.push({
        track: t,
        kind,
        fixed: fixed[j],
        ...(lanes !== undefined && rowsFill(t) === undefined
          ? { content: Math.max(other, lanesHeight(lanes + (labels ? 1 : 0), config)) }
          : {}),
        ...(t.height !== undefined ? { height: t.height } : {}),
      });
    }
    const name = section.params.heights;
    let sizing = typeof name === 'string' ? (sizings[name] ?? input.sizings?.[name]) : ContentSizing;
    let issue: Issue | undefined;
    if (!sizing) {
      issue = issueOf(
        RenderIssues,
        'seq.render.heights',
        `The track sizing ${name} is not registered here: tracks fit their content`,
        { address: section.address, hint: RenderIssues.codes['seq.render.heights'].hint },
      );
      sizing = ContentSizing;
    }
    let named = sizing.heights(sized, config);
    if (!heightsFit(named, sized.length)) named = ContentSizing.heights(sized, config);
    const ownHeights = named.map((h, j) => sized[j].height ?? h);
    const hook = input.trackSizing?.(section, sized, ownHeights);
    const heights = heightsFit(hook, sized.length) ? [...hook] : ownHeights;
    return { ...base, heights, contents, ...(issue ? { issue } : {}) };
  };
  const axes = new Map<string, AxisTransform | undefined>();
  const visiting = new Set<string>();
  /** A following ring's window: its leader's range on its own axis (made when a ring follows). */
  let followed: Map<string, Viewport> | undefined;
  const axisOf = (section: ResolvedSection): AxisTransform | undefined => {
    if (axes.has(section.address)) return axes.get(section.address);
    const at = section.address;
    let made: { pieces: readonly number[]; viewport?: Viewport } | undefined;
    if (section.follow) {
      const leaderSection = bySection.get(section.follow.leader);
      if (visiting.has(at) || !leaderSection) {
        const cycle = visiting.has(at);
        const code = cycle ? 'seq.render.follow.cycle' : 'seq.render.follow.leader';
        report(
          issueOf(
            RenderIssues,
            code,
            cycle ? 'The section follows itself through its leaders' : 'The section’s leader is not in the view',
            {
              address: at,
              hint: RenderIssues.codes[code].hint,
            },
          ),
        );
      } else if (
        (multiSystem(section) || multiSystem(leaderSection)) &&
        section.axis.system !== leaderSection.axis.system
      ) {
        // which block's map a follow would carry is undecided (_spec/README.md "Render", Blocks); a grid over a section
        // of blocks is on the blocks' own units, its follow the identity ("The grid layout", A grid over blocks)
        report(
          issueOf(RenderIssues, 'seq.render.follow.blocks', 'A section of blocks on several systems cannot follow', {
            address: at,
            hint: RenderIssues.codes['seq.render.follow.blocks'].hint,
          }),
        );
      } else {
        visiting.add(at);
        const leader = axisOf(leaderSection);
        visiting.delete(at);
        const own = ownPieces(section, input.model, input.data);
        if (leader) made = followAxis(section, own, leader, leaderWindow(leaderSection, leader), config);
        // a follower that draws its whole axis with a window (a ring): the leader's range is its window
        if (made && own && layoutOf(section).window) {
          followed ??= new Map();
          followed.set(at, followedWindow(section, own, made));
          made = { pieces: own };
        }
      }
    }
    if (!made) {
      const own = ownPieces(section, input.model, input.data);
      // a layout that reads its viewport its own way (text: an anchor; a ring: its window) places over the whole axis
      const viewport = layoutOf(section).followWindow ? undefined : (input.viewports?.[at] ?? input.viewport);
      if (own) made = { pieces: turned(section, own, viewport), ...(viewport ? { viewport } : {}) };
    }
    const axis = made?.pieces.length
      ? AxisTransform.on(section.axis, made.pieces, dataX, dataW, dpr, made.viewport)
      : undefined;
    if (!axes.has(at)) axes.set(at, axis);
    return axes.get(at);
  };
  for (const section of input.view.sections) axisOf(section);
  /**
   * A grid's y (_spec/README.md "Render", The grid layout): x's pieces and window when it follows x on x's system; a
   * leader's image through its path (x's own, or a section's); else its own window (`viewport.y`), all of it by default.
   */
  const yOf = (section: ResolvedSection, axis: AxisTransform | undefined): PlaceContext['y'] => {
    const y = section.y;
    const own = y?.axis.ranges;
    if (!y || !own || !axis) return undefined;
    if (y.followsX && !y.follow)
      return { axis: { ...y.axis, dir: axis.dir }, pieces: axis.pieces, viewport: { start: axis.u0, end: axis.u1 } };
    if (y.follow) {
      const leaderSection = y.followsX ? section : bySection.get(y.follow.leader);
      const leader = y.followsX ? axis : leaderSection && axes.get(leaderSection.address);
      if (leaderSection && leader) {
        const window = y.followsX ? { start: axis.u0, end: axis.u1 } : leaderWindow(leaderSection, leader);
        const made = followAxis({ ...section, axis: y.axis, follow: y.follow }, own, leader, window, config);
        if (made) return { axis: y.axis, pieces: made.pieces, ...(made.viewport ? { viewport: made.viewport } : {}) };
      }
      return { axis: y.axis, pieces: own };
    }
    const v = (input.viewports?.[section.address] ?? input.viewport)?.y;
    return { axis: y.axis, pieces: own, ...(v ? { viewport: v } : {}) };
  };
  const labelSize = fontSize(config.labelFont);
  const sections: SectionLayout[] = [];
  const used = new Set<string>();
  const showing = new Set<string>();
  // what a section's reuse compares besides its own viewport and reads (_spec/README.md "Render", The key)
  const stamp = viewStamp(input, config, dpr, overlay, state);
  const sinks: FrameSinks = { wants, report, theme, used, showing };
  const colors = defaultColors(input.view);
  let y = config.inset;
  for (const section of input.view.sections) {
    if (sections.length) y += config.sectionGap;
    const axis = axes.get(section.address);
    if (!axis)
      report(
        issueOf(RenderIssues, 'seq.render.axis.pending', 'The section’s length is not known yet', {
          address: section.address,
          hint: RenderIssues.codes['seq.render.axis.pending'].hint,
        }),
      );
    const layout = layoutOf(section);
    const viewport = input.viewports?.[section.address] ?? input.viewport;
    // an unchanged section is not laid out again: its placement and lays are drawn again (_spec/README.md "Render",
    // The key: per section)
    const kept = state.sections.get(section.address);
    // a grid's y: part of its key (a leader it follows may move it alone)
    const gridY = section.y ? yOf(section, axis) : undefined;
    if (kept && reusable(kept, section, layout, stamp, y, axis, viewport, input.data) && sameGridY(kept.gridY, gridY)) {
      kept.replay(sinks);
      if (layout.window) setWindow(layout, kept.layout.placement, section.address, viewport, followed, input, config);
      sections.push(kept.layout);
      y += kept.layout.height;
      continue;
    }
    // the kept layout's placement, read before the cache starts again
    const prevPlacement = kept?.layout?.layout === layout ? kept.layout.placement : undefined;
    const cache = kept ?? new SectionCache();
    cache.begin(input.data, sinks);
    const mine = cache.reads;
    const sectionReport = cache.report;
    const sectionTheme = cache.theme;
    const want = cache.want;
    const sectionTop = y;
    // a ruler rail on a linear section (tracks, not a stack of blocks); a slim row of its own when no ruler is shown
    const railable = layout === TracksLayout && !(section.blocks && section.params.arrange === 'stack');
    const railOn = input.rails?.[section.address] !== undefined && railable;
    const top = y;
    // a template is expanded here, over the overlay, into a track per node it selects; its source's region in view is
    // wanted, so the nodes it brings arrive
    const shown: ResolvedTrack[] = [];
    // each shown track's index in the section's (a template's instance: −1), what a block's layers are kept by, and
    // the section's track it came from (an instance's template)
    const trackOf: number[] = [];
    const sourceOf: ResolvedTrack[] = [];
    for (let ti = 0; ti < section.tracks.length; ti++) {
      const track = section.tracks[ti];
      if (!track.template) {
        trackOf[shown.length] = ti;
        sourceOf[shown.length] = track;
        shown.push(track);
        continue;
      }
      if (axis) wantTemplate(input, mine, section, track.template, axis, overlay, want);
      const instances = expandTemplate(input.model, overlay, section, track, config.maxTemplateTracks);
      if (instances.length >= config.maxTemplateTracks && templateNodes(track, overlay).length > instances.length)
        sectionReport(
          issueOf(
            RenderIssues,
            'seq.render.template.limit',
            `The template selects more than ${instances.length} nodes`,
            {
              address: track.address,
              hint: RenderIssues.codes['seq.render.template.limit'].hint,
            },
          ),
        );
      for (const t of instances) {
        // an instance whose node's region is out of view is not laid out: its track goes, its hold with it
        if (axis && !instanceInView(input, mine, section, track.template, axis, overlay, t, state)) continue;
        trackOf[shown.length] = -1;
        sourceOf[shown.length] = track;
        shown.push(t);
      }
    }
    const byType = layout.heights?.(config) ?? config.heights;
    const sizedBands = bandHeights(section, shown, layout, false, sectionReport);
    // the kept heights, copied only where this layout changes them (a rows band's)
    let heights = sizedBands.heights as number[];
    const trackContents = sizedBands.contents;
    // each track's padding inside its well, `[top, bottom]` a track: its own, else the config's (_spec/README.md
    // "Render", Tracks layout)
    const pads: number[] = [];
    let padded = false;
    for (const t of shown) {
      const p = t.padding ?? config.trackPadding;
      pads.push(p.top, p.bottom);
      if (p.top > 0 || p.bottom > 0) padded = true;
    }
    const ruler = railOn ? rulerTrack(shown, heights) : -1;
    // the rail draws just below the ruler's bottom border: the ruler's well grows by the rail's row, so it never covers
    // the ticks or labels; only where a rail is drawn live (the canvas with rails on): an export has no row, a figure's
    // rails (navigators included) lie on the border (_spec/README.md "Render", Navigators)
    if (ruler >= 0 && input.figure !== true) {
      pads[2 * ruler + 1] += SlimRailPx;
      padded = true;
    }
    // a rows band that fills (_spec/README.md "Render", The rows layer): the section's room left on the page, at least
    // its default height (a figure, laid out for its content, takes the default)
    // a band of every row in a layout whose bands do not scroll (text: an alignment figure's lines repeat its rows)
    if (layout !== TracksLayout)
      for (let j = 0; j < shown.length; j++) {
        if (shown[j].height !== undefined || rowsFill(shown[j]) === undefined) continue;
        const rows = rowsCount(shown[j], input.model);
        if (rows === undefined) continue;
        if (heights === sizedBands.heights) heights = heights.slice();
        heights[j] = rows * byType.rows + titleRowOf(shown[j], config);
      }
    const fill = shown.findIndex((t) => t.height === undefined && rowsFill(t) === 'fill');
    if (fill >= 0 && !input.figure && layout === TracksLayout) {
      let used = 0;
      for (let j = 0; j < heights.length; j++)
        used += (j !== fill ? heights[j] + config.trackGap : 0) + pads[2 * j] + pads[2 * j + 1];
      if (heights === sizedBands.heights) heights = heights.slice();
      heights[fill] = Math.max(heights[fill], Math.floor(height - top - used - config.inset));
    }
    const nb = shown.length;
    const lays: LayerLayout[][] = [];
    const stripBreaks: number[][] = [];
    // one layer's layout into its kept lay (its cache key: the layer's, by strip where a layout has several)
    const layOne = (
      track: ResolvedTrack,
      layer: ResolvedLayer,
      type: LayerType<any, any>,
      cacheKey: string,
      ctx: LayoutContext,
    ): LayerLayout => {
      let prep = state.preps.get(layer);
      if (prep === undefined) state.preps.set(layer, (prep = type.prepare(layer, input.model)));
      const keptLay = state.lays.get(cacheKey);
      const lay = type.layout(prep, ctx, mine, overlay, keptLay?.type === type ? keptLay.lay : undefined);
      state.lays.set(cacheKey, { type, lay });
      cache.addUsed(cacheKey);
      // what a shown layer draws and goes through (a node, an item system)
      if (layer.data) cache.addShowing(layer.data.address, layer.data.system);
      const edges = layer.path?.edges;
      if (edges) for (let e = 0; e < edges.length; e++) cache.addShowing(edges[e].from, edges[e].to);
      return { key: ctx.key, layer, type, lay };
    };
    /** A layer's type: the first of its fallbacks drawn here; none is said. */
    const typeOf = (layer: ResolvedLayer, key: string): LayerType<any, any> | undefined => {
      const type = firstType(layer, types);
      if (!type)
        sectionReport(
          issueOf(RenderIssues, 'seq.render.layer.type', 'No type of the layer’s fallbacks is drawn here', {
            address: key,
            hint: RenderIssues.codes['seq.render.layer.type'].hint,
          }),
        );
      return type;
    };
    /** Band `j`'s label lane, its layers' labels beside sharing it (kept in the state, one a band). */
    const labelLaneOf = (j: number): LabelLane => (state.labelLanes[j] ??= new LabelLane());
    /** A band's own scroll and row zoom (the viewport's `bands`, by track address). */
    const bandOf = (track: ResolvedTrack) => viewport?.bands?.[track.address];
    // one rows layer a section (the cursor has one row order): the first shown, a second said and not drawn
    let firstRows: string | undefined;
    for (const t of shown)
      for (const l of t.layers)
        if (firstRows === undefined && l.type === 'rows' && l.placed && l.params.hidden !== true)
          firstRows = layerKey(t, l);
    // the lays of a strip's bands (_spec/README.md "Render", Section layouts: the lay cache bounded by the strips in view);
    // a `section` layer is laid out once a section, after
    const layStrip = (placement: Placement, i: number, labels?: LabelSink): boolean => {
      let pill = false;
      const strips = placement.strips;
      const a = strips.axes[i];
      const keyed = layout !== TracksLayout || strips.n > 1 || section.blocks !== undefined;
      // a block's strip draws its layers resolved onto the block's system (a path a block)
      const block = section.blocks?.[strips.block[i]];
      const span = stripSpan(strips, i);
      for (let j = 0; j < nb; j++) {
        const k = i * nb + j;
        const out: LayerLayout[] = [];
        lays[k] = out;
        strips.used[k] = 0;
        const h = strips.bandH[k];
        // a band below the page is not laid out: nothing of it is drawn or asked for
        if (!(strips.bandY[k] + strips.dy[i] < height) || (h <= 0 && layout !== TracksLayout)) continue;
        const track = shown[j];
        const band = bandOf(track);
        if (trackContents[j]) labelLaneOf(j).reset();
        const box: Box = { x: a.x0, y: strips.bandY[k], width: a.width, height: h };
        const layers = block && trackOf[j] >= 0 ? block.layers[trackOf[j]] : track.layers;
        for (const layer of layers) {
          if (!layer.placed || layer.params.hidden === true) continue; // a hidden layer is not drawn
          const key = layerKey(track, layer);
          const type = typeOf(layer, key);
          // a layout of two axes lays every layer out a strip at a time, `section` layers too (pairs as dots)
          if (!type || (type.scope === 'section' && !layout.axesY)) continue;
          if (type.name === 'rows' && key !== firstRows) {
            sectionReport(
              issueOf(RenderIssues, 'seq.render.rows.second', 'A section draws one rows layer', {
                address: key,
                hint: RenderIssues.codes['seq.render.rows.second'].hint,
              }),
            );
            continue;
          }
          const ck = colorKey(sourceOf[j], layer);
          const colorIndex = ck === undefined ? undefined : colors.get(ck);
          const ctx: LayoutContext = {
            box,
            axis: a,
            section,
            key,
            dpr,
            config,
            tokens: input.tokens,
            text: input.text,
            settled: input.settled ?? false,
            theme: sectionTheme,
            want,
            report: sectionReport,
            viewWidth: width,
            ...(colorIndex !== undefined ? { colorIndex } : {}),
            ...(trackContents[j] ? { content: trackContents[j], labelLane: labelLaneOf(j) } : {}),
            strip: i,
            span: Array.isArray(layer.params.across) ? acrossSpan(strips, i, shown, layer.params.across) : span,
            place: placement,
            labels: labels ?? 'inside',
            ...(band ? { bandTop: band.bandTop } : {}),
            ...(band?.rowScale !== undefined ? { rowScale: band.rowScale } : {}),
            ...(band?.folds !== undefined
              ? { folds: band.folds }
              : input.folds?.[track.address]
                ? { folds: input.folds[track.address] }
                : {}),
            ...(input.images ? { images: input.images } : {}),
            ...(strips.ys[i] ? { axisY: strips.ys[i] } : {}),
            // a chip's pill: its label elided to its slot (a position's px on its strip: the whole chip's)
            ...(block && strips.chip[i] ? { chip: chipText(input.text, config, block.label, a.scale) } : {}),
            ...(grants && type.name === 'rows' && grants.has(key) ? { rowColumns: grants.get(key)! } : {}),
          };
          // kept by the strip's first unit, and its block's (a stack's blocks all start at 0)
          let at = block ? `${key}#${strips.block[i]}@${strips.first[i]}` : `${key}@${strips.first[i]}`;
          // a grid over blocks: a strip a pair of blocks, kept by its y block too
          if (block && strips.ys[i]) at += `|${strips.yBlock[i]}@${strips.yFirst[i]}`;
          const l = layOne(track, layer, type, keyed ? at : key, ctx);
          const u = type.used ? type.used(l.lay) : h;
          if (u > strips.used[k]) strips.used[k] = u;
          out.push(l);
          if (type.name === 'letters' && strips.chip[i]) {
            const ll = l.lay as LettersLay;
            if (ll.chip !== '' && ll.n > 0) pill = true;
          }
        }
      }
      return pill;
    };
    /** The `section` layers: once a section, over its box on the page, the placement their reader. */
    const laySection = (placement: Placement, a: AxisTransform): LayerLayout[] => {
      let out = NoSectionLays;
      if (layout.axesY) return out;
      for (let band = 0; band < shown.length; band++)
        for (const layer of shown[band].layers) {
          const track = shown[band];
          if (!layer.placed || layer.params.hidden === true) continue;
          const key = layerKey(track, layer);
          const type = firstType(layer, types);
          if (type?.scope !== 'section') continue;
          const ctx: LayoutContext = {
            box: placement.box,
            axis: a,
            section,
            key,
            dpr,
            config,
            tokens: input.tokens,
            text: input.text,
            settled: input.settled ?? false,
            theme: sectionTheme,
            want,
            report: sectionReport,
            span: placement.box,
            place: placement,
            labels: 'inside',
            band,
          };
          if (out === NoSectionLays) out = [];
          out.push(layOne(track, layer, type, `${key}@section`, ctx));
        }
      return out;
    };
    const sectionBox: Box = {
      x: dataX,
      y: top,
      width: Math.max(1, width - dataX - config.inset - right),
      height: Math.max(0, height - top),
    };
    let placement: Placement;
    let sectionLays: LayerLayout[] = NoSectionLays;
    if (axis) {
      const ctx: PlaceContext = {
        width: layout === TextLayout ? textWidth : width,
        height,
        dpr,
        config,
        tokens: input.tokens,
        text: input.text,
        model: input.model,
        box: sectionBox,
        axis,
        bands: shown,
        heights,
        ...(padded ? { pads } : {}),
        layStrip,
        report: sectionReport,
        ...(input.figure ? { figure: true } : {}),
        ...(input.windows ? { windows: true } : {}),
        ...(gridY ? { y: gridY } : {}),
      };
      placement = layout.place(section, ctx, viewport, prevPlacement);
      if (layout.window) setWindow(layout, placement, section.address, viewport, followed, input, config);
      for (let i = 0; i < placement.strips.n; i++) stripBreaks.push(placement.strips.axes[i].breaks([]));
      sectionLays = laySection(placement, axis);
    } else placement = pendingPlacement(top, sectionBox, heights, pads, config);
    // the section's tracks as bands of its first strip (on the page: tracks)
    const tracks: TrackLayout[] = [];
    const strips = placement.strips;
    let by = top;
    for (let j = 0; j < nb; j++) {
      if (j) by += config.trackGap;
      const track = shown[j];
      // a grid's track is its plot (over blocks: its cells together), its well the plot's grown by its padding
      const plot = placement instanceof GridPlacement && strips.n ? placement : undefined;
      const pt = plot ? plot.padTop : strips.n ? strips.padT[j] : pads[2 * j];
      const pb = plot ? plot.padBottom : strips.n ? strips.padB[j] : pads[2 * j + 1];
      const box: Box = {
        x: dataX,
        y: plot ? plot.plotTop : strips.n ? strips.bandY[j] + strips.dy[0] : by + pt,
        width: sectionBox.width,
        // a layout of two axes overlays its bands on the plot
        height: plot ? plot.plotHeight : heights[j],
      };
      const well: Box =
        pt > 0 || pb > 0 ? { x: box.x, y: box.y - pt, width: box.width, height: box.height + pt + pb } : box;
      by += pt + heights[j] + pb;
      const layers = lays[j] ?? NoLays;
      let filter: TrackLayout['filter'];
      // its layers', then its `section` layers' (pairs: laid out once a section)
      for (const own of [layers, sectionLays]) {
        for (const l of own) {
          if (own === sectionLays && !track.layers.includes(l.layer)) continue;
          const count = l.type.filtered?.(l.lay) ?? 0;
          if (!count) continue;
          const text = `${count} filtered`;
          const fw = input.text.measure(config.labelFont, text);
          const layer = own === sectionLays ? layerKey(track, l.layer) : l.key;
          // a track of two lines' height: the chip under the label, which keeps its width; else at the column's end,
          // the label fading before it
          const dy = box.height >= 2 * (labelSize + 6) ? labelSize + 4 : 0;
          const x = dy ? config.inset : config.inset + config.labelWidth - fw;
          filter = { layer, count, text, x, width: fw, dy };
          break;
        }
        if (filter) break;
      }
      // the label's room: the column, less a chip beside it (a chip on a line of its own leaves it whole)
      const room = filter && !filter.dy ? Math.max(0, filter.x - config.inset - 6) : config.labelWidth;
      const labelText = elidedLabel(state, input.text, config.labelFont, track.label, room);
      tracks.push({
        track,
        box,
        well,
        labelText,
        labelWidth: input.text.measure(config.labelFont, labelText),
        // centred on its well (its padding with it), not on its data
        labelY: well.y + well.height / 2 + labelSize * 0.35,
        layers,
        ...(filter ? { filter } : {}),
      });
    }
    y = top + placement.box.height;
    // the slim rail's row, under the tracks
    const slimTop = y;
    if (railOn && ruler < 0) y += SlimRailPx;
    // the column the emphasis spans with `emphasis.selectionShape: 'column'`: the tracks laid out on the page and the
    // gaps between them
    let bandEnd = -1;
    for (const t of tracks) if (t.box.y < height) bandEnd = Math.max(bandEnd, t.box.y + t.box.height);
    const first = tracks[0]?.box;
    const band =
      first && bandEnd > first.y
        ? { x: first.x, y: first.y, width: first.width, height: bandEnd - first.y }
        : undefined;
    // each band's box and each strip's band column, in the strip's frame (tracks: its first strip's are the page's)
    const boxes: Box[] = [];
    const columns: (Box | undefined)[] = [];
    for (let i = 0; i < strips.n; i++) {
      const a = strips.axes[i];
      for (let j = 0; j < nb; j++) {
        const k = i * nb + j;
        boxes.push(
          i === 0 && layout.fillBands && !layout.axesY
            ? tracks[j].box
            : { x: a.x0, y: strips.bandY[k], width: a.width, height: strips.bandH[k] },
        );
      }
      columns.push(i === 0 && layout.fillBands && !layout.axesY ? band : stripColumn(strips, i, height));
    }
    const text = tracks.length ? directionMark(section, config) : undefined;
    const mark = text
      ? {
          text,
          x: config.inset + tracks[0].labelWidth + 6,
          width: input.text.measure(config.labelFont, text),
        }
      : undefined;
    let rail: Rail | undefined;
    if (railOn && axis) {
      rail = kept?.layout?.rail ?? new Rail();
      const own = ownPieces(section, input.model, input.data) ?? axis.pieces;
      placeRail(
        rail,
        section,
        sectionId(section),
        axis,
        own,
        ruler >= 0 ? tracks[ruler].box : undefined,
        slimTop,
        viewport,
        config,
        input.figure === true,
      );
    }
    // each rows layer's rail, upright at the right edge (beside its last numbers), along its band's scrolled rows
    let rowsRails: Rail[] | undefined;
    if (input.rails && layout === TracksLayout)
      for (let j = 0; j < tracks.length; j++) {
        const t = tracks[j];
        if (input.rails[t.track.address] === undefined) continue;
        const l = t.layers.find((x) => x.lay instanceof RowsLay);
        const lay = l?.lay as RowsLay | undefined;
        if (!lay || lay.figure || lay.total === undefined) continue;
        const r = kept?.layout?.rowsRails?.[rowsRails?.length ?? 0] ?? new Rail();
        placeRowsRail(
          r,
          lay,
          t.track.address,
          sectionId(section),
          width - config.inset - SlimRailPx,
          input.figure === true,
        );
        (rowsRails ??= []).push(r);
      }
    const sectionLayout: SectionLayout = {
      section,
      id: sectionId(section),
      layout,
      placement,
      ...(axis ? { axis } : {}),
      ...(band ? { band } : {}),
      ...(rail ? { rail } : {}),
      ...(rowsRails ? { rowsRails } : {}),
      y: sectionTop,
      height: y - sectionTop,
      tracks,
      breaks: axis && band ? (stripBreaks[0] ?? NoBreaks) : NoBreaks,
      ...(mark ? { mark } : {}),
      lays,
      stripBreaks,
      boxes,
      columns,
      sectionLays,
      pageHeight: height,
      tokens: input.tokens,
      config,
    };
    sections.push(sectionLayout);
    cache.layout = sectionLayout;
    cache.section = section;
    cache.stamp = stamp;
    cache.y = sectionTop;
    cache.axis = axis;
    cache.gridY = gridY;
    // a copy: a caller may hand the same record again with new numbers (the viewer reuses its own)
    cache.viewport = viewport ? { ...viewport } : undefined;
    state.sections.set(section.address, cache);
  }
  // a section gone from the view lets its cache go
  if (state.sections.size > input.view.sections.length) {
    const here = new Set(input.view.sections.map((s) => s.address));
    for (const k of [...state.sections.keys()]) if (!here.has(k)) state.sections.delete(k);
  }
  // a layer gone from the view lets its last layout go
  for (const key of [...state.lays.keys()]) if (!used.has(key)) state.lays.delete(key);
  // open scales: settled, the union of this layout's extents is held (laid out once more with it where it moved);
  // moving, a scale seen for the first time takes this layout's extent for the next frame
  if (input.domains) {
    // the domains are another layout's: none held here
  } else if (input.settled) {
    let moved = extents.size !== state.domains.size;
    for (const [k, e] of extents) {
      const h = state.domains.get(k);
      if (!h || h[0] !== e[0] || h[1] !== e[1]) moved = true;
    }
    if (moved) {
      state.domains = extents;
      state.domainsVersion++;
      if (!again) return layoutView(input, state, true);
    }
  } else
    for (const [k, e] of extents)
      if (!state.domains.has(k)) {
        state.domains.set(k, e);
        state.domainsVersion++;
      }
  for (const i of overlay.issues) report(i);
  const shownNodes = config.nodes.holdShown && overlay.nodes.length ? overlay.holds(showing) : [];
  return { width, height, dpr, tokens: input.tokens, config, sections, wants, issues, shown: shownNodes };
}

/**
 * A template source's region in view, on the section's system as the source names it, at the level that fits; and for
 * a template below its source, each nested source's systems' regions in view (through the overlay's path from the
 * section), so a nested source that sends its nodes by region is asked for them.
 */
function wantTemplate(
  input: LayoutInput,
  data: StoreReader,
  section: ResolvedSection,
  template: TrackTemplate,
  axis: AxisTransform,
  overlay: Overlay,
  out: (w: Want) => void,
) {
  const pieces = axis.visible([]);
  const at = templateSource(input, data, section, template, axis);
  if (at)
    for (let i = 0; i + 1 < pieces.length; i += 2)
      out({
        source: template.source,
        system: at.system,
        start: Math.floor(pieces[i]),
        end: Math.ceil(pieces[i + 1]),
        level: at.level,
      });
  if (template.depth < 1) return;
  for (const nested of overlay.nested(new Map([[template.source, template.depth]]))) {
    const cov = data.coverage(nested.address);
    for (const s of cov?.systems ?? []) {
      const target = overlay.systemNamed(s.name, nested.parent);
      if (!target) continue;
      const path =
        target.address === section.axis.system ? undefined : overlay.path(section.axis.system, target.address).path;
      if (!path && target.address !== section.axis.system) continue;
      const ranges = HalfOpenRanges.normalize(path ? projectRanges(path, pieces) : pieces);
      const ratio = path?.step === 1 ? 3 : path?.step === -1 ? 1 / 3 : 1;
      // bounded as a layer's level is: a long nested system at the fit is one coarse tile, never 81,307 at level 1
      const level = boundedLevel(cov!.levels, ratio / axis.devicePx, (axis.u1 - axis.u0) * ratio, data);
      for (let i = 0; i + 1 < ranges.length; i += 2)
        out({
          source: nested.address,
          system: s.name,
          start: Math.floor(ranges[i]),
          end: Math.ceil(ranges[i + 1]),
          level,
        });
    }
  }
}

/**
 * A template source's system as it names the section's, its coverage, and the level the template wants of it: the one
 * that fits a device pixel, coarser while the view would want more than the store's tiles a want (`boundedLevel`).
 */
function templateSource(
  input: LayoutInput,
  data: StoreReader,
  section: ResolvedSection,
  template: TrackTemplate,
  axis: AxisTransform,
): { readonly coverage: Coverage; readonly system: string; readonly level: number } | undefined {
  const coverage = data.coverage(template.source);
  const name = input.model.systems.get(section.axis.system)?.name;
  const info = coverage?.systems.find((s) => s.name === name || (name !== undefined && s.aliases?.includes(name)));
  if (!coverage || !info) return undefined;
  return {
    coverage,
    system: info.name,
    level: boundedLevel(coverage.levels, 1 / axis.devicePx, axis.u1 - axis.u0, data),
  };
}

/**
 * Whether a template's instance is in view (_spec/README.md "Render", Shown nodes): its node's top node (the one of
 * the template's source it came under) is listed up front in that source's coverage, or carried by a loaded tile of it
 * over the view at any level whose tiles there are within the store's bound; while the template's tiles over the view
 * are pending (a zoom's new level), an instance shown last time stays. Out of view it is not laid out: its track goes
 * and its hold with it, and its node is evicted as any other once the store needs room. Where the template's source
 * does not name the section's system it cannot be told: shown.
 */
function instanceInView(
  input: LayoutInput,
  data: StoreReader,
  section: ResolvedSection,
  template: TrackTemplate,
  axis: AxisTransform,
  overlay: Overlay,
  instance: ResolvedTrack,
  state: LayoutState,
): boolean {
  const at = templateSource(input, data, section, template, axis);
  let address: string | undefined;
  for (const l of instance.layers) if (l.data && !address) address = l.data.address;
  let root = address ? overlay.node(address) : undefined;
  while (root && root.depth > 0) root = overlay.node(root.source);
  if (!at || !root) return true;
  const top = root;
  if (at.coverage.nodes?.some((x) => x.key === top.key)) return true;
  const key = nodeKey(top.source, top.key);
  const pieces = axis.visible([]);
  let pending = false;
  for (const level of at.coverage.levels)
    for (let i = 0; i + 1 < pieces.length; i += 2) {
      const want = {
        source: template.source,
        system: at.system,
        start: Math.floor(pieces[i]),
        end: Math.ceil(pieces[i + 1]),
        level,
      };
      // bounded by what it reads: a level whose tiles over the view pass the store's bound is not read
      if (data.tooMany?.(want)) continue;
      for (const k of data.keysOf(want)) {
        const e = data.get(k);
        if (e?.status === 'loaded') {
          if (e.nodes.includes(key)) return remember(state, instance.address, true);
        } else if (level === at.level && e?.status !== 'empty' && e?.status !== 'failed') pending = true;
      }
    }
  return remember(state, instance.address, pending && state.instancesShown.has(instance.address));
}

/** Records whether an instance was shown, for the next frame's pending rule; returns it. */
function remember(state: LayoutState, address: string, shown: boolean): boolean {
  if (shown) state.instancesShown.add(address);
  else state.instancesShown.delete(address);
  return shown;
}

/** The first shown track drawing a ruler (with height: a hidden ruler is no rail's), or −1. */
function rulerTrack(shown: readonly ResolvedTrack[], heights: readonly number[]): number {
  for (let j = 0; j < shown.length; j++)
    if (heights[j] > 0 && shown[j].layers.some((l) => l.type === 'ruler' && l.placed && l.params.hidden !== true))
      return j;
  return -1;
}

/** How a track's rows layer takes its height: `fill`, a number of rows, or none (no rows layer). */
function rowsFill(track: ResolvedTrack): 'fill' | number | undefined {
  for (const l of track.layers) {
    if (l.type !== 'rows') continue;
    const h = l.params.height;
    return typeof h === 'number' ? h : 'fill';
  }
  return undefined;
}

/** The rows a track's rows layer shows, when its alignment's length is known (an inline one's). */
function rowsCount(track: ResolvedTrack, model: StaticModel): number | undefined {
  for (const l of track.layers) {
    if (l.type !== 'rows' || !l.data) continue;
    const rows = model.annotations.get(l.data.address)?.rows;
    return rows === undefined ? undefined : model.systems.get(rows)?.length;
  }
  return undefined;
}

function trackHeight(track: ResolvedTrack, config: RenderConfig, heights: RenderConfig['heights']): number {
  if (track.height !== undefined) return track.height;
  const rows = rowsFill(track);
  if (typeof rows === 'number') return rows * RowHeight + titleRowOf(track, config);
  let h = 0;
  for (const l of track.layers) h = Math.max(h, heights[l.type as keyof RenderConfig['heights']] ?? 0);
  return h || heights.ranges;
}

const NoBreaks: readonly number[] = [];

/** Draws the layout's passes (all, or those named, e.g. the emphasis alone) into the lists, cleared first. */
export function drawView(
  layout: ViewLayout,
  lists: DrawLists,
  emphasis?: Emphasis,
  passes: readonly PassName[] = Passes,
): void {
  lists.reset(passes);
  lists.width = layout.width;
  lists.height = layout.height;
  // indexed loops: an iterator a loop is an allocation where the compiler keeps it (nothing allocated once warm)
  for (let pi = 0; pi < passes.length; pi++) {
    const pass = passes[pi];
    drawPass(ViewDraws, layout, lists.writer(pass), emphasis);
  }
}

/** The view's own draw of each pass, through `drawPass` like a layer type's. */
const ViewDraws = passDraws<ViewLayout, Emphasis | undefined>({
  background: drawBackground,
  data: drawTracks,
  window: drawTracks,
  highlight: drawTracks,
  selection: drawTracks,
  hover: drawTracks,
  cursor: drawTracks,
  labels: drawTracks,
  tip: drawTracks,
});

/**
 * The background pass: the page, each track's well (a layout that fills its bands: a squircle by the squircle rule,
 * `trackShape`), and with the tints under the data (`emphasis.selectionPlace: 'under'`) the highlight channels' tints
 * and then the selection's above them, a band a section or a box a track (`emphasis.selectionShape`: inside the well's
 * clip), a strip at a time through its warp.
 */
function drawBackground(layout: ViewLayout, out: PassWriter, emphasis: Emphasis | undefined) {
  const { tokens, config } = layout;
  out.fill(0, 0, layout.width, layout.height, tokens.colors['msui-bg']);
  const under = config.emphasis.selectionPlace === 'under' && emphasis !== undefined;
  const column = config.emphasis.selectionShape === 'column';
  for (let si = 0; si < layout.sections.length; si++) {
    const s = layout.sections[si];
    if (s.layout.fillBands) {
      for (let ti = 0; ti < s.tracks.length; ti++) {
        const t = s.tracks[ti];
        const { box, well } = t;
        if (box.y >= layout.height) continue;
        // the well through `args` (nothing boxed): its squircle by its height
        const g = out.args;
        g[0] = well.x;
        g[1] = well.y;
        g[2] = well.width;
        g[3] = well.height;
        wellCorners(g, 7, config);
        g[4] = tokens.colors['seq-track'];
        g[5] = NoColor;
        g[6] = 0;
        out.rectArgs();
        if (under && !column && s.axis) {
          const strips = s.placement.strips;
          if (strips.n) {
            // in the page's frame: the first strip's (a layout that fills its bands translates it by nothing)
            wellClipArgs(out.args, s, 0, ti, box);
            out.args[0] += strips.dx[0];
            out.args[1] += strips.dy[0];
            out.clipSquircleArgs();
            EmphasisPart = Part.tints;
          }
          drawTints(out, s, s.axis, box, emphasis!, tokens, config);
          if (strips.n) {
            out.unclip();
            EmphasisPart = Part.outlines;
            drawTints(out, s, s.axis, box, emphasis!, tokens, config);
            EmphasisPart = Part.all;
          }
        }
      }
      if (under && column && s.axis && s.band) drawTints(out, s, s.axis, s.band, emphasis!, tokens, config);
      continue;
    }
    if (!under) continue;
    const strips = s.placement.strips;
    const pb = s.placement.box;
    if (s.placement.clip === true) clipTo(out, pb);
    for (let i = 0; i < strips.n; i++) {
      const w = openStrip(s, i, out);
      const a = strips.axes[i];
      const col = s.columns[i];
      if (column && col) drawTints(w, s, a, col, emphasis!, tokens, config);
      else if (!column)
        for (let j = 0; j < strips.bands; j++) {
          const box = s.boxes[i * strips.bands + j];
          if (box.height > 0) drawTints(w, s, a, box, emphasis!, tokens, config);
        }
      closeStrip(s, i, out);
    }
    if (s.placement.clip === true) out.unclip();
  }
}

/** The tints under the data: the highlight channels', then the selection's. */
function drawTints(
  out: PassWriter,
  s: SectionLayout,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  if (emphasis.highlights) drawHighlights(out, axis, box, emphasis, config);
  if (EmphasisPart !== Part.outlines)
    drawSelectionTint(out, axis, box, emphasis, tokens, config, projectedSides(s, axis, emphasis));
}

/**
 * Which of the emphasis a box draws: all of it, its tints alone (inside a track's well clip), or its outlines alone
 * (after the clip: an outline around a mark is never clipped).
 */
const Part = { all: 0, tints: 1, outlines: 2 } as const;
let EmphasisPart: number = Part.all;

/** Whether a section draws the selection's projected parts: not a selected rect's sides where it shows both systems. */
function projectedSides(s: SectionLayout, axis: AxisTransform, emphasis: Emphasis): boolean {
  return !rectOver(emphasis.selection, axis.system, s.section.y?.axis.system);
}

/** The strips' translation (one at a time: drawing is not reentrant). */
const Translate = new TranslateWarp();

/** The writer a strip draws into for `out`'s pass: its polar warp's, else `out` with the translation marked. */
/** A clip to a box through the writer's `args` (its fractional numbers passed to a call would be boxed). */
function clipTo(out: PassWriter, b: Box): void {
  const g = out.args;
  g[0] = b.x;
  g[1] = b.y;
  g[2] = b.width;
  g[3] = b.height;
  out.clipArgs();
}

function openStrip(s: SectionLayout, i: number, out: PassWriter): PassWriter {
  const strips = s.placement.strips;
  if (strips.warp[i] === WarpKind.polar) return s.placement.warpOf!(i).writer(out);
  Translate.open(out);
  return out;
}

/** What a strip drew, onto the page through its warp. */
function closeStrip(s: SectionLayout, i: number, out: PassWriter): void {
  const strips = s.placement.strips;
  if (strips.warp[i] === WarpKind.polar) s.placement.warpOf!(i).close(out);
  else Translate.close(out, strips.dx, strips.dy, i);
}

/**
 * Any other pass (the writer's): per section its layout's chrome, then per strip (through its warp) its emphasis marks
 * of the pass over its band column (`column`) or each band (`tracks`), each band's chrome (its label) and its layers'
 * marks, and its breaks.
 */
function drawTracks(layout: ViewLayout, out: PassWriter, emphasis: Emphasis | undefined) {
  const { config } = layout;
  const pass = out.pass.name;
  const column = config.emphasis.selectionShape === 'column';
  for (let si = 0; si < layout.sections.length; si++) {
    const s = layout.sections[si];
    // the tip is a rows layer's (a hovered cell read whole): a section with none has nothing to walk, as the window pass
    if (pass === 'tip' && !hasRows(s) && s.layout.name !== 'rings') continue;
    s.layout.chrome(s, pass, out, emphasis, -1, -1);
    if (pass === 'window' && s.rail) drawRail(s.rail, out, emphasis, s.tokens, s.config);
    if (pass === 'window' && s.rowsRails)
      for (let i = 0; i < s.rowsRails.length; i++) drawRail(s.rowsRails[i], out, emphasis, s.tokens, s.config);
    const strips = s.placement.strips;
    if (!strips.n) {
      // pending: its tracks' chrome alone
      for (let ti = 0; ti < s.tracks.length; ti++)
        if (s.tracks[ti].box.y < layout.height) s.layout.chrome(s, pass, out, emphasis, 0, ti);
      continue;
    }
    drawStrips(layout, s, out, emphasis, pass, column);
    // what crosses strips, on the page; a section that scrolls inside its box is clipped to it
    const sl = s.sectionLays;
    if (!sl.length) continue;
    const clipped = s.placement.clip === true;
    if (clipped) clipTo(out, s.placement.box);
    for (let li = 0; li < sl.length; li++) sl[li].type.draw(sl[li].lay, pass, out, emphasis);
    if (clipped) out.unclip();
  }
}

/** Whether a section has a rows layer (kept by the layout: its layers do not change within it). */
const RowsIn = new WeakMap<SectionLayout, boolean>();
function hasRows(s: SectionLayout): boolean {
  let has = RowsIn.get(s);
  if (has === undefined) {
    has = false;
    for (const band of s.lays) for (const l of band ?? []) if (l.type.name === 'rows') has = true;
    RowsIn.set(s, has);
  }
  return has;
}

/** A section's strips of a pass, each through its warp (`drawTracks`). */
function drawStrips(
  layout: ViewLayout,
  s: SectionLayout,
  out: PassWriter,
  emphasis: Emphasis | undefined,
  pass: PassName,
  column: boolean,
) {
  const { tokens } = layout;
  const strips = s.placement.strips;
  const nb = strips.bands;
  // a section that scrolls inside its box: each strip clipped to it, in the strip's frame (its warp moves the clip)
  const clipped = s.placement.clip === true;
  const pb = s.placement.box;
  {
    for (let i = 0; i < strips.n; i++) {
      const w = openStrip(s, i, out);
      const a = strips.axes[i];
      const col = s.columns[i];
      const dy = strips.dy[i];
      if (clipped) {
        const g = w.args;
        g[0] = pb.x - strips.dx[i];
        g[1] = pb.y - dy;
        g[2] = pb.width;
        g[3] = pb.height;
        w.clipArgs();
      }
      if (column && emphasis && col) drawSectionEmphasis(w, s, a, col, emphasis, layout);
      for (let j = 0; j < nb; j++) {
        const k = i * nb + j;
        const box = s.boxes[k];
        if (box.y + dy >= layout.height) continue;
        if (box.height <= 0 && !s.layout.fillBands) continue;
        s.layout.chrome(s, pass, w, emphasis, i, j);
        if (pass !== 'labels' && !column && emphasis) {
          // a box a track: its tints inside its well (the well's shape), the outlines around marks unclipped
          const tint =
            s.layout.fillBands &&
            ((pass === 'highlight' && emphasis.highlights !== undefined) ||
              (pass === 'selection' &&
                (emphasis.selection !== undefined || emphasis.projected?.selection !== undefined)));
          if (tint) {
            // the tints inside the well's clip, then the outlines over them unclipped
            wellClipArgs(w.args, s, i, j, box);
            w.clipSquircleArgs();
            EmphasisPart = Part.tints;
            drawSectionEmphasis(w, s, a, box, emphasis, layout);
            w.unclip();
            EmphasisPart = Part.outlines;
            drawSectionEmphasis(w, s, a, box, emphasis, layout);
            EmphasisPart = Part.all;
          } else drawSectionEmphasis(w, s, a, box, emphasis, layout);
        }
        const lays = s.lays[k];
        // the window pass is the layouts' own (their chrome): no layer draws there, so a hover's redraw walks none
        if (!lays || !lays.length || pass === 'window') continue;
        // the data clipped to its well: one clip a track (its squircle in a layout that fills its bands)
        const clip = pass === 'data';
        if (clip) {
          wellClipArgs(w.args, s, i, j, box);
          w.clipSquircleArgs();
        }
        for (let li = 0; li < lays.length; li++) {
          const l = lays[li];
          l.type.draw(l.lay, pass, w, emphasis);
        }
        if (clip) w.unclip();
      }
      const breaks = s.stripBreaks[i];
      if (pass === 'data' && col && breaks.length) drawBreaks(w, breaks, a, col, tokens);
      if (clipped) w.unclip();
      closeStrip(s, i, out);
    }
  }
}

/**
 * A strip's breaks (_spec/README.md "Render", Tracks layout): a hairline a device pixel wide in the break ink across
 * the band column, over the cells either side (they touch), and a notch above it. Zero units wide: the hairline has no
 * width of its own on the axis.
 */
function drawBreaks(out: PassWriter, b: readonly number[], axis: AxisTransform, band: Box, tokens: SeqTokens) {
  // through `args` (a number passed to a call that is not inlined is boxed: nothing allocated once warm)
  const g = out.args;
  const dpr = axis.dpr;
  const px = 1 / dpr;
  for (let i = 0; i + 2 < b.length; i += 3) {
    const x = Math.round(b[i] * dpr) / dpr;
    for (let k = 0; k < 2; k++) {
      g[0] = k === 0 ? x - px / 2 : x - Notch.width / 2;
      g[1] = k === 0 ? band.y : band.y - Notch.height;
      g[2] = k === 0 ? px : Notch.width;
      g[3] = k === 0 ? band.height : Notch.height;
      g[4] = tokens.colors['seq-break'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
  }
}

/**
 * A section's emphasis marks of the writer's pass over `box` (a strip's band column, or one band's box) on the strip's
 * axis: the highlights' tints (over the data), the selection's tint and outline, the hover carried here from a linked
 * section (never on the hover's own section: its layer outlines it), the keyboard cursor.
 */
function drawSectionEmphasis(
  out: PassWriter,
  s: SectionLayout,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  layout: ViewLayout,
) {
  const { tokens, config } = layout;
  const over = config.emphasis.selectionPlace === 'over';
  switch (out.pass.name) {
    case 'highlight':
      // over the data, below the selection's pass; under it, in the background pass
      if (over && emphasis.highlights) drawHighlights(out, axis, box, emphasis, config);
      break;
    case 'selection': {
      const sides = projectedSides(s, axis, emphasis);
      if (over && EmphasisPart !== Part.outlines) drawSelectionTint(out, axis, box, emphasis, tokens, config, sides);
      if (EmphasisPart !== Part.tints) drawSelectionOutline(out, axis, box, emphasis, tokens, config, sides);
      break;
    }
    case 'hover':
      if (emphasis.projected?.hover && emphasis.projected.hoverSection !== s.section.address)
        drawProjectedHover(out, axis, box, emphasis.projected.hover, tokens, config);
      break;
    case 'cursor':
      // a cursor on a cell of two axes is its layer's on its own section (it alone knows its y); the others mark its
      // column, unless it stands on a row's header (the column it left is no place of it)
      if (
        emphasis.cursor?.system === axis.system &&
        emphasis.cursor.part !== 'header' &&
        (emphasis.cursor.ySystem === undefined || emphasis.cursor.section !== s.section.address)
      )
        drawCursor(out, axis, box, emphasis.cursor, tokens, config);
      break;
  }
}

// The cursor and the emphasis marks take the layout, not its numbers: a number passed to a call that is not inlined
// is boxed (nothing allocated once warm). What a shared mark function needs besides its objects (a fill, a stroke, a
// width, the corners, an offset) goes through `Paint`, a typed array: its numbers are never boxed; colors at an
// opacity are made inline (`withAlpha` would box its result).

/**
 * The keyboard cursor: a thin outline of its cell over the box, in a neutral ink (`emphasis.cursorInk`); a squircle by
 * the cell over the box (every emphasis a squircle), inside it by half its width (concentric), square on a side the
 * view cuts and below the minimum.
 */
function drawCursor(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  cursor: { readonly position: number },
  tokens: SeqTokens,
  config: RenderConfig,
) {
  if (!axis.cell(cursor.position)) return;
  const w = config.outline.cursor;
  // `snap` inline, the rect through `args`: an emphasis redraw allocates nothing once warm
  const dpr = axis.dpr;
  let x0 = Math.round(axis.sx0 * dpr) / dpr;
  let x1 = Math.round(axis.sx1 * dpr) / dpr;
  // a side past the view's edge is cut there: square
  const round = (x0 >= axis.x0 ? Round.left : 0) | (x1 <= axis.x0 + axis.width ? Round.right : 0);
  x0 = Math.max(x0, axis.x0);
  x1 = Math.min(x1, axis.x0 + axis.width);
  const cw = x1 - x0;
  const h = box.height;
  // `squircleRadius` inline (a number returned from a call is boxed); the outline's less half its width
  const r =
    cw < config.squircleMin || h < config.squircleMin
      ? 0
      : Math.min(h * config.squircleRatio, config.squircleMax, cw / 2, h / 2);
  const g = out.args;
  g[0] = x0 + w / 2;
  g[1] = box.y + w / 2;
  g[2] = Math.max(cw - w, 0.5);
  g[3] = h - w;
  g[4] = NoColor;
  g[5] = tokens.colors[config.emphasis.cursorInk];
  g[6] = w;
  g[7] = r > 0 && round !== 0 ? Math.max(r - w / 2, 0) : 0;
  g[8] = round;
  out.rectArgs();
}

/**
 * What `paintRanges` draws with: the fill, the stroke and its width, whether the ends are squircles where they really
 * end (1) or square (0), how far the stroke stands outside the range (0: inside it, concentric with the fill), and
 * which of the range's ends are no real end (a highlight's run cut where the stripes change: 1 its start, 2 its end).
 */
const Paint = new Float64Array(6);
/** Whether a range wholly in a hidden stretch lights its break (a highlight's later runs of it do not: once). */
let LitBreaks = true;

/**
 * A part's ranges on a box: per range and piece of the axis in view, a rect over the box's height snapped to the device
 * grid, filled with `Paint[0]` (none: `NoColor`), then outlined with `Paint[1]` at `Paint[2]` px, inside it or
 * `Paint[4]` px outside it, concentric. With `Paint[3]`, an end the range really has is a squircle (a third of the
 * height, capped) and an end cut by the view or a piece is square (the frame rule).
 */
function paintRanges(out: PassWriter, axis: AxisTransform, box: Box, part: LocationPart, config: RenderConfig) {
  const fill = Paint[0];
  const stroke = Paint[1];
  const w = Paint[2];
  const rounded = Paint[3] !== 0;
  const outset = Paint[4];
  const dpr = axis.dpr;
  const g = out.args;
  const h = box.height;
  for (let i = 0; i + 1 < part.ranges.length; i += 2) {
    axis.piecesIn(part.ranges[i], part.ranges[i + 1]);
    const ka = axis.ka;
    const kb = axis.kb;
    let drew = false;
    for (let k = ka; k <= kb; k++) {
      if (!axis.spanIn(k, part.ranges[i], part.ranges[i + 1])) continue;
      drew = true;
      const x0 = Math.round(axis.sx0 * dpr) / dpr;
      const x1 = Math.max(Math.round(axis.sx1 * dpr) / dpr, x0 + 1 / dpr);
      // a run's ends that are no real end are cut (on the screen's sides: a reversed axis's start on the right)
      const open = Paint[5];
      let cut = axis.scut;
      if ((open & 1) !== 0) cut |= axis.dir === -1 ? Cut.right : Cut.left;
      if ((open & 2) !== 0) cut |= axis.dir === -1 ? Cut.left : Cut.right;
      const bw = x1 - x0;
      // `squircleRadius` inline (a number returned from a call is boxed)
      const r =
        !rounded || bw < config.squircleMin || h < config.squircleMin
          ? 0
          : Math.min(h * config.squircleRatio, config.squircleMax, bw / 2, h / 2);
      const round = r > 0 ? ((cut & Cut.left) === 0 ? Round.left : 0) | ((cut & Cut.right) === 0 ? Round.right : 0) : 0;
      if (fill !== NoColor) {
        g[0] = x0;
        g[1] = box.y;
        g[2] = bw;
        g[3] = h;
        g[4] = fill;
        g[5] = NoColor;
        g[6] = 0;
        g[7] = round ? r : 0;
        g[8] = round;
        out.rectArgs();
      }
      if (stroke !== NoColor) {
        // inside: inset by half the width; outside: out by the offset and half the width; the radius follows
        const d = outset > 0 ? -(outset + w / 2) : w / 2;
        g[0] = x0 + d;
        g[1] = box.y + d;
        g[2] = Math.max(bw - 2 * d, 0.5);
        g[3] = Math.max(h - 2 * d, 0.5);
        g[4] = NoColor;
        g[5] = stroke;
        g[6] = w;
        g[7] = round ? Math.max(r - d, 0) : 0;
        g[8] = round;
        out.rectArgs();
      }
    }
    // a range wholly in a hidden stretch lights its break: what a break hides is not dropped
    if (!drew && LitBreaks) litBreak(out, axis, box, part.ranges[i]);
  }
}

/**
 * The break hiding position `p`, when in view: a bar two device pixels wide across the box in the paint's fill, else its
 * stroke (read from `Paint`: a color passed to a call that is not inlined is boxed).
 */
function litBreak(out: PassWriter, axis: AxisTransform, box: Box, p: number) {
  const b = axisBreaks(axis);
  const i = b.length ? breakHolding(axis, p) : -1;
  if (i < 0 || b[i] < axis.u0 || b[i] > axis.u1) return;
  const px = 1 / axis.dpr;
  const g = out.args;
  g[0] = Math.round(axis.xOfUnit(b[i]) * axis.dpr) / axis.dpr - px;
  g[1] = box.y;
  g[2] = 2 * px;
  g[3] = box.height;
  g[4] = Paint[0] !== NoColor ? Paint[0] : Paint[1];
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** A Location's part on a system, found with an indexed loop (no iterator, no closure). */
function partOn(loc: Location | undefined, system: string): LocationPart | undefined {
  if (!loc) return undefined;
  for (let i = 0; i < loc.parts.length; i++) if (loc.parts[i].system === system) return loc.parts[i];
  return undefined;
}

/**
 * The selection's accent tint over a box's positions: its own part at `emphasis.selectionTint`, else the part
 * projected here at that times `emphasis.projected` (the same look, lighter). Squircle ends where it really ends.
 */
function drawSelectionTint(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  tokens: SeqTokens,
  config: RenderConfig,
  projected: boolean,
) {
  const own = partOn(emphasis.selection, axis.system);
  const part = own ?? (projected ? partOn(emphasis.projected?.selection, axis.system) : undefined);
  if (!part) return;
  const accent = tokens.colors['msui-accent'];
  const a = config.emphasis.selectionTint * (own ? 1 : config.emphasis.projected);
  Paint[0] = ((accent & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;
  Paint[1] = NoColor;
  Paint[2] = 0;
  Paint[3] = 1;
  Paint[4] = 0;
  Paint[5] = 0;
  paintRanges(out, axis, box, part, config);
}

/**
 * The selection's one thin outline around each range on a box, inside its tint: the accent at
 * `emphasis.selectionOutline`, its own part, else the projected one, lighter.
 */
function drawSelectionOutline(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  emphasis: Emphasis,
  tokens: SeqTokens,
  config: RenderConfig,
  projected: boolean,
) {
  const own = partOn(emphasis.selection, axis.system);
  const part = own ?? (projected ? partOn(emphasis.projected?.selection, axis.system) : undefined);
  if (!part) return;
  const accent = tokens.colors['msui-accent'];
  const a = ((accent & 0xff) / 255) * config.emphasis.selectionOutline * (own ? 1 : config.emphasis.projected);
  Paint[0] = NoColor;
  Paint[1] = ((accent & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;
  Paint[2] = config.outline.selection;
  Paint[3] = 1;
  Paint[4] = 0;
  Paint[5] = 0;
  paintRanges(out, axis, box, part, config);
}

/**
 * The hover carried onto a box from a linked section: the hover's outline around each range, `emphasis.hoverOffset`
 * outside it (so it never touches the selection's), lighter.
 */
function drawProjectedHover(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  projected: Location,
  tokens: SeqTokens,
  config: RenderConfig,
) {
  const part = partOn(projected, axis.system);
  if (!part) return;
  const hover = tokens.colors['seq-hover'];
  const a = ((hover & 0xff) / 255) * config.emphasis.projected;
  Paint[0] = NoColor;
  Paint[1] = ((hover & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;
  Paint[2] = config.outline.hover;
  // every emphasis a squircle (_spec/README.md "Render", Shapes): round where it really ends
  Paint[3] = 1;
  Paint[4] = config.emphasis.hoverOffset;
  Paint[5] = 0;
  paintRanges(out, axis, box, part, config);
}

/** A mutable box: one stripe of a band (reused: a stripe makes no object). */
class StripeBox implements Box {
  x = 0;
  y = 0;
  width = 0;
  height = 0;
}
const Stripe = new StripeBox();
/** The one `[start, end)` range of a highlight run, as a part (reused). */
const RunRanges = new Int32Array(2);
const RunPart: LocationPart = { system: '', ranges: RunRanges };
/** The channels with a part on this axis that take a stripe (fills), and those that do not (outlines): by index. */
const HParts: (LocationPart | undefined)[] = [];
const HIndex: number[] = [];
const HOutlines: number[] = [];
/** What `emitRun` draws on (set before `walkStripes`: it is given no closure). */
let HOut!: PassWriter;
let HAxis!: AxisTransform;
let HBox!: Box;
let HList!: readonly EmphasisHighlight[];
let HConfig!: RenderConfig;

/**
 * The highlight channels on a box, in order (_spec/README.md "Render", Emphasis): where one channel alone lights a
 * position it takes the whole band; where several do, the band is split into equal stripes, one a channel present, in
 * channel order, at most `emphasis.highlightStripes` (the rest share the last one). A channel's `style` is drawn in its
 * stripe: `tint` (its color at `emphasis.highlightTint`) or `band` (a bar along the stripe's foot); an `outline` takes
 * no stripe (it has no fill): one outline over the channel's own ranges on the whole band, above the fills.
 */
function drawHighlights(out: PassWriter, axis: AxisTransform, box: Box, emphasis: Emphasis, config: RenderConfig) {
  const hs = emphasis.highlights;
  if (!hs) return;
  let n = 0;
  let o = 0;
  for (let i = 0; i < hs.length; i++) {
    const part = partOn(hs[i].location, axis.system);
    if (!part) continue;
    if (hs[i].style === 'outline') HOutlines[o++] = i;
    else {
      HParts[n] = part;
      HIndex[n] = i;
      n++;
    }
  }
  // the fills (tints, bands) and the outlines drawn apart when a well's clip holds the fills alone
  const fills = EmphasisPart === Part.outlines ? 0 : n;
  if (EmphasisPart === Part.tints) o = 0;
  if (fills === 1) paintHighlight(out, axis, box, hs[HIndex[0]], HParts[0]!, 0, 1, config);
  else if (fills > 1) {
    HOut = out;
    HAxis = axis;
    HBox = box;
    HList = hs;
    HConfig = config;
    walkStripes(HParts, n, Math.max(1, config.emphasis.highlightStripes | 0), emitRun);
  }
  for (let k = 0; k < o; k++) {
    const h = hs[HOutlines[k]];
    paintHighlight(out, axis, box, h, partOn(h.location, axis.system)!, 0, 1, config);
  }
  for (let c = 0; c < n; c++) HParts[c] = undefined;
}

/** The ends of the run `emitRun` paints that are no real ends (`Paint[5]`); 0 for a channel's whole ranges. */
let RunOpen = 0;

/** `walkStripes`' run: a channel's `[start, end)` in its stripe. A hidden range lights its break once, at its first run. */
function emitRun(
  channel: number,
  start: number,
  end: number,
  slot: number,
  count: number,
  first: number,
  last: number,
) {
  RunRanges[0] = start;
  RunRanges[1] = end;
  LitBreaks = first === 1;
  // a run's ends where the stripes change are no real ends: square there
  RunOpen = (first === 1 ? 0 : 1) | (last === 1 ? 0 : 2);
  paintHighlight(HOut, HAxis, HBox, HList[HIndex[channel]], RunPart, slot, count, HConfig);
  RunOpen = 0;
  LitBreaks = true;
}

/** `v` snapped to the device pixels of an axis. */
const snap = (v: number, dpr: number) => Math.round(v * dpr) / dpr;

/** One channel's part in stripe `slot` of `stripes` of the box, in its style; stripe edges on device pixels. */
function paintHighlight(
  out: PassWriter,
  axis: AxisTransform,
  box: Box,
  h: EmphasisHighlight,
  part: LocationPart,
  slot: number,
  stripes: number,
  config: RenderConfig,
) {
  const em = config.emphasis;
  const dpr = axis.dpr;
  Stripe.x = box.x;
  Stripe.width = box.width;
  Stripe.y = box.y;
  Stripe.height = box.height;
  if (stripes > 1) {
    // one edge between two stripes, so they neither overlap nor leave a seam
    const y0 = snap(box.y + (slot * box.height) / stripes, dpr);
    const y1 = snap(box.y + ((slot + 1) * box.height) / stripes, dpr);
    if (y1 <= y0) return;
    Stripe.y = y0;
    Stripe.height = y1 - y0;
  }
  const c = h.color;
  const style = h.style ?? 'tint';
  if (style === 'band') {
    const bottom = Stripe.y + Stripe.height;
    const bh = Math.min(Math.max(snap(em.highlightBand, dpr), 1 / dpr), Stripe.height);
    Stripe.y = bottom - bh;
    Stripe.height = bh;
    Paint[0] = c;
    Paint[1] = NoColor;
    Paint[2] = 0;
  } else if (style === 'outline') {
    Paint[0] = NoColor;
    Paint[1] = c;
    Paint[2] = config.outline.selection;
  } else {
    const a = ((c & 0xff) / 255) * em.highlightTint;
    Paint[0] = ((c & 0xffffff00) | Math.round(Math.max(0, Math.min(1, a)) * 255)) >>> 0;
    Paint[1] = NoColor;
    Paint[2] = 0;
  }
  // every style a squircle by its stripe's height (a `band` bar below the minimum: square), round where it really ends
  Paint[3] = 1;
  Paint[4] = 0;
  Paint[5] = RunOpen;
  paintRanges(out, axis, Stripe, part, config);
}

/** The point `hitView` places (reused: a hover makes no object). */
const Point = new PlacePoint();

/**
 * What is under a point (_spec/README.md "Render", Lanes and Hits): a break within `BreakSlop` px of its notch (above
 * a strip's band column, or a hidden items' tick a layer hits as the break), the section layout's own marks, else the
 * point placed (`Placement.at`: its strip, band and the point in the strip's frame) and the topmost layer's hit in that
 * band. On data the cells and items beside a break win: the hairline has no width of its own.
 */
/** The point a layer's hit reads (`LayerType.hit`): one array, so no number is boxed on the call. */
const HitAt = new Float64Array(2);
/** The page point a placement's `at` reads (`Placement.at`: from an array, never passed). */
const PageAt = new Float64Array(2);

/** `HitAt` set to the placed point (in its strip's frame). */
function pointAt(): Float64Array {
  HitAt[0] = Point.x;
  HitAt[1] = Point.y;
  return HitAt;
}

export function hitView(layout: ViewLayout, x: number, y: number, touch = false): Hit | undefined {
  for (const s of layout.sections) {
    const strips = s.placement.strips;
    for (let i = 0; i < strips.n; i++) {
      const col = s.columns[i];
      const breaks = s.stripBreaks[i];
      // a polar strip's notches are its layout's to hit (the seam: `hitChrome`, through the warp's inverse)
      if (!col || !breaks.length || strips.warp[i] !== WarpKind.translate) continue;
      // the page point passed as it came (no arithmetic on it here: a number computed is boxed on each call)
      const notch = breakAt(s, i, col, breaks, x, y);
      if (notch) return notch;
    }
    // the ruler rail's grab zone, over the ruler (an idle rail is no target: a press falls through)
    const rail = s.rail ? hitRail(s.rail, x, y, s.config, touch) : undefined;
    if (rail) return rail;
    if (s.rowsRails)
      for (let i = 0; i < s.rowsRails.length; i++) {
        const h = hitRail(s.rowsRails[i], x, y, s.config, touch);
        if (h) return h;
      }
    // what a layer draws outside its band (a rows layer's column titles, above the rows), on its first strip
    for (let j = 0; j < s.tracks.length; j++) {
      const lays = s.tracks[j].layers;
      for (let i = 0; i < lays.length; i++) {
        const t = lays[i].type;
        if (!t.hitOutside) continue;
        // in the first strip's frame (a text figure's first line)
        HitAt[0] = strips.n ? x - strips.dx[0] : x;
        HitAt[1] = strips.n ? y - strips.dy[0] : y;
        const h = t.hitOutside(lays[i].lay, HitAt);
        if (h) return h;
      }
    }
    const own = s.layout.hitChrome?.(s, x, y);
    if (own) return own;
    // a `section` layer's marks, on the page, over the strips'
    for (let i = s.sectionLays.length - 1; i >= 0; i--) {
      const l = s.sectionLays[i];
      HitAt[0] = x;
      HitAt[1] = y;
      const h = l.type.hit(l.lay, HitAt);
      if (h?.element) return h;
    }
    PageAt[0] = x;
    PageAt[1] = y;
    if (!s.placement.at(PageAt, Point)) continue;
    // a band's header in the label column (a rows layer's tree and labels): its layers of two axes hit it
    if (Point.kind === 'gutter' && Point.band >= 0) {
      const lays = s.lays[Point.strip * strips.bands + Point.band];
      for (let i = lays ? lays.length - 1 : -1; i >= 0; i--) {
        if (lays![i].type.axes !== 2) continue;
        const h = lays![i].type.hit(lays![i].lay, pointAt());
        if (h) return withBlock(s, Point.strip, h);
      }
    }
    if (Point.kind !== 'data') continue;
    // a layer drawn across other bands of its strip (a frame's outline, `across`) is hit there, over their layers
    for (let j = 0; j < strips.bands; j++) {
      const own = s.lays[Point.strip * strips.bands + j];
      for (let i = own ? own.length - 1 : -1; i >= 0; i--) {
        if (!Array.isArray(own![i].layer.params.across)) continue;
        const h = own![i].type.hit(own![i].lay, pointAt());
        if (h?.element) return h;
      }
    }
    const lays = s.lays[Point.strip * strips.bands + Point.band];
    if (!lays) return undefined;
    // the topmost layer with an element under the point, else the topmost with a position
    let position: Hit | undefined;
    for (let i = lays.length - 1; i >= 0; i--) {
      const l = lays[i];
      const h = l.type.hit(l.lay, pointAt());
      if (h?.element || h?.target === 'break') return withBlock(s, Point.strip, h);
      position ??= h;
    }
    return position && withBlock(s, Point.strip, position);
  }
  return undefined;
}

/**
 * A placement's window (a ring's), out of its key: set on a kept placement too. Active while something follows the
 * section or it follows (its window then its leader's range); shown when active or `rings.window: 'always'`, in a
 * figure only with `windows`. A function of the module, not a closure a section a frame.
 */
function setWindow(
  layout: SectionLayoutType,
  placement: Placement,
  at: string,
  viewport: Viewport | undefined,
  followed: Map<string, Viewport> | undefined,
  input: LayoutInput,
  config: RenderConfig,
): void {
  const own = followed?.get(at);
  let active = own !== undefined;
  const sections = input.view.sections;
  for (let i = 0; i < sections.length && !active; i++) if (sections[i].follow?.leader === at) active = true;
  const shown = (active || config.rings.window === 'always') && (!input.figure || input.windows === true);
  layout.window!(placement, own ?? viewport, active, shown);
}

/**
 * A hit with its block (its index in the section) when the section has blocks (_spec/README.md "Render", Lanes and
 * Hits): written into the layer's hit, no copy (a layer's hit is its lay's, and a lay is one strip's: one block's).
 */
function withBlock(s: SectionLayout, strip: number, h: Hit): Hit {
  if (s.section.blocks) {
    const st = s.placement.strips;
    (h as { block?: number }).block = st.block[strip];
    // a grid over blocks: the strip's y block too
    if (st.ys[strip]) (h as { yBlock?: number }).yBlock = st.yBlock[strip];
  }
  return h;
}

/** The break whose notch is within `BreakSlop` px of a point above the band column (in the strip's frame); else undefined. */
/** A notch's hit at a page point over strip `i`. */
function breakAt(
  s: SectionLayout,
  i: number,
  band: Box,
  b: readonly number[],
  px: number,
  py: number,
): Hit | undefined {
  const strips = s.placement.strips;
  const axis = strips.axes[i];
  const x = px - strips.dx[i];
  const y = py - strips.dy[i];
  if (y < band.y - Notch.height - BreakSlop || y >= band.y) return undefined;
  for (let k = 0; k + 2 < b.length; k += 3)
    if (Math.abs(x - b[k]) <= BreakSlop) return breakHit(s, axis, b[k], b[k + 1], b[k + 2]);
  return undefined;
}
