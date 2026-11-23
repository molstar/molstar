/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import type { Want } from '@molstar/seqx-core/data';
import type { OverlayNode, StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer, ResolvedSection, ResolvedTrack, ResolvedView } from '@molstar/seqx-core/view';
import { AxisTransform, type Viewport } from './axis.js';
import type { RenderConfig } from './config.js';
import type { PassName, PassWriter } from './draw-lists.js';
import type { Box, Emphasis, Hit, LayerType, Overlay, StoreReader } from './layer-type.js';
import type { ThemeContext } from './channels.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';
import { LabelLane, type LayerContent, type TrackSizing, type TrackSizingHook } from './sizing.js';
import type { Placement, SectionLayoutType } from './layout/types.js';
import { Rail, type RowsZoom } from './rail.js';
import type { ImageCache } from './images.js';
import type { NarrowColumns } from './narrow.js';
import type { SectionCache } from './tracks-cache.js';
import type { HeaderLine, HeaderPartsOf, SectionHeaderLine } from './headers.js';

// The tracks layout's records (_spec/README.md "Render", Tracks layout): what it is given (`LayoutInput`), what it
// gives (the view's, sections', tracks' and layers' layouts) and what it keeps between frames (`LayoutState`).

export interface LayerLayout {
  readonly key: string;
  readonly layer: ResolvedLayer;
  readonly type: LayerType<any, any>;
  readonly lay: unknown;
  /**
   * The layer's wants together were refused (past `maxTilesPerLayer` or `maxBytesPerLayer`): its band hatched with
   * "Zoom in to read: ⌘/Ctrl+scroll, or +", which no wait resolves (pending does).
   */
  readonly zoomIn?: boolean;
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
   * Its header's line in the label column (_spec/README.md "Render", Track headers): its slots, its states (the "n
   * filtered" chip among them, `seq:filtered`), the "⋯"; in a layout that draws a label column of track headers.
   */
  readonly header?: HeaderLine;
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
  /**
   * Each track with no items (`track.empty`): its one-line row, with a note (`note`) or without (`collapse`); undefined
   * for any other. The sizing's kept list (a layout makes none).
   */
  readonly empties: readonly ('note' | 'collapse' | undefined)[];
  /** Its header line above its first track (`section.header`), when drawn: its hidden tracks' chip, its slots. */
  readonly header?: SectionHeaderLine;
  /** The breaks in view, `[x, hiddenStart, hiddenEnd, …]` (x unsnapped): a hairline each across the band. */
  readonly breaks: readonly number[];
  /** Its ruler rail, when it draws one (_spec/README.md "Render", Navigators): the first of `rails`. */
  readonly rail?: Rail;
  /** Its positions rails, one a place of `navigator.railPlace` (each on the section's one window). */
  readonly rails?: readonly Rail[];
  /** Its rows layers' rails (upright, at the right edge of each band), when they draw them. */
  readonly rowsRails?: readonly Rail[];
  /** The rails of its tracks whose lanes scroll inside their band (upright, plain: the wheel moves them). */
  readonly lanesRails?: readonly Rail[];
  /**
   * Its own rail while its content scrolls inside its box (a text's lines, a tracks section's tracks below its ruler:
   * _spec/README.md "Render", Section scroll), upright at the box's right edge.
   */
  readonly sectionRail?: Rail;
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

/**
 * A bridge laid out (_spec/README.md "Render", Bridges; the bridge chunk's): its band between its sides, its layers'
 * lays, its draw of a pass and its hit (a page point in `at`).
 */
export interface BridgeLayout {
  /** Its index in `ResolvedView.bridges` (a hit's `bridge`), its address, id and label, its sides' addresses. */
  readonly index: number;
  readonly address: string;
  readonly id: string;
  readonly label: string;
  readonly a: string;
  readonly b: string;
  readonly box: Box;
  readonly lays: readonly LayerLayout[];
  /** The ribbons drawn (its ARIA description's count). */
  readonly count: number;
  draw(pass: PassName, out: PassWriter, emphasis: Emphasis | undefined): void;
  hit(at: Float64Array): Hit | undefined;
  /**
   * Each ribbon drawn, in drawing order (the topmost last), its ends on the bridge's sides: what the viewer's align,
   * select across and cursor read (_spec/README.md "Render", Bridges). Made on each call: a command's, never a frame's.
   */
  ribbons(): BridgeRibbon[];
  /** The topmost ribbon drawn of a pair's key (a hit's element): one built, never the list. */
  ribbon(key: string | number): BridgeRibbon | undefined;
}

/**
 * A ribbon drawn in a bridge, by its sides: on the bridge's `a` side its system, its range of positions and its drawn
 * extent in page px (`ax0 ≤ ax1`), the same on its `b` side; `reversed` when it is drawn crossed (its strand); its
 * pair's key; `length` its aligned positions (its narrower end's).
 */
export interface BridgeRibbon {
  readonly aSystem: string;
  readonly aStart: number;
  readonly aEnd: number;
  readonly ax0: number;
  readonly ax1: number;
  readonly bSystem: string;
  readonly bStart: number;
  readonly bEnd: number;
  readonly bx0: number;
  readonly bx1: number;
  readonly reversed: boolean;
  readonly key: string | number;
  readonly length: number;
  /** A stub's unplaced end (its side shows no place of it). */
  readonly stub?: 'a' | 'b';
  /** A bin pair's level and its count (absent for an item). */
  readonly bin?: number;
  readonly count?: number;
}

/**
 * A view's bridges laid out (the bridge chunk's): each one, and what the layout's draw, hit and export ask of them all,
 * so the eager layout holds one call each.
 */
export interface BridgesLayout {
  readonly list: readonly BridgeLayout[];
  draw(pass: PassName, out: PassWriter, emphasis: Emphasis | undefined): void;
  /** The hit at a page point (`at`): undefined off every band, null in a band over no ribbon. */
  hit(at: Float64Array): Hit | null | undefined;
  /** The bridge whose band holds a page point (its index in `ResolvedView.bridges`), −1 for none. */
  at(x: number, y: number): number;
  /** The export's document with each bridge's group (`data-bridge`, `data-a`, `data-b`, …) after the sections'. */
  svg(document: string): string;
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
  /** The columns a narrow frame laid out with, the plot's width and the frame's least (a view with a tracks section). */
  readonly narrow?: NarrowColumns;
  /** The frame's foot the pinned rails take, px (`pinned`): past the last section, kept free of the fill band. */
  readonly pinnedPx?: number;
  /** Its placed bridges laid out (absent: none, or the bridge chunk not landed yet). */
  readonly bridges?: BridgesLayout | undefined;
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
   * Each header's parts from its facts (_spec/README.md "Render", Track headers: a host's entries made into slots and
   * states), once a layout a header; `builtInHeaderStates` when absent. Part of the sections' key by identity.
   */
  readonly headers?: HeaderPartsOf;
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
  /** The host's rows zoom (the rows rails' bounds and grips): part of the key by identity; absent: the defaults. */
  readonly rowsZoom?: RowsZoom;
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
  /** A section of the last layout scrolled past its cap: the right gutter keeps a rail's slot for it. */
  sectionSlot = false;
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
export interface SizedHeights {
  readonly tracks: readonly ResolvedTrack[];
  readonly fit: number;
  readonly config: RenderConfig;
  readonly layout: SectionLayoutType<any>;
  readonly version: number;
  readonly sizings: LayoutInput['sizings'];
  readonly hook: LayoutInput['trackSizing'];
  /** The count of source lanes grown it was sized at (`knownLanesVersion`). */
  readonly known: number;
  readonly heights: readonly number[];
  /** Each track's content (one lane block a track), undefined for a track not sized by its content. */
  readonly contents: readonly (LayerContent | undefined)[];
  /** Each track with no items: how it is drawn (`note` or `collapse`), undefined for any other. */
  readonly empties: readonly ('note' | 'collapse' | undefined)[];
  readonly issue?: Issue;
}
