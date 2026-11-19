/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import type { Coverage, DataEntry, Want } from '@molstar/seqx-core/data';
import type { ElementKey, Elements, Location } from '@molstar/seqx-core/query';
import type { Overlay, StaticModel } from '@molstar/seqx-core/resolve';
import type { Legend } from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer, ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import type { ImageCache } from './images.js';
import type { ThemeContext } from './channels.js';
import type { RenderConfig } from './config.js';
import type { PassName, PassWriter } from './draw-lists.js';
import type { Placement, WindowPart } from './layout/types.js';
import type { TextMetrics } from './text.js';
import type { SeqTokens } from './tokens.js';

// The layer-type contract (_spec/README.md "Render", Layer types): a handful of members; built-ins are written to it
// as extensions will be.

/** What a layout reads of the data store: `DataStore` is one. */
export interface StoreReader {
  coverage(source: string): Coverage | undefined;
  keysOf(want: Want): string[];
  /** Whether a want is past the store's bound on its tiles (\`keysOf\` then lists none), where the reader says it. */
  tooMany?(want: Want): boolean;
  get(key: string): DataEntry | undefined;
  /** The store's tiling and its bound on a want's tiles (a layer picks a level within it), where it has them. */
  readonly settings?: { readonly tileBins: number; readonly maxTilesPerWant: number };
}

/** The dynamic nodes of tier 2 (the resolve pass's overlay): its version is part of the draw lists' key. */
export type { Overlay } from '@molstar/seqx-core/resolve';
export { EmptyOverlay } from '@molstar/seqx-core/resolve';

export interface Box {
  readonly x: number;
  readonly y: number;
  readonly width: number;
  readonly height: number;
}

export interface LayoutContext {
  /** The track's data box, CSS px. */
  readonly box: Box;
  /** The section's positions to x. */
  readonly axis: AxisTransform;
  readonly section: ResolvedSection;
  /** The layer's key: its address, else its track's address and `:` its index. */
  readonly key: string;
  readonly dpr: number;
  readonly config: RenderConfig;
  readonly tokens: SeqTokens;
  readonly text: TextMetrics;
  /** The view has settled (no data pending, no gesture): open domains are computed, sparse lanes repacked. */
  readonly settled: boolean;
  /** What the themes read: schemes, measures, scales, held domains, the reader's substitution. */
  readonly theme: ThemeContext;
  /** Asks for a region the layer reads (the store's want, by priority). */
  want(w: Want): void;
  /** Says what went wrong drawing the layer (a color it cannot read, data it cannot place). */
  report(issue: Issue): void;
  /** The strip laid out (_spec/README.md "Render", Section layouts): `box` and `axis` are in its frame. */
  readonly strip?: number;
  /**
   * The view's width, CSS px, in the page's frame (a strip's `box` is in the strip's own, translated by its `dx`: read
   * with that in mind): where a hover's text stops (a rows layer's cell read whole). The tracks layout's strips are not
   * translated, which is where a rows layer is laid out.
   */
  readonly viewWidth?: number;
  /** The union of the strip's bands (a layer with `across`, frames, outlines it). */
  readonly span?: Box;
  /** The section's placement, for a `section` layer (it draws in page coordinates). */
  readonly place?: Placement;
  /** Item labels drawn `inside` marks, or a sink a layout collects them into (rings: placed outside). */
  readonly labels?: 'inside' | LabelSink;
  /** A band's own scroll (display rows below its pinned band) and row zoom (px a row), when it scrolls its rows. */
  readonly bandTop?: number;
  readonly rowScale?: number;
  /** The clades a rows layer folds (session state, by name), when the band says; its `collapsed` when absent. */
  readonly folds?: readonly string[];
  /** The shared image cache a layer keeps its image tiles in (the renderer's). */
  readonly images?: ImageCache;
  /**
   * A rows layer's row columns as the view grants them (_spec/README.md "Render", Row columns): each column's width
   * after its yielding to a narrow view (0: hidden), in the order written, and the x its right columns start at.
   */
  readonly rowColumns?: { readonly widths: readonly number[]; readonly rightX: number };
  /**
   * The strip's second transform (`strip.y`: a grid's y, down the strip's frame), on a strip that has one
   * (_spec/README.md "Render", The grid layout).
   */
  readonly axisY?: AxisTransform;
  /** A `section` layer's band: its track's index in the section's bands (where its arcs stand). */
  readonly band?: number;
  /**
   * The strip is a chip (a one-position block, _spec/README.md "Render", Blocks): its label, which a letters layer draws
   * as a pill over its one position.
   */
  readonly chip?: string;
}

/** Where a layer gives its items' labels when its layout places them (rings: outside, with leaders). */
export interface LabelSink {
  /**
   * An item's label: its text, its extent `[start, end)` as positions on the strip's system (past each other over a
   * cyclic origin), its element, the hit its placed label gives (the item's, never a second hit at an anchor, EC-108)
   * and its order (what a layout that cannot place every label keeps first: lower first).
   */
  label(text: string, start: number, end: number, element: ElementKey, hit?: Hit, order?: number): void;
}

/**
 * What is under a point: a layer, its section (the section's id), and the element (an item's key, with its extent when
 * the layer knows it) or the position there, on the section's system.
 */
export interface Hit {
  readonly layer: string;
  readonly section: string;
  readonly element?: ElementKey;
  /** The element's extent, on its data's system. */
  readonly extent?: Location;
  /**
   * The element's hull on the section's system, positions `[start, end)`: its extent carried along the layer's path
   * (all its pieces), what a focus on it shows.
   */
  readonly hull?: readonly [number, number];
  readonly position?: number;
  readonly system: string;
  /**
   * A position on the second axis (_spec/README.md "Render", Lanes and Hits): a row of the rows system, a grid's y.
   * Flat, beside `position`, so a hit of two axes allocates nothing.
   */
  readonly ySystem?: string;
  readonly yPosition?: number;
  /**
   * `break`: a collapsed break (within 3 CSS px of its notch, or on a tick of the items it hides): `position` the
   * section's position at its edge (the last shown before it in display order), `extent` the hidden stretch on the
   * section's system, `layer` the section's address, `count` the hidden items a tick stands for.
   */
  readonly target?: HitTarget;
  readonly count?: number;
  /** A `window` hit's part: its body (dragged: pans) or an end (dragged: resizes about the other end). */
  readonly part?: WindowPart;
  /**
   * A navigator's `window` hit: the navigator's id, `rail` for a section's ruler rail (_spec/README.md "Render",
   * Navigators); `position` is then the position under the point at the navigator's scale.
   */
  readonly navigator?: string;
  /** The block (its index in the section) when the section has blocks. */
  readonly block?: number;
  /** A grid over blocks: the cell's y block (its index in the section; `block` is its x block). */
  readonly yBlock?: number;
  /**
   * What is under the point within the position, when a layer says: a logo's letter; a rows layer's column cell, what
   * it names (`title: value`, on its row's `row` hit); a column title's field (`column-title`).
   */
  readonly state?: string;
}

/**
 * What a hit is, when it is more than an element or a position (_spec/README.md "Render", Lanes and Hits): `break`,
 * `filter` (a track's chip), `window` (a ring's), `row` (a gutter's row, a rows header), `cell` (a cell of two axes:
 * `ySystem` and `yPosition` set), `tree-node` (a tree's node: its element the clade), `summary-row` (a folded
 * clade's row: its element the clade, `position` the column), `column-title` (a rows layer's column title: `state`
 * its field; a click orders the rows by it), `columns-hidden` (a rows layer's hidden-columns chip: `state` its text;
 * a hover names the columns, a click does nothing).
 */
export type HitTarget =
  | 'break'
  | 'filter'
  | 'window'
  | 'row'
  | 'cell'
  | 'tree-node'
  | 'summary-row'
  | 'column-title'
  | 'columns-hidden';

/** A section's id: the last step of its address. */
export function sectionId(section: { readonly address: string }): string {
  const a = section.address;
  return a.slice(a.lastIndexOf('/') + 1);
}

/** A highlight channel as drawn: its Location (on every system it shows on) and its color, tinted. */
export interface EmphasisHighlight {
  readonly location: Location;
  /** `0xRRGGBBAA`: drawn at `emphasis.highlightTint` opacity. */
  readonly color: number;
  /** How it is drawn in its stripe: a `tint` (the default), a `band` along the foot, an `outline`. */
  readonly style?: HighlightStyle;
}

/** A highlight channel's look (_spec/README.md "Render", Emphasis). */
export const HighlightStyles = ['tint', 'band', 'outline'] as const;
export type HighlightStyle = (typeof HighlightStyles)[number];

/**
 * What the emphasis passes draw (_spec/README.md "Render", Emphasis): the hover, the selection (a Location), the
 * keyboard cursor, the highlight channels, and the hover and the selection carried onto linked sections.
 */
export interface Emphasis {
  readonly hover?: Hit;
  readonly selection?: Location;
  /**
   * The selection's element keys: a layer draws those its Location cannot hold (a boundary, an item of zero length, as a
   * hairline, EC-127).
   */
  readonly elements?: Elements;
  /**
   * The keyboard cursor; on a cell of two axes (a rows cell, a grid's) its row too, drawn by the layer of two axes on
   * its own section (`section`, its address: the others mark its column); on rows its `part`, a cell or the header.
   */
  readonly cursor?: {
    readonly system: string;
    readonly position: number;
    readonly section?: string;
    readonly ySystem?: string;
    readonly yPosition?: number;
    readonly part?: 'cell' | 'header';
  };
  /** The highlight channels in the order they were set: each a tint of its color (the `highlight` pass). */
  readonly highlights?: readonly EmphasisHighlight[];
  /** The window being dragged (its hit): a rail is drawn grown while its drag runs, wherever the pointer goes. */
  readonly drag?: Hit;
  /**
   * The hover and the selection carried onto the linked sections' systems (`projectOntoView`), one Location each, keyed
   * by system: drawn in the same look, lighter (`emphasis.projected`). The projected hover is drawn on every section
   * but `hoverSection` (its layer outlines it there); the projected selection where the selection has no part.
   */
  readonly projected?: {
    readonly hover?: Location;
    readonly selection?: Location;
    /** The hover's own section's address: the projected hover is drawn on every other section. */
    readonly hoverSection?: string;
  };
}

/** An item a layer shows, as an export names it and the viewer announces it (`LayerType.items`). */
export interface ItemDescription {
  readonly key: string | number;
  /** Its label (its key when it has none). */
  readonly label: string;
  /** Its extent on its data's system, `[start, end)`. */
  readonly start: number;
  readonly end: number;
  /** Where it is drawn, in the section's x (clamped to the data box) and the layer's band. */
  readonly x0: number;
  readonly x1: number;
  /** Its band on the page when it is not the track's whole (a rows layer's row): top and height. */
  readonly y?: number;
  readonly height?: number;
  /** The system its extent is on, when it is not the layer's data's (a rows layer's rows system). */
  readonly system?: string;
}

/**
 * A layer type: `prepare` is tier 1 (cached per resolved layer), `layout` per viewport (given its last layout to reuse
 * its buffers and keep its lanes), `draw` one pass into the lists, `hit` from the same layout.
 */

export interface LayerType<Prep = unknown, Lay = unknown> {
  readonly name: string;
  /** The layer type a target without this one draws. */
  readonly fallback?: string;
  prepare(layer: ResolvedLayer, model: StaticModel): Prep;
  layout(prep: Prep, ctx: LayoutContext, data: StoreReader, overlay: Overlay, prev?: Lay): Lay;
  draw(lay: Lay, pass: PassName, out: PassWriter, emphasis?: Emphasis): void;
  /**
   * The mark under a point, `at[0]` and `at[1]` (in the strip's frame; a `section` layer's on the page). The point comes
   * in an array: two fractional numbers passed to a call that is not inlined are boxed, and a hit allocates nothing.
   */
  hit(lay: Lay, at: Float64Array): Hit | undefined;
  /**
   * What it draws outside its band that a point hits (a rows layer's column titles, in the band above the rows), on the
   * page; asked before the bands' layers, so it is quick to refuse a point not on it.
   */
  hitOutside?(lay: Lay, at: Float64Array): Hit | undefined;
  /** The legends of its themed channels as last laid out (_spec/README.md "Themes", Legends), titled by its track. */
  legend?(lay: Lay, options: { readonly title: string }): readonly Legend[];
  /** The y domain it drew (a values layer's), recorded by an export. */
  domain?(lay: Lay): readonly [number, number] | undefined;
  /**
   * The items it shows, as the accessible layer of an export and the viewer's announcements name them (`ItemDescription`):
   * each item once, its pieces' hull. Absent: the layer has no items to name (letters, values).
   */
  items?(lay: Lay): readonly ItemDescription[];
  /** The rows its data's filter left out of what it read (the track's "n filtered" chip). */
  filtered?(lay: Lay): number;
  /**
   * The height it used in its band on this strip, 0 when it drew nothing there (text collapses a band no layer used);
   * the band's height when absent.
   */
  used?(lay: Lay): number;
  /**
   * `strip` (the default): laid out once a strip, drawn in its frame; `section`: once a section with the placement
   * (`ctx.place`), drawn in page coordinates (what crosses strips).
   */
  readonly scope?: 'strip' | 'section';
  /**
   * 2: it has a second axis of its own (a grid's y, rows): the emphasis passes call its `draw` with the emphasis, whose
   * parts on its y systems it paints through its own y transform; the strip paints x.
   */
  readonly axes?: 1 | 2;
}
