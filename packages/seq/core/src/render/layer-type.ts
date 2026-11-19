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
  get(key: string): DataEntry | undefined;
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
  /** The union of the strip's bands (a layer with `across`, frames, outlines it). */
  readonly span?: Box;
  /** The section's placement, for a `section` layer (it draws in page coordinates). */
  readonly place?: Placement;
  /** Item labels drawn `inside` marks, or a sink a layout collects them into (rings: placed outside). */
  readonly labels?: 'inside' | LabelSink;
  /** A band's own vertical scroll and row zoom (px a row), when it scrolls its rows. */
  readonly bandTop?: number;
  readonly rowScale?: number;
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
  /** The block (its index in the section) when the section has blocks. */
  readonly block?: number;
}

/**
 * What a hit is, when it is more than an element or a position (_spec/README.md "Render", Lanes and Hits): `break`,
 * `filter` (a track's chip), `window` (a ring's), `row` (a gutter's row, a rows header), `cell` (a cell of two axes:
 * `ySystem` and `yPosition` set), `tree-node` (a tree's node: its element the clade), `summary-row` (a folded
 * clade's row: its element the clade, `position` the column).
 */
export type HitTarget = 'break' | 'filter' | 'window' | 'row' | 'cell' | 'tree-node' | 'summary-row';

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
}

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
  readonly cursor?: { readonly system: string; readonly position: number };
  /** The highlight channels in the order they were set: each a tint of its color (the `highlight` pass). */
  readonly highlights?: readonly EmphasisHighlight[];
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
  hit(lay: Lay, x: number, y: number): Hit | undefined;
  /** The legends of its themed channels as last laid out (_spec/README.md "Themes", Legends), titled by its track. */
  legend?(lay: Lay, options: { readonly title: string }): readonly Legend[];
  /** The y domain it drew (a values layer's), recorded by an export. */
  domain?(lay: Lay): readonly [number, number] | undefined;
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
