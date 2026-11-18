/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Issue } from '@molstar/core/repr/tree';
import type { Coverage, DataEntry, Want } from '@molstar/seqx-core/data';
import type { ElementKey, Elements, Location } from '@molstar/seqx-core/query';
import type { Overlay, StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer, ResolvedSection } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import type { RenderConfig } from './config.js';
import type { PassName, PassWriter } from './draw-lists.js';
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
  /** Asks for a region the layer reads (the store's want, by priority). */
  want(w: Want): void;
  /** Says what went wrong drawing the layer (a color it cannot read, data it cannot place). */
  report(issue: Issue): void;
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
  readonly position?: number;
  readonly system: string;
}

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

/** A layer's legend: a ramp over a domain, or categories with their colors. */
export interface LegendSpec {
  readonly kind: 'ramp' | 'categories';
  readonly domain?: readonly [number, number];
  readonly labels?: readonly string[];
  readonly colors: readonly number[];
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
  legend?(lay: Lay): LegendSpec | undefined;
}
