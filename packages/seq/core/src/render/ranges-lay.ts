/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { fontSize } from './config.js';
import { type ItemExtents, type ItemsBase, type ItemsPrep, newItemsBase } from './items.js';
import type { Channel } from '@molstar/seqx-core/theme/index';
import { Extents, type LayerChannels, RowDatum } from './channels.js';
import { LaneIndex, StickyLanes } from './lanes.js';
import { type GeneScratch, newGeneScratch } from './ranges-gene.js';
import type { LayoutContext } from './layer-type.js';

// What a ranges layer keeps frame to frame (_spec/README.md "Render", Layer types): its prepared form, its items in
// view and their lanes, its marks and hit marks, its label lane's labels, and the scratch its layout writes over.

/**
 * How an item is drawn: a block (0, the default), or a line (`line` param): solid, dashed, dotted, wavy, double, a box
 * a highlight, or a gene model (`gene`: `ranges-gene.ts`). The dash and the wave are anchored to the content: they start at the item's own start (a piece's),
 * so they move with it when the view pans. A style that is none of these is drawn solid and said once.
 */
export const LineStyles = [
  'blocks',
  'solid',
  'dashed',
  'dotted',
  'wavy',
  'double',
  'box',
  'highlight',
  'gene',
] as const;
export const Solid = 1;
/** The gene line's index in `LineStyles`. */
export const Gene = 8;

export interface RangesPrep extends ItemsPrep {
  /** The rows of every item when the data is inline (an item's ends are known), else undefined: no arrowhead. */
  readonly extents?: ItemExtents;
  /** The style asked for, as written: undefined for blocks. */
  readonly line?: string;
  /** An arrowhead on the item's `end`, `start`, or both: bits 1 and 2 (none: 0). */
  readonly arrow: number;
  /** Whether an item's pieces are joined across the positions between them (the `connector` param: `line` or `none`). */
  readonly connector: boolean;
  readonly cyclic: boolean;
  /** Its channels (`color`, `size`, `opacity`), compiled per layout's tokens. */
  readonly channels: LayerChannels;
}

/**
 * The part of the density ramp a strip in a band under a lane's height reaches: a thin band (a text line's underline)
 * reads as a band in the light end, not as a dark rule under the line.
 */
export const ThinStripShare = 0.3;

type Mode = 'items' | 'bars' | 'density';

export interface RangesLay extends ItemsBase {
  mode: Mode;
  /** The style a mark is drawn in (an index of `LineStyles`; in `bars` a run is a plain line or block) and the arrows. */
  style: number;
  arrow: number;
  /** The items holding the end and the start of their item in reading order (what an arrowhead goes on). */
  first: Uint8Array;
  last: Uint8Array;
  /**
   * The items in view (`items` of them; the arrays are kept longer, written over frame to frame): element key, lane id,
   * extent on the data's system, label, fill.
   */
  items: number;
  keys: (string | number)[];
  ids: string[];
  starts: number[];
  ends: number[];
  labels: (string | undefined)[];
  fills: number[];
  /** Each item's color group: its legend's order for a palette's value (listed first), else a distinct fill's first appearance. */
  group: number[];
  /** The groups' ids by value (a palette's value or a fill), and a palette value's name by id: kept, cleared each read. */
  groupIds: Map<string | number, number>;
  groupNames: string[];
  /** Each item's size channel: a share of its lane's height (1 when unthemed); its hit box stays the lane's. */
  sizes: number[];
  /** The channels this layout read, and the datum they read each item through. */
  color?: Channel;
  size?: Channel;
  datum: RowDatum;
  extents: Extents;
  lanes: Int32Array;
  laneCount: number;
  /** Items that did not fit the lanes: drawn in the density row (or the layer's density). */
  overflow: number;
  /** Marks: drawn extents (snapped), lane, item (−1 for a merged run), round ends, label width. */
  m: number;
  mx0: Float32Array;
  mx1: Float32Array;
  mlane: Int32Array;
  mitem: Int32Array;
  mround: Uint8Array;
  mlabel: Float32Array;
  mcut: Int32Array;
  /** Each mark's squircle radius (0: square) and its label's ink, made at layout: the draw computes no number. */
  mradius: Float64Array;
  mink: Uint32Array;
  /** Each mark's arrowheads (1: tip on the left, 2: on the right) and, for a line, where its content starts, px. */
  mtip: Uint8Array;
  manchor: Float64Array;
  /** Each mark's thinness (the gene line: 1 a thin box, outside the item's thick part). */
  mthin: Uint8Array;
  /**
   * The hit marks (`hn` of them): each item's own marks, before `bars` merges a lane's runs (a run of several items is
   * drawn as one mark, hit and outlined by the item under the pointer); the drawn marks at `items`. The arrays are the
   * drawn marks' or the merge's scratch, set at each layout.
   */
  hn: number;
  hx0: Float32Array;
  hx1: Float32Array;
  hlane: Int32Array;
  hitem: Int32Array;
  hround: Uint8Array;
  hradius: Float64Array;
  /** How far past a hit mark a point still hits it, px: a device pixel at `bars` (a run's gaps), else 0. */
  hslack: number;
  index: LaneIndex;
  /** `mergeRuns`'s scratch: the marks' order and their fields in it (and the radii of the marks), kept while large enough. */
  merge: {
    order: Int32Array;
    x0: Float32Array;
    x1: Float32Array;
    lane: Int32Array;
    item: Int32Array;
    round: Uint8Array;
    radius: Float64Array;
    tip: Uint8Array;
    /** A run's tally of its items by color group (zeroed after each run). */
    count: Float64Array;
  };
  /** The heat strip's band: the whole box, or the bottom row past the lanes. */
  densityY: number;
  densityH: number;
  sticky: StickyLanes;
  font: string;
  fontSize: number;
  radius: { ratio: number; max: number; min: number };
  /**
   * Each item's value of the color legend's field (text: a legend's `painted` entries), stored per item and read into a
   * set only when a legend is derived (on settle, never per frame).
   */
  painted: (string | undefined)[];
  /** A table's colors by dictionary code (kept, grown). */
  byCode: Uint32Array;
  /** Lane ids by key and occurrence, kept across frames (a layout makes no string an item); cleared when large. */
  idCache: Map<string | number, string[]>;
  /** The items by start, kept. */
  order: Int32Array;
  /** The most lanes the layer has shown, its whole data's not known (a tiled source's): what its lanes are centred by. */
  seen: number;
  /**
   * The band's lanes scroll inside it (_spec/README.md "Render", Track sizing): the lanes past the band, the whole lanes
   * it shows (`room`), the first lane in view in lanes (`top`, fractional) and the most it may be (`scrollMax`: 0 while
   * every lane fits); `laneY` is shifted by it, `laneY0` the unshifted top of the lane block.
   */
  scrolls: boolean;
  room: number;
  top: number;
  scrollMax: number;
  laneY0: number;
  /** The label lane's top (−1: none) and the labels placed there (`ln`): x, width, where their item starts, item. */
  labelY: number;
  ln: number;
  lx: Float32Array;
  lw: Float32Array;
  lwant: Float32Array;
  litem: Int32Array;
  /**
   * The connectors (`cn` of them; axis.md, An item split by a path): the hairline between two pieces of one item with
   * shown positions between them, its x extent clipped to the box, lane and item, with the line's width, ink and alpha
   * as the layout read them (`cwidth`, `cink`, `calpha`: 0 alpha draws none).
   */
  cn: number;
  cx0: Float32Array;
  cx1: Float32Array;
  clane: Int32Array;
  citem: Int32Array;
  /** Each connector's chevrons (the gene line's strand on the screen): 0 none, 1 pointing right, 2 left. */
  cdir: Uint8Array;
  /** Where each gene line's hairline starts on the screen, unclipped: its chevrons run from it, so a pan moves them. */
  canchor: Float64Array;
  cwidth: number;
  cink: number;
  calpha: number;
  /** The gene line's scratch: the items' rows grouped, their extents, the thick part's columns. */
  gene: GeneScratch;
  /** Scratch: an item's leftmost mark's x, the candidates by priority. */
  left: Float32Array;
  right: Float32Array;
  cand: Int32Array;
  /** Scratch: the keys of the items labeled inside a mark (kept, cleared each layout). */
  labeled: Set<string | number>;
  /** Scratch: where the last item of a label's run ends, px (kept, cleared each layout). */
  runEnds: Map<string, number>;
}

/** An item's lane id: its key and which occurrence of it in a tile, made once and kept. */
export function laneId(cache: Map<string | number, string[]>, key: string | number, k: number): string {
  let ids = cache.get(key);
  if (!ids) cache.set(key, (ids = []));
  return (ids[k] ??= `${typeof key === 'number' ? 'n' : 's'}${key}#${k}`);
}

export function newLay(ctx: LayoutContext, prep: RangesPrep): RangesLay {
  const base = newItemsBase(ctx, prep);
  return {
    ...base,
    mode: 'items',
    items: 0,
    // the rows' own arrays: one reading of the tiles for every items layer
    keys: base.rows.keys,
    ids: [],
    starts: base.rows.starts,
    ends: base.rows.ends,
    labels: base.rows.labels,
    fills: [],
    group: [],
    groupIds: new Map(),
    groupNames: [],
    sizes: [],
    datum: new RowDatum(),
    extents: new Extents(),
    lanes: new Int32Array(64),
    laneCount: 0,
    overflow: 0,
    m: 0,
    mx0: new Float32Array(64),
    mx1: new Float32Array(64),
    mlane: new Int32Array(64),
    mitem: new Int32Array(64),
    mround: new Uint8Array(64),
    mlabel: new Float32Array(64),
    mcut: new Int32Array(64),
    mradius: new Float64Array(64),
    mink: new Uint32Array(64),
    mtip: new Uint8Array(64),
    manchor: new Float64Array(64),
    mthin: new Uint8Array(64),
    style: 0,
    arrow: 0,
    first: new Uint8Array(64),
    last: new Uint8Array(64),
    hn: 0,
    hx0: new Float32Array(0),
    hx1: new Float32Array(0),
    hlane: new Int32Array(0),
    hitem: new Int32Array(0),
    hround: new Uint8Array(0),
    hradius: new Float64Array(0),
    hslack: 0,
    index: new LaneIndex(),
    merge: {
      order: new Int32Array(0),
      x0: new Float32Array(0),
      x1: new Float32Array(0),
      lane: new Int32Array(0),
      item: new Int32Array(0),
      round: new Uint8Array(0),
      radius: new Float64Array(0),
      tip: new Uint8Array(0),
      count: new Float64Array(0),
    },
    densityY: 0,
    densityH: 0,
    sticky: new StickyLanes(),
    font: ctx.config.itemFont,
    fontSize: fontSize(ctx.config.itemFont),
    radius: { ratio: ctx.config.squircleRatio, max: ctx.config.squircleMax, min: ctx.config.squircleMin },
    painted: [],
    byCode: new Uint32Array(16),
    idCache: new Map(),
    order: new Int32Array(64),
    seen: 0,
    scrolls: false,
    room: 0,
    top: 0,
    scrollMax: 0,
    laneY0: 0,
    cn: 0,
    cx0: new Float32Array(8),
    cx1: new Float32Array(8),
    clane: new Int32Array(8),
    citem: new Int32Array(8),
    cdir: new Uint8Array(8),
    canchor: new Float64Array(8),
    cwidth: 1,
    cink: 0,
    calpha: 0,
    labelY: -1,
    ln: 0,
    lx: new Float32Array(16),
    lw: new Float32Array(16),
    lwant: new Float32Array(16),
    litem: new Int32Array(16),
    gene: newGeneScratch(),
    left: new Float32Array(64),
    right: new Float32Array(64),
    cand: new Int32Array(64),
    labeled: new Set(),
    runEnds: new Map(),
  };
}

/** An item's piece on its own system, and its pieces through a map: kept (a layout runs to its end). */
export const itemPiece: number[] = [0, 0];
export const mappedPieces: number[] = [];
