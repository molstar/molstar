/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArrayRange } from '@molstar/core/data/util/sort';
import { CategoryColumn, issueOf } from '@molstar/core/repr/tree';
import type { RangesDensity } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { type AxisTransform, Cut } from './axis.js';
import { fontSize } from './config.js';
import {
  Anchor,
  drawPass,
  grown,
  type ImageData,
  NoColor,
  type PassName,
  type PassWriter,
  passDraws,
  Round,
} from './draw-lists.js';
import {
  beginItems,
  drawTicks,
  type ItemExtents,
  type ItemRows,
  type ItemsBase,
  type ItemsPrep,
  itemExtents,
  itemHit,
  markEnds,
  newItemsBase,
  RowImages,
  readItems,
  tickHit,
  tickOn,
  itemsUsed,
  describeItems,
} from './items.js';
import { RenderIssues } from './issues.js';
import type { Channel } from '@molstar/seqx-core/theme/index';
import { channelForms, channelLegends, Extents, type LayerChannels, RowDatum } from './channels.js';
import { LaneIndex, StickyLanes } from './lanes.js';
import { type LayerSource, layerSource, toSection } from './layer-data.js';
import type { Emphasis, Hit, LabelSink, LayerType, LayoutContext, Overlay, StoreReader } from './layer-type.js';
import { arrowArgs, drawFailed, drawSpans } from './marks.js';
import { squircleRadius } from './squircle.js';
import { inkOn, mixColors, type SeqTokens, withAlpha } from './tokens.js';

// the anchor a number read once: a draw reads no property of an imported object (that site's feedback, cleared, made
// a warm draw deoptimize and allocate without end once the 150,000-piece scenario had run first)
const AnchorStart: number = Anchor.start;

// The `ranges` layer type (_spec/README.md "Render", Layer types): items in sticky lanes. Each item is carried onto
// the section along the layer's path (a domain on its gene: a piece an exon) and drawn as a squircle where it is big
// enough, square below, with square ends where the view or a piece's edge cuts it (the frame rule). Below a device
// pixel a position, the marks of a lane merge into one bar a run; at the density step the layer is one heat strip, a
// device column a pixel. Items past the lanes collapse into a density row at the bottom (`ranges.overflow`:
// `extra-lanes`), or the whole layer does (`layer`). An item's key is its `id`, else its `row`, else its tile and
// index; a key and start seen in an earlier tile is not drawn again (a range repeated in the tiles it touches), and a
// key repeated with another start is another piece of one item (a wrapped row, the exons of a transcript).

/**
 * How an item is drawn: a block (0, the default), or a line (`line` param): solid, dashed, dotted, wavy, double, a box
 * or a highlight. The dash and the wave are anchored to the content: they start at the item's own start (a piece's),
 * so they move with it when the view pans. A style that is none of these is drawn solid and said once.
 */
export const LineStyles = ['blocks', 'solid', 'dashed', 'dotted', 'wavy', 'double', 'box', 'highlight'] as const;
const Solid = 1;

/** The dash of a dashed and a dotted line (on, then as long off), the wave's half period and amplitude, CSS px. */
const Dash = 6;
const Dot = 2;
const WaveStep = 3;
const WaveAmplitude = 2;
const LineWidth = 2;
/** A mark this wide (CSS px) or wider draws its item's arrowhead, at any step: narrower, the head would be all of it. */
const ArrowMinPx = 4;
/** A highlight's alpha: 0.35 of 255, the channel's default. */
const HighlightAlpha = 89;

export interface RangesPrep extends ItemsPrep {
  /** The rows of every item when the data is inline (an item's ends are known), else undefined: no arrowhead. */
  readonly extents?: ItemExtents;
  /** The style asked for, as written: undefined for blocks. */
  readonly line?: string;
  /** An arrowhead on the item's `end`, `start`, or both: bits 1 and 2 (none: 0). */
  readonly arrow: number;
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
  /** Each mark's squircle radius (0: square) and its label's ink, made at layout: the draw computes no number. */
  mradius: Float64Array;
  mink: Uint32Array;
  /** Each mark's arrowheads (1: tip on the left, 2: on the right) and, for a line, where its content starts, px. */
  mtip: Uint8Array;
  manchor: Float64Array;
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
}

/** An item's lane id: its key and which occurrence of it in a tile, made once and kept. */
function laneId(cache: Map<string | number, string[]>, key: string | number, k: number): string {
  let ids = cache.get(key);
  if (!ids) cache.set(key, (ids = []));
  return (ids[k] ??= `${typeof key === 'number' ? 'n' : 's'}${key}#${k}`);
}

function newLay(ctx: LayoutContext, prep: RangesPrep): RangesLay {
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
    mradius: new Float64Array(64),
    mink: new Uint32Array(64),
    mtip: new Uint8Array(64),
    manchor: new Float64Array(64),
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
  };
}

/**
 * The next mark's index, its arrays grown to hold it and its label, radius and ink cleared; the caller writes its extent,
 * lane, item and ends (a number passed to a call that is not inlined is boxed: an allocation a mark).
 */
function nextMark(lay: RangesLay): number {
  const j = lay.m++;
  if (j >= lay.mx0.length) {
    lay.mx0 = grown(lay.mx0, j + 1);
    lay.mx1 = grown(lay.mx1, j + 1);
    lay.mlane = grown(lay.mlane, j + 1);
    lay.mitem = grown(lay.mitem, j + 1);
    lay.mround = grown(lay.mround, j + 1);
    lay.mlabel = grown(lay.mlabel, j + 1);
    lay.mradius = grown(lay.mradius, j + 1);
    lay.mink = grown(lay.mink, j + 1);
    lay.mtip = grown(lay.mtip, j + 1);
    lay.manchor = grown(lay.manchor, j + 1);
  }
  lay.mlabel[j] = 0;
  lay.mradius[j] = 0;
  lay.mink[j] = 0;
  lay.mtip[j] = 0;
  return j;
}

const images = new RowImages();

/** A wave's scratch points (x, y pairs), kept and grown at layout. */
let wave = new Float64Array(64);
const segment = new Float64Array(4);

/** The style `line` names (an index of `LineStyles`; 0: blocks); another value is solid, said once. */
function lineStyle(line: string | undefined, ctx: LayoutContext): number {
  if (line === undefined) return 0;
  const i = LineStyles.indexOf(line as (typeof LineStyles)[number]);
  if (i > 0) return i;
  ctx.report(
    issueOf(RenderIssues, 'seq.render.line-style', `The line style ${JSON.stringify(line)} is not drawn`, {
      address: ctx.key,
      hint: RenderIssues.codes['seq.render.line-style'].hint,
    }),
  );
  return Solid;
}

/**
 * Where the content of piece `a` of the axis starts on the screen, unclipped: the left edge of its low boundary on a
 * forward axis, the right edge on a reversed one (the dash and the wave run from it).
 */
function anchorOf(axis: AxisTransform, a: number, s: number, e: number): number {
  const ua = axis.xOfUnit(axis.unitOf(a, Math.max(s, axis.pieces[2 * a])));
  const ub = axis.xOfUnit(axis.unitOf(a, Math.min(e, axis.pieces[2 * a + 1])));
  return axis.dir === 1 ? Math.min(ua, ub) : Math.max(ua, ub);
}

/**
 * Every shown item's label to a layout's sink (_spec/README.md "Render", Layer types: `labels`), an item once: its
 * extent on the section's system (its image's first and last edge: past each other over a cyclic origin), its element
 * and its hit, ordered by its band's place (the layer's order).
 */
function sinkLabels(lay: RangesLay, src: LayerSource, sink: LabelSink, order: number) {
  // an item's rows (its pieces) by key, on the section's system
  const byKey = new Map<string | number, { text: string; row: number; pieces: number[] }>();
  for (let i = 0; i < lay.items; i++) {
    const text = lay.labels[i];
    if (!text || lay.lanes[i] < 0) continue;
    const key = lay.keys[i];
    let g = byKey.get(key);
    if (!g) byKey.set(key, (g = { text, row: i, pieces: [] }));
    if (src.path) g.pieces.push(...toSection(src, lay.starts[i], lay.ends[i], []));
    else g.pieces.push(lay.starts[i], lay.ends[i]);
  }
  const base = { layer: lay.key, section: lay.section, system: lay.system };
  const length = lay.axis.cyclic ? lay.axis.length : undefined;
  for (const [key, g] of byKey) {
    const [s, e] = itemArc(g.pieces, length);
    const hit = itemHit(base, lay, key, lay.starts[g.row], lay.ends[g.row]);
    sink.label(g.text, s, e, { owner: lay.owner, key }, hit, order);
  }
}

/**
 * An item's extent from its pieces: their hull, or on a cyclic system of `length` the circle less its widest gap
 * (an item over the origin, two rows meeting there, is one arc: its start past its end).
 */
function itemArc(pieces: number[], length: number | undefined): [number, number] {
  const n = pieces.length >> 1;
  const order = Array.from({ length: n }, (_, k) => k).sort((a, b) => pieces[2 * a] - pieces[2 * b]);
  const first = pieces[2 * order[0]];
  let last = 0;
  for (let k = 0; k < n; k++) last = Math.max(last, pieces[2 * k + 1]);
  if (length === undefined) return [first, last];
  let best = length - last + first;
  let arc: [number, number] = [first, last];
  let reach = pieces[2 * order[0] + 1];
  for (let k = 1; k < n; k++) {
    const s = pieces[2 * order[k]];
    if (s - reach > best) {
      best = s - reach;
      arc = [s, reach];
    }
    reach = Math.max(reach, pieces[2 * order[k] + 1]);
  }
  return arc;
}

export const RangesType: LayerType<RangesPrep, RangesLay> = {
  name: 'ranges',

  prepare(layer: ResolvedLayer, model: StaticModel): RangesPrep {
    const src = layerSource(model, layer);
    const line = layer.params.line;
    const arrow = layer.params.arrow;
    return {
      ...(typeof line === 'string' ? { line } : {}),
      arrow: arrow === 'end' ? 1 : arrow === 'start' ? 2 : arrow === 'both' ? 3 : 0,
      cyclic: model.systems.get(layer.data?.system ?? '')?.cyclic ?? false,
      ...(itemExtents(model, layer) ? { extents: itemExtents(model, layer) } : {}),
      ...(src ? { src } : {}),
      owner: layer.data?.address ?? '',
      channels: channelForms(layer.theme),
    };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: RangesLay): RangesLay {
    const lay = prev ?? newLay(ctx, prep);
    const { axis, box, config, tokens } = ctx;
    beginItems(lay, ctx, prep);
    lay.font = config.itemFont;
    lay.fontSize = fontSize(config.itemFont);
    lay.items = 0;
    lay.m = 0;
    lay.hn = 0;
    lay.overflow = 0;
    lay.laneCount = 0;
    const src = prep.src;
    if (!src) return lay;

    lay.style = lineStyle(prep.line, ctx);
    lay.arrow = prep.arrow;
    // the items in view and in the stretches its breaks hide (an item there is its break's tick)
    const densities = readItems(lay, ctx, data, src, ['category']);
    const step = lay.step;
    const rows = lay.rows;
    lay.items = rows.n;
    themeItems(lay, prep, ctx);
    for (let j = 0; j < rows.n; j++) lay.ids[j] = laneId(lay.idCache, rows.keys[j], rows.occ[j]);
    const n = lay.items;
    if (lay.idCache.size > 4 * n + 4096) lay.idCache.clear();
    if (lay.arrow) {
      if (lay.first.length < n) {
        lay.first = new Uint8Array(Math.max(n, 2 * lay.first.length));
        lay.last = new Uint8Array(lay.first.length);
      }
      markEnds(rows, prep.extents, prep.cyclic, src.length ?? Number.POSITIVE_INFINITY, lay.first, lay.last);
    }
    const maxLanes = Math.max(1, Math.floor(box.height / config.laneHeight));
    lay.laneHeight = box.height / maxLanes;
    let row = false;
    if (step !== 'density' && !densities.length) {
      lay.order = grown(lay.order, n);
      const order = lay.order.subarray(0, n);
      for (let i = 0; i < n; i++) order[i] = i;
      const starts = lay.starts;
      sortArrayRange(order, 0, n, (o, a, b) => starts[o[a]] - starts[o[b]] || o[a] - o[b]);
      lay.lanes = grown(lay.lanes, n);
      const lanes = lay.lanes.subarray(0, n);
      const assign = (cap: number) => {
        lay.laneCount = lay.sticky.assign(lay.ids, lay.starts, lay.ends, order, cap, ctx.settled, lanes);
        lay.overflow = 0;
        for (let i = 0; i < n; i++) if (lanes[i] < 0) lay.overflow++;
      };
      assign(maxLanes);
      // past the lanes: the last lane becomes a density row of what does not fit; the rest stay items
      if (lay.overflow && config.ranges.overflow === 'extra-lanes' && maxLanes > 1) {
        assign(maxLanes - 1);
        row = true;
      }
    }
    lay.mode =
      step === 'density' || densities.length || (lay.overflow && !row) ? 'density' : step === 'bars' ? 'bars' : 'items';
    // below a position a pixel a run has no ends, no dash and no wave: a plain line or block
    if (lay.mode !== 'items' && lay.style >= 2 && lay.style <= 5) lay.style = Solid;
    // an arrowhead is an item's, by its own width (below): a run of bars keeps the arrows of the items wide enough
    if (lay.mode === 'density') lay.arrow = 0;
    if (lay.style === 4) wave = grown(wave, 2 * Math.ceil(box.width / WaveStep + 8));

    if (lay.mode === 'density') {
      lay.densityY = box.y;
      lay.densityH = box.height;
      densityStrip(lay, src, axis, densities, tokens, false, box.height < config.laneHeight ? ThinStripShare : 1);
      return lay;
    }
    if (row) {
      lay.densityY = box.y + (maxLanes - 1) * lay.laneHeight;
      lay.densityH = lay.laneHeight;
      densityStrip(lay, src, axis, [], tokens, true, 1);
    }

    // marks: each item's pieces on the section, square where cut, round where the item really ends (on its own
    // system an item is its one piece: no mapping and no array an item)
    const px = 1 / axis.dpr;
    for (let i = 0; i < n; i++) {
      const lane = lay.lanes[i];
      if (lane < 0) continue;
      let pieces = itemPiece;
      if (src.path) pieces = toSection(src, lay.starts[i], lay.ends[i], mappedPieces);
      else {
        itemPiece[0] = lay.starts[i];
        itemPiece[1] = lay.ends[i];
      }
      let lo = Number.POSITIVE_INFINITY;
      let hi = Number.NEGATIVE_INFINITY;
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        lo = Math.min(lo, pieces[k]);
        hi = Math.max(hi, pieces[k + 1]);
      }
      const marks = lay.m;
      // the arrowheads go on the piece holding the end (the start) the item really has, in reading order
      const atEnd = (lay.arrow & 1) !== 0 && lay.last[i] === 1;
      const atStart = (lay.arrow & 2) !== 0 && lay.first[i] === 1;
      if (atEnd || atStart) images.set(src, lay.starts[i], lay.ends[i]);
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        axis.piecesIn(pieces[k], pieces[k + 1]);
        const ka = axis.ka;
        const kb = axis.kb;
        for (let a = ka; a <= kb; a++) {
          if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
          // an end the item really has is round; on a reversed axis its start is on the right
          const atLo = pieces[k] === lo;
          const atHi = pieces[k + 1] === hi;
          let round = 0;
          if (!(axis.scut & Cut.left) && (axis.dir === -1 ? atHi : atLo)) round |= Round.left;
          if (!(axis.scut & Cut.right) && (axis.dir === -1 ? atLo : atHi)) round |= Round.right;
          const x0 = axis.snap(axis.sx0);
          const j = nextMark(lay);
          lay.mx0[j] = x0;
          lay.mx1[j] = Math.max(axis.snap(axis.sx1), x0 + px);
          lay.mlane[j] = lane;
          lay.mitem[j] = i;
          lay.mround[j] = round;
          if (lay.style > 0) lay.manchor[j] = anchorOf(axis, a, pieces[k], pieces[k + 1]);
          if ((atEnd || atStart) && lay.mx1[j] - x0 >= ArrowMinPx) {
            let tips = 0;
            const pa = Math.max(pieces[k], axis.pieces[2 * a]);
            const pb = Math.min(pieces[k + 1], axis.pieces[2 * a + 1]);
            const t = atEnd ? images.tip(axis, box, pa, pb, true) : 0;
            const u = atStart ? images.tip(axis, box, pa, pb, false) : 0;
            if (t) tips |= t === 1 ? 2 : 1;
            if (u) tips |= u === 1 ? 2 : 1;
            lay.mtip[j] = tips;
          }
        }
      }
      // an item wholly in a hidden stretch: its break's tick in its lane
      if (lay.m === marks && pieces.length >= 2) tickOn(lay, axis, pieces[0], lane);
    }
    if (lay.mode === 'bars') mergeRuns(lay, px);
    // a layout's label sink (rings: placed outside with leaders) takes every shown item's label, an item once, with its
    // hit; else labels inside marks wide enough, faded past their width (drawn in the labels pass)
    const sink = ctx.labels !== undefined && ctx.labels !== 'inside' ? ctx.labels : undefined;
    if (sink) sinkLabels(lay, src, sink, box.y);
    else if (lay.mode === 'items' && lay.laneHeight >= 10)
      for (let j = 0; j < lay.m; j++) {
        const label = lay.mitem[j] >= 0 ? lay.labels[lay.mitem[j]] : undefined;
        if (label && lay.mx1[j] - lay.mx0[j] >= 24) {
          lay.mlabel[j] = ctx.text.measure(lay.font, label);
          lay.mink[j] = lay.style ? tokens.colors['msui-text'] : inkOn(lay.fills[lay.mitem[j]], tokens);
        }
      }
    // a mark cut at both ends is square
    lay.gap = lay.laneHeight >= 6 ? 2 : 0;
    const h = lay.laneHeight - lay.gap;
    for (let j = 0; j < lay.m; j++)
      if (lay.mround[j])
        lay.mradius[j] = squircleRadius(lay.mx1[j] - lay.mx0[j], h, lay.radius.ratio, lay.radius.max, lay.radius.min);
    // the hit marks: an item's own marks, also inside a run merged with others (_spec/README.md "Render", Lanes and Hits)
    if (lay.mode === 'bars') {
      const s = lay.merge;
      for (let k = 0; k < lay.hn; k++)
        s.radius[k] = s.round[k]
          ? squircleRadius(s.x1[k] - s.x0[k], h, lay.radius.ratio, lay.radius.max, lay.radius.min)
          : 0;
      lay.hx0 = s.x0;
      lay.hx1 = s.x1;
      lay.hlane = s.lane;
      lay.hitem = s.item;
      lay.hround = s.round;
      lay.hradius = s.radius;
      lay.hslack = px;
    } else {
      lay.hn = lay.m;
      lay.hx0 = lay.mx0;
      lay.hx1 = lay.mx1;
      lay.hlane = lay.mlane;
      lay.hitem = lay.mitem;
      lay.hround = lay.mround;
      lay.hradius = lay.mradius;
      lay.hslack = 0;
    }
    lay.index.build(lay.hn, lay.hlane, lay.hx0, lay.hx1, lay.laneCount);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(RangesDraws, lay, out, emphasis);
  },

  items: describeItems,

  hit(lay, at): Hit | undefined {
    const x = at[0];
    const y = at[1];
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height || x < box.x || x >= box.x + box.width) return undefined;
    const position = lay.axis.positionAt(x);
    const base = { layer: lay.key, section: lay.section, system: lay.system, ...(position >= 0 ? { position } : {}) };
    const lane = Math.floor((y - box.y) / lay.laneHeight);
    // a hidden items' tick is hit as its break (within the break's slop), the count of what it hides with it
    const tick = tickHit(lay, x, lane);
    if (tick) return tick;
    if (lay.mode === 'density' || !lay.hn) return base;
    const j = lay.index.near(lane, x, lay.hslack);
    if (j < 0) return base;
    const item = lay.hitem[j];
    if (item < 0) return base;
    return itemHit(base, lay, lay.keys[item], lay.starts[item], lay.ends[item]);
  },

  filtered: (lay) => lay.rows.filtered,
  used: itemsUsed,

  legend(lay, options) {
    // the density strip is a count, not a theme: its key is the strip itself
    if (lay.mode === 'density') return [];
    const painted = new Set<string>();
    if (lay.color?.legend?.kind === 'palette')
      for (let j = 0; j < lay.items; j++) {
        const v = lay.painted[j];
        if (v !== undefined) painted.add(v);
      }
    return channelLegends([lay.color, lay.size], { ...options, layer: lay.key, painted });
  },
};

// Each pass draws in a function of its own, called through a table (`drawPass`: nothing allocated once warm).

/** The data pass: placeholders, the density strip, the marks. */
function drawMarks(lay: RangesLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  const gap = lay.gap;
  // in a function of its own: this one stays small enough to have its calls inlined (a call that is not boxes the
  // numbers it is passed)
  if (lay.tn) drawTicks(lay, out);
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  const a = out.args;
  if (lay.density) {
    a[0] = box.x;
    a[1] = lay.densityY;
    a[2] = box.width;
    a[3] = lay.densityH;
    out.imageArgs(lay.density);
  }
  if (lay.mode === 'density') return;
  if (lay.style !== 0 || lay.arrow) {
    drawStyled(lay, out);
    return;
  }
  // items and merged runs alike are discrete marks: squircles where big enough, round only where a run really ends,
  // so the base mark has the shape its emphasis outlines (_spec/README.md "Render", Shapes)
  const full = lay.laneHeight - gap;
  for (let j = 0; j < lay.m; j++) {
    const item = lay.mitem[j];
    // the size channel shrinks the mark about its lane's middle; its lane and hit box stay (EC-137)
    const h = item >= 0 ? full * lay.sizes[item] : full;
    a[0] = lay.mx0[j];
    a[1] = box.y + lay.mlane[j] * lay.laneHeight + gap / 2 + (full - h) / 2;
    a[2] = lay.mx1[j] - lay.mx0[j];
    a[3] = h;
    a[4] = item >= 0 ? lay.fills[item] : tokens.colors['seq-range'];
    a[5] = NoColor;
    a[6] = 0;
    a[7] = lay.mradius[j];
    a[8] = lay.mround[j];
    out.rectArgs();
  }
}

/**
 * The marks of a layer with a line style or arrowheads: a block with a pointed end where it holds an end the item really
 * has, or the line (its end shortened to where its triangle begins). The dash and the wave run from the content's start
 * (`manchor`), so a pan moves them with the data. No closure, nothing made per mark: the numbers go through `args`.
 */
function drawStyled(lay: RangesLay, out: PassWriter) {
  const { box, tokens } = lay;
  const gap = lay.gap;
  const a = out.args;
  const style = lay.style;
  const forward = lay.axis.dir === 1;
  for (let j = 0; j < lay.m; j++) {
    const item = lay.mitem[j];
    const fill = item >= 0 ? lay.fills[item] : tokens.colors['seq-range'];
    const full = lay.laneHeight - gap;
    const h = item >= 0 ? full * lay.sizes[item] : full;
    const y = box.y + lay.mlane[j] * lay.laneHeight + gap / 2 + (full - h) / 2;
    const cy = y + h / 2;
    const x0 = lay.mx0[j];
    const x1 = lay.mx1[j];
    const tips = lay.mtip[j];
    const head = Math.min(h * 0.9, x1 - x0);
    if (style === 0 || style === 7) {
      // a block (a highlight is a translucent one); a tip makes it an arrow, both tips two overlapping ones
      const color = style === 0 ? fill : ((fill & 0xffffff00) | HighlightAlpha) >>> 0;
      a[0] = x0;
      a[1] = y;
      a[2] = x1 - x0;
      a[3] = h;
      a[4] = color;
      a[5] = NoColor;
      a[6] = 0;
      if (!tips) {
        a[7] = lay.mradius[j];
        a[8] = lay.mround[j];
        out.rectArgs();
        continue;
      }
      if (tips & 2) {
        a[7] = head;
        a[8] = 1;
        arrowArgs(out);
        a[0] = x0;
        a[2] = x1 - x0;
        a[4] = color;
        a[5] = NoColor;
        a[6] = 0;
      }
      if (tips & 1) {
        a[7] = head;
        a[8] = -1;
        arrowArgs(out);
      }
      continue;
    }
    // a line: from where the content starts, over what is in view, ending where a triangle begins
    const l0 = tips & 1 ? x0 + head : x0;
    const l1 = tips & 2 ? x1 - head : x1;
    // the start's triangle, on the side the content starts from: the dash and the wave begin at the line, not under it
    const startTip = forward ? tips & 1 : tips & 2;
    if (style === Solid) out.fill(l0, cy - LineWidth / 2, l1 - l0, LineWidth, fill);
    else if (style === 2 || style === 3) {
      const dash = style === 2 ? Dash : Dot;
      const period = 2 * dash;
      const anchor = lay.manchor[j];
      const m = Math.max(0, (startTip ? Math.ceil : Math.floor)((forward ? l0 - anchor : anchor - l1) / period));
      segment[0] = forward ? anchor + m * period : anchor - m * period;
      segment[1] = cy;
      segment[2] = forward ? l1 : l0;
      segment[3] = cy;
      a[4] = NoColor;
      a[5] = fill;
      a[6] = LineWidth;
      a[7] = dash;
      out.pathArgs(segment, 2, false);
    } else if (style === 4) {
      // the wave from the content's start, a corner every half period, phase kept when it starts before the view
      const period = 2 * WaveStep;
      const anchor = lay.manchor[j];
      const m = Math.max(0, (startTip ? Math.ceil : Math.floor)((forward ? l0 - anchor : anchor - l1) / period));
      const xs = forward ? anchor + m * period : anchor - m * period;
      const xe = forward ? l1 : l0;
      const cap = wave.length / 2;
      const count = Math.min(cap, Math.floor(Math.abs(xe - xs) / WaveStep) + 2);
      for (let i = 0; i < count; i++) {
        const x = forward ? xs + i * WaveStep : xs - i * WaveStep;
        wave[2 * i] = forward ? Math.min(x, xe) : Math.max(x, xe);
        wave[2 * i + 1] = cy + (i & 1 ? WaveAmplitude : -WaveAmplitude);
      }
      a[4] = NoColor;
      a[5] = fill;
      a[6] = 1;
      a[7] = 0;
      out.pathArgs(wave, count, false);
    } else if (style === 5) {
      out.fill(l0, cy - LineWidth, l1 - l0, 1, fill);
      out.fill(l0, cy + LineWidth - 1, l1 - l0, 1, fill);
    } else {
      // a box: the lane's outline, inside the mark
      a[0] = x0 + 0.5;
      a[1] = y + 0.5;
      a[2] = x1 - x0 - 1;
      a[3] = h - 1;
      a[4] = NoColor;
      a[5] = fill;
      a[6] = 1;
      a[7] = 0;
      a[8] = 0;
      out.rectArgs();
    }
    if (tips) {
      // the triangle: a head with no body, over the line's end
      a[1] = y;
      a[3] = h;
      a[4] = fill;
      a[5] = NoColor;
      a[6] = 0;
      if (tips & 2) {
        a[0] = x1 - head;
        a[2] = head;
        a[7] = head;
        a[8] = 1;
        arrowArgs(out);
      }
      if (tips & 1) {
        a[0] = x0;
        a[2] = head;
        a[4] = fill;
        a[5] = NoColor;
        a[6] = 0;
        a[7] = head;
        a[8] = -1;
        arrowArgs(out);
      }
    }
  }
}

/** The labels pass: the labels inside the marks wide enough, faded past their width. */
function drawLabels(lay: RangesLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box } = lay;
  const gap = lay.gap;
  const g = out.args;
  for (let j = 0; j < lay.m; j++) {
    if (!(lay.mlabel[j] > 0)) continue;
    const y = box.y + lay.mlane[j] * lay.laneHeight + gap / 2;
    const h = lay.laneHeight - gap;
    const baseline = y + h / 2 + lay.fontSize * 0.35;
    if (lay.style >= 1 && lay.style <= 5) {
      // a line runs through its label's place: the page's color under the text, so the label reads (the line stops at it)
      g[0] = lay.mx0[j] + 1;
      g[1] = baseline - lay.fontSize * 0.85;
      g[2] = lay.mlabel[j] + 4;
      g[3] = lay.fontSize;
      g[4] = lay.tokens.colors['msui-bg'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
    }
    g[0] = lay.mx0[j] + 3;
    g[1] = baseline;
    g[2] = lay.mlabel[j];
    g[3] = lay.mx1[j] - lay.mx0[j] - 6;
    g[4] = lay.mink[j];
    out.textArgs(lay.labels[lay.mitem[j]] as string, AnchorStart, lay.font);
  }
}

/**
 * The hover pass: an outline around each of the hovered item's own marks (its hit marks: inside a merged run too),
 * concentric.
 */
function drawHover(lay: RangesLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  const gap = lay.gap;
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || !h.element) return;
  // outside the mark by the offset, concentric: it never touches an outline drawn on the mark
  const half = lay.hoverWidth / 2 + lay.hoverOffset;
  const a = out.args;
  for (let j = 0; j < lay.hn; j++) {
    const item = lay.hitem[j];
    if (item < 0 || lay.keys[item] !== h.element.key) continue;
    const r = lay.hradius[j];
    // concentric: the outline's radius is the mark's plus its offset
    a[0] = lay.hx0[j] - half;
    a[1] = box.y + lay.hlane[j] * lay.laneHeight + gap / 2 - half;
    a[2] = lay.hx1[j] - lay.hx0[j] + 2 * half;
    a[3] = lay.laneHeight - gap + 2 * half;
    a[4] = NoColor;
    a[5] = tokens.colors['seq-hover'];
    a[6] = lay.hoverWidth;
    a[7] = r > 0 ? r + half : 0;
    a[8] = lay.hround[j];
    out.rectArgs();
  }
}

/**
 * The selection pass: a selected item of zero length (a boundary, an insertion point) is a hairline at its place over
 * its lane in the accent (EC-127), which the selection's Location cannot hold; the rest is the view's.
 */
function drawSelected(lay: RangesLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const els = emphasis?.elements;
  if (!els || !els.keys.length) return;
  const { box, tokens } = lay;
  const dpr = lay.axis.dpr;
  const a = out.args;
  for (let j = 0; j < lay.m; j++) {
    const item = lay.mitem[j];
    if (item < 0 || lay.starts[item] !== lay.ends[item]) continue;
    let selected = false;
    for (let k = 0; k < els.keys.length && !selected; k++)
      selected = els.keys[k].owner === lay.owner && els.keys[k].key === lay.keys[item];
    if (!selected) continue;
    const x = Math.round(((lay.mx0[j] + lay.mx1[j]) / 2) * dpr) / dpr;
    a[0] = x - 1;
    a[1] = box.y + lay.mlane[j] * lay.laneHeight;
    a[2] = 2;
    a[3] = lay.laneHeight;
    a[4] = tokens.colors['msui-accent'];
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    a[8] = 0;
    out.rectArgs();
  }
}

const RangesDraws = passDraws<RangesLay, Emphasis | undefined>({
  data: drawMarks,
  labels: drawLabels,
  selection: drawSelected,
  hover: drawHover,
});

/**
 * Adds the strip's `v` over the device columns data positions `[s, e)` land on in view: the positions are integers, the
 * count in the strip, `column` inlined (a number made and passed to a call that is not inlined is boxed: an allocation
 * a bin).
 */
export function addToStrip(strip: StripBuffers, src: LayerSource, axis: AxisTransform, s: number, e: number) {
  let pieces = itemPiece;
  if (src.path) pieces = toSection(src, s, e, mappedPieces);
  else {
    itemPiece[0] = s;
    itemPiece[1] = e;
  }
  const diff = strip.diff;
  const v = strip.v;
  const dpr = axis.dpr;
  const last = axis.columns - 1;
  for (let k = 0; k + 1 < pieces.length; k += 2) {
    axis.piecesIn(pieces[k], pieces[k + 1]);
    const ka = axis.ka;
    const kb = axis.kb;
    for (let a = ka; a <= kb; a++) {
      if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
      const f0 = Math.floor((axis.sx0 - axis.x0) * dpr);
      const f1 = Math.floor((axis.sx1 - 1e-6 - axis.x0) * dpr);
      const c0 = f0 < 0 ? 0 : f0 >= last ? last : f0;
      const c1 = Math.max(c0, f1 < 0 ? 0 : f1 >= last ? last : f1);
      diff[c0] += v;
      diff[c1 + 1] -= v;
    }
  }
}

/** An item's piece on its own system, and its pieces through a map: kept (a layout runs to its end). */
const itemPiece: number[] = [0, 0];
const mappedPieces: number[] = [];

/**
 * Merges a lane's marks that touch within a device pixel into one bar a run: at most one mark a device column. The marks
 * are sorted by lane and left edge (ties in the order they were made) into the layout's kept scratch, then written back.
 */
function mergeRuns(lay: RangesLay, px: number) {
  const m = lay.m;
  const s = lay.merge;
  if (s.order.length < m) {
    const n = Math.max(m, 2 * s.order.length);
    s.order = new Int32Array(n);
    s.x0 = new Float32Array(n);
    s.x1 = new Float32Array(n);
    s.lane = new Int32Array(n);
    s.item = new Int32Array(n);
    s.round = new Uint8Array(n);
    s.radius = new Float64Array(n);
    s.tip = new Uint8Array(n);
  }
  const order = s.order;
  for (let j = 0; j < m; j++) order[j] = j;
  const { mlane, mx0 } = lay;
  sortArrayRange(order, 0, m, (o, i, j) => mlane[o[i]] - mlane[o[j]] || mx0[o[i]] - mx0[o[j]] || o[i] - o[j]);
  const { x0, x1, lane, item, round, tip } = s;
  for (let k = 0; k < m; k++) {
    const j = order[k];
    x0[k] = lay.mx0[j];
    x1[k] = lay.mx1[j];
    lane[k] = lay.mlane[j];
    item[k] = lay.mitem[j];
    round[k] = lay.mround[j];
    tip[k] = lay.mtip[j];
  }
  // the item's own marks stay in the scratch: the hit marks (`hn`)
  lay.hn = m;
  lay.m = 0;
  for (let k = 0; k < m; ) {
    let end = x1[k];
    let it = item[k];
    let r = round[k];
    let e = k + 1;
    // an arrow is its item's alone: never merged into a run
    while (!tip[k] && e < m && !tip[e] && lane[e] === lane[k] && x0[e] <= end + px) {
      end = Math.max(end, x1[e]);
      if (item[e] !== it) it = -1;
      r = (r & Round.left) | (round[e] & Round.right);
      e++;
    }
    const j = nextMark(lay);
    lay.mx0[j] = x0[k];
    lay.mx1[j] = end;
    lay.mlane[j] = lane[k];
    lay.mitem[j] = it;
    lay.mround[j] = r;
    lay.mtip[j] = tip[k];
    k = e;
  }
}

/**
 * One heat strip: items (all, or with `overflowOnly` those past the lanes) and a source's density bins, counted a
 * device column, colored low to high. Its buffers are kept while the box's columns do not change.
 */
function densityStrip(
  lay: RangesLay,
  src: LayerSource,
  axis: AxisTransform,
  densities: readonly { table: RangesDensity; level: number }[],
  tokens: SeqTokens,
  overflowOnly: boolean,
  share: number,
) {
  const strip = lay.strip;
  resetStrip(strip, axis.columns);
  strip.v = 1;
  for (let i = 0; i < lay.items; i++)
    if (!overflowOnly || lay.lanes[i] < 0)
      addToStrip(strip, src, axis, lay.starts[i], Math.max(lay.ends[i], lay.starts[i] + 1));
  for (let d = 0; d < densities.length; d++) {
    const { table, level } = densities[d];
    for (let i = 0; i < table.bin.length; i++) {
      strip.v = table.count[i];
      addToStrip(strip, src, axis, table.bin[i] * level, (table.bin[i] + 1) * level);
    }
  }
  lay.density = paintStrip(strip, axis.columns, tokens, share);
}

/** A heat strip's buffers: counts by difference over the device columns, and the image they are painted into. */
export interface StripBuffers {
  image?: ImageData;
  diff: Float64Array;
  counts: Float64Array;
  v: number;
}

/** The strip emptied for `cols` device columns (its buffers kept while their number does not change). */
export function resetStrip(strip: StripBuffers, cols: number) {
  if (strip.diff.length !== cols + 1) {
    strip.diff = new Float64Array(cols + 1);
    strip.counts = new Float64Array(cols);
  } else strip.diff.fill(0);
  if (!strip.image || strip.image.width !== cols)
    strip.image = { width: cols, height: 1, data: new Uint8Array(4 * cols) };
  else strip.image.data.fill(0);
}

/**
 * The strip's counts (`addToStrip` made them) painted low to high in the density tokens: its image. `share` is the part
 * of the ramp the busiest column reaches (1: the whole ramp; a thin band uses the light end, `ThinStripShare`).
 */
export function paintStrip(strip: StripBuffers, cols: number, tokens: SeqTokens, share = 1): ImageData {
  const { diff, counts } = strip;
  const data = strip.image!.data;
  let max = 0;
  let run = 0;
  for (let c = 0; c < cols; c++) {
    run += diff[c];
    counts[c] = run;
    if (run > max) max = run;
  }
  const low = tokens.colors['seq-density-low'];
  const high = tokens.colors['seq-density-high'];
  for (let c = 0; c < cols; c++) {
    if (counts[c] <= 0) continue;
    const color = mixColors(low, high, max > 0 ? (share * counts[c]) / max : 0);
    data[4 * c] = color >>> 24;
    data[4 * c + 1] = (color >>> 16) & 0xff;
    data[4 * c + 2] = (color >>> 8) & 0xff;
    data[4 * c + 3] = color & 0xff;
  }
  return strip.image!;
}

/**
 * The items' channels (_spec/README.md "Themes"): each item's fill, size and painted value, read through its table's
 * row. The channels are compiled once a layout (kept while the tokens and held domains agree), by the first table's
 * declared column types; a field theme over a dictionary-coded column is a color per code, once a table.
 */
function themeItems(lay: RangesLay, prep: RangesPrep, ctx: LayoutContext) {
  const rows = lay.rows;
  const datum = lay.datum;
  const missing = ctx.tokens.colors['seq-missing'];
  datum.bind(rows.n ? rows.fields[0] : undefined);
  const extra = datum.env;
  const color = prep.channels.get('color', 'token:range', ctx, extra);
  const size = prep.channels.get('size', 1, ctx, extra);
  const opacity = prep.channels.get('opacity', 1, ctx, extra);
  lay.color = color;
  lay.size = size;
  lay.extents.reset([color, size, opacity]);
  const painted = color.legend?.kind === 'palette' ? color.legend.field : undefined;
  const cconst = color.constant;
  const sconst = size.constant;
  const oconst = opacity.constant;
  let bound: ItemRows['fields'][number] | null = null;
  let codes: ArrayLike<number> | undefined;
  let dict: readonly string[] | undefined;
  for (let j = 0; j < rows.n; j++) {
    const fields = rows.fields[j];
    if (fields !== bound) {
      bound = fields;
      datum.bind(fields);
      const coded = cconst === undefined && color.field !== undefined ? datum.column(color.field) : undefined;
      codes = coded instanceof CategoryColumn ? coded.codes : undefined;
      dict = coded instanceof CategoryColumn ? coded.dict : undefined;
      if (dict) lay.byCode = codeFills(color, dict, missing, lay.byCode);
    }
    const i = rows.row[j];
    datum.row = i;
    datum.start = rows.starts[j];
    datum.end = rows.ends[j];
    let fill: number;
    if (codes) {
      const code = codes[i];
      fill = code < 0 ? missing : lay.byCode[code];
    } else {
      const c = cconst ?? color.at(datum);
      fill = typeof c === 'number' ? c : missing;
    }
    const o = oconst ?? opacity.at(datum);
    if (typeof o === 'number' && o < 1) fill = withAlpha(fill, ((fill & 0xff) / 255) * Math.max(0, o));
    const z = sconst ?? size.at(datum);
    lay.sizes[j] = typeof z === 'number' ? Math.max(0, Math.min(1, z)) : 1;
    if (painted !== undefined) {
      if (codes) lay.painted[j] = codes[i] < 0 ? undefined : dict![codes[i]];
      else {
        const v = datum.field(painted);
        lay.painted[j] = typeof v === 'string' ? v : undefined;
      }
    }
    if (lay.extents.any) lay.extents.add(datum);
    lay.fills[j] = fill;
  }
  lay.extents.flush(ctx.theme);
}

/** The colors of a dictionary's values through a channel that reads one field: index a code. */
function codeFills(channel: Channel, dict: readonly string[], missing: number, out: Uint32Array): Uint32Array {
  if (out.length < dict.length) out = new Uint32Array(Math.max(dict.length, 2 * out.length));
  let v: string | undefined;
  const d = { value: Number.NaN, letter: '', start: 0, end: 0, field: () => v };
  for (let k = 0; k < dict.length; k++) {
    v = dict[k];
    const c = channel.at(d);
    out[k] = typeof c === 'number' ? c : missing;
  }
  return out;
}
