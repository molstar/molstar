/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArrayRange } from '@molstar/core/data/util/sort';
import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { RangesDensity, RangesList } from '@molstar/seqx-core/data';
import { rangesLocation } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { type AxisTransform, Cut } from './axis.js';
import { type DensityStep, densityStep, fontSize } from './config.js';
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
import { LaneIndex, StickyLanes } from './lanes.js';
import {
  dataRegions,
  type LayerSource,
  layerLevel,
  layerSource,
  placeholderSpans,
  readTiles,
  toSection,
} from './layer-data.js';
import {
  type Box,
  type Emphasis,
  type Hit,
  type LayerType,
  type LayoutContext,
  type LegendSpec,
  type Overlay,
  type StoreReader,
  sectionId,
} from './layer-type.js';
import { reportColor } from './letters.js';
import { drawFailed, drawSpans } from './marks.js';
import { squircleRadius } from './squircle.js';
import { categoryToken, inkOn, mixColors, resolveColor, type SeqTokens } from './tokens.js';

// The `ranges` layer type (_spec/README.md "Render", Layer types): items in sticky lanes. Each item is carried onto
// the section along the layer's path (a domain on its gene: a piece an exon) and drawn as a squircle where it is big
// enough, square below, with square ends where the view or a piece's edge cuts it (the frame rule). Below a device
// pixel a position, the marks of a lane merge into one bar a run; at the density step the layer is one heat strip, a
// device column a pixel. Items past the lanes collapse into a density row at the bottom (`ranges.overflow`:
// `extra-lanes`), or the whole layer does (`layer`). An item's key is its `id`, else its `row`, else its tile and
// index; a key and start seen in an earlier tile is not drawn again (a range repeated in the tiles it touches), and a
// key repeated with another start is another piece of one item (a wrapped row, the exons of a transcript).

export interface RangesPrep {
  readonly src?: LayerSource;
  /** The annotation's address: the owner of its items' element keys. */
  readonly owner: string;
  readonly color?: unknown;
}

type Mode = 'items' | 'bars' | 'density';

export interface RangesLay {
  key: string;
  system: string;
  section: string;
  owner: string;
  /** The data's system: an element's extent is on it. */
  dataSystem: string;
  box: Box;
  axis: AxisTransform;
  step: DensityStep;
  mode: Mode;
  tokens: SeqTokens;
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
  lanes: Int32Array;
  laneCount: number;
  laneHeight: number;
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
  /** The gap between lanes, px. */
  gap: number;
  index: LaneIndex;
  /** `mergeRuns`'s scratch: the marks' order and their fields in it, kept while large enough. */
  merge: {
    order: Int32Array;
    x0: Float32Array;
    x1: Float32Array;
    lane: Int32Array;
    item: Int32Array;
    round: Uint8Array;
  };
  /** The heat strip and its band: the whole box, or the bottom row past the lanes. */
  density?: ImageData;
  densityY: number;
  densityH: number;
  /** The strip's buffers, kept while the box's columns do not change. */
  strip: { image?: ImageData; diff: Float64Array; counts: Float64Array; v: number };
  pending: number[];
  failed: number[];
  sticky: StickyLanes;
  font: string;
  fontSize: number;
  hoverWidth: number;
  hoverOffset: number;
  radius: { ratio: number; max: number; min: number };
  categories: Map<string, number>;
  /** Lane ids by key and occurrence, kept across frames (a layout makes no string an item); cleared when large. */
  idCache: Map<string | number, string[]>;
  /** The items by start, kept. */
  order: Int32Array;
}

/** A table's field by name (looked up once a table, not once an item). */
function fieldOf(fields: RangesList['fields'], name: string): Column | undefined {
  if (!fields) return undefined;
  for (let k = 0; k < fields.length; k++) if (fields[k].name === name) return fields[k].column;
  return undefined;
}

/** A field's value at a row, as text. */
function valueAt(c: Column | undefined, i: number): string | undefined {
  if (!c) return undefined;
  const v = c instanceof CategoryColumn ? c.get(i) : c[i];
  return v === null || v === undefined ? undefined : String(v);
}

/** An item's lane id: its key and which occurrence of it in a tile, made once and kept. */
function laneId(cache: Map<string | number, string[]>, key: string | number, k: number): string {
  let ids = cache.get(key);
  if (!ids) cache.set(key, (ids = []));
  return (ids[k] ??= `${typeof key === 'number' ? 'n' : 's'}${key}#${k}`);
}

function newLay(ctx: LayoutContext, prep: RangesPrep): RangesLay {
  return {
    key: ctx.key,
    system: ctx.axis.system,
    section: sectionId(ctx.section),
    owner: prep.owner,
    dataSystem: prep.src?.dataSystem ?? '',
    box: ctx.box,
    axis: ctx.axis,
    step: 'cells',
    mode: 'items',
    tokens: ctx.tokens,
    items: 0,
    keys: [],
    ids: [],
    starts: [],
    ends: [],
    labels: [],
    fills: [],
    lanes: new Int32Array(64),
    laneCount: 0,
    laneHeight: 0,
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
    gap: 0,
    index: new LaneIndex(),
    merge: {
      order: new Int32Array(0),
      x0: new Float32Array(0),
      x1: new Float32Array(0),
      lane: new Int32Array(0),
      item: new Int32Array(0),
      round: new Uint8Array(0),
    },
    densityY: 0,
    densityH: 0,
    strip: { diff: new Float64Array(0), counts: new Float64Array(0), v: 0.5 },
    pending: [],
    failed: [],
    sticky: new StickyLanes(),
    font: ctx.config.itemFont,
    fontSize: fontSize(ctx.config.itemFont),
    hoverWidth: ctx.config.outline.hover,
    hoverOffset: ctx.config.emphasis.hoverOffset,
    radius: { ratio: ctx.config.squircleRatio, max: ctx.config.squircleMax, min: ctx.config.squircleMin },
    categories: new Map(),
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
  }
  lay.mlabel[j] = 0;
  lay.mradius[j] = 0;
  lay.mink[j] = 0;
  return j;
}

export const RangesType: LayerType<RangesPrep, RangesLay> = {
  name: 'ranges',

  prepare(layer: ResolvedLayer, model: StaticModel): RangesPrep {
    const src = layerSource(model, layer);
    return {
      ...(src ? { src } : {}),
      owner: layer.data?.address ?? '',
      ...(layer.theme.color !== undefined ? { color: layer.theme.color } : {}),
    };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: RangesLay): RangesLay {
    const lay = prev ?? newLay(ctx, prep);
    const { axis, box, config, tokens } = ctx;
    lay.key = ctx.key;
    lay.system = axis.system;
    lay.section = sectionId(ctx.section);
    lay.owner = prep.owner;
    lay.dataSystem = prep.src?.dataSystem ?? '';
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.font = config.itemFont;
    lay.fontSize = fontSize(config.itemFont);
    lay.hoverWidth = config.outline.hover;
    lay.hoverOffset = config.emphasis.hoverOffset;
    lay.items = 0;
    lay.pending.length = 0;
    lay.failed.length = 0;
    lay.m = 0;
    lay.overflow = 0;
    lay.laneCount = 0;
    lay.density = undefined;
    lay.categories.clear();
    const src = prep.src;
    if (!src) return lay;

    let fixed: number | undefined;
    if (prep.color !== undefined) {
      fixed = resolveColor(prep.color, tokens);
      if (fixed === undefined) reportColor(ctx);
    }
    const step = densityStep(axis.scale, axis.dpr, config);
    lay.step = step;
    const regions = dataRegions(src, axis);
    const level = layerLevel(src, data, axis);
    const read = readTiles(ctx, data, src, regions, level);
    placeholderSpans(src, axis, read.pending, lay.pending);
    placeholderSpans(src, axis, read.failed, lay.failed);

    // the items: a key and start seen in an earlier tile is that range again (over the tile's edge); a key repeated
    // with another start, in this tile or a later one, is another piece of it (an exon of one transcript)
    const earlier = new Map<string | number, number[]>();
    const occurrences = new Map<string | number, number>();
    const densities: { table: RangesDensity; level: number }[] = [];
    for (const e of read.entries) {
      const t = e.tables.ranges;
      if (!t) continue;
      if (t.form === 'density') {
        densities.push({ table: t, level: e.region.level });
        continue;
      }
      // this tile's items are `from` to `lay.items`: their keys join the earlier ones after it
      const from = lay.items;
      const labels = fieldOf(t.fields, 'label');
      const categories = fieldOf(t.fields, 'category');
      for (let i = 0; i < t.start.length; i++) {
        const id = t.id?.get(i);
        const key: string | number =
          id !== null && id !== undefined ? id : t.row ? t.row[i] : `${e.region.level}/${e.region.tile}:${i}`;
        if (earlier.get(key)?.includes(t.start[i])) continue;
        const k = occurrences.get(key) ?? 0;
        occurrences.set(key, k + 1);
        const label = valueAt(labels, i) ?? id ?? undefined;
        const category = valueAt(categories, i);
        const j = lay.items++;
        lay.keys[j] = key;
        lay.ids[j] = laneId(lay.idCache, key, k);
        lay.starts[j] = t.start[i];
        lay.ends[j] = t.end[i];
        lay.labels[j] = label;
        const fill = fixed ?? tokens.colors[category === undefined ? 'seq-range' : categoryToken(category)];
        if (category !== undefined) lay.categories.set(category, fill);
        lay.fills[j] = fill;
      }
      for (let j = from; j < lay.items; j++) {
        const starts = earlier.get(lay.keys[j]);
        if (starts) starts.push(lay.starts[j]);
        else earlier.set(lay.keys[j], [lay.starts[j]]);
      }
    }
    const n = lay.items;
    if (lay.idCache.size > 4 * n + 4096) lay.idCache.clear();
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

    if (lay.mode === 'density') {
      lay.densityY = box.y;
      lay.densityH = box.height;
      densityStrip(lay, src, axis, densities, tokens, false);
      return lay;
    }
    if (row) {
      lay.densityY = box.y + (maxLanes - 1) * lay.laneHeight;
      lay.densityH = lay.laneHeight;
      densityStrip(lay, src, axis, [], tokens, true);
    }

    // marks: each item's pieces on the section, square where cut, round where the item really ends (on its own
    // system an item is its one piece: no mapping and no array an item)
    const px = 1 / axis.dpr;
    const nk = axis.offsets.length;
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
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        for (let a = 0; a < nk; a++) {
          if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
          let round = 0;
          if (!(axis.scut & Cut.left) && pieces[k] === lo) round |= Round.left;
          if (!(axis.scut & Cut.right) && pieces[k + 1] === hi) round |= Round.right;
          const x0 = axis.snap(axis.sx0);
          const j = nextMark(lay);
          lay.mx0[j] = x0;
          lay.mx1[j] = Math.max(axis.snap(axis.sx1), x0 + px);
          lay.mlane[j] = lane;
          lay.mitem[j] = i;
          lay.mround[j] = round;
        }
      }
    }
    if (lay.mode === 'bars') mergeRuns(lay, px);
    // labels inside marks wide enough, faded past their width (drawn in the labels pass)
    if (lay.mode === 'items' && lay.laneHeight >= 10)
      for (let j = 0; j < lay.m; j++) {
        const label = lay.mitem[j] >= 0 ? lay.labels[lay.mitem[j]] : undefined;
        if (label && lay.mx1[j] - lay.mx0[j] >= 24) {
          lay.mlabel[j] = ctx.text.measure(lay.font, label);
          lay.mink[j] = inkOn(lay.fills[lay.mitem[j]], tokens);
        }
      }
    // a mark cut at both ends is square
    lay.gap = lay.laneHeight >= 6 ? 2 : 0;
    const h = lay.laneHeight - lay.gap;
    for (let j = 0; j < lay.m; j++)
      if (lay.mround[j])
        lay.mradius[j] = squircleRadius(lay.mx1[j] - lay.mx0[j], h, lay.radius.ratio, lay.radius.max, lay.radius.min);
    lay.index.build(lay.m, lay.mlane, lay.mx0, lay.mx1, lay.laneCount);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(RangesDraws, lay, out, emphasis);
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height || x < box.x || x >= box.x + box.width) return undefined;
    const position = lay.axis.positionAt(x);
    const base = { layer: lay.key, section: lay.section, system: lay.system, ...(position >= 0 ? { position } : {}) };
    if (lay.mode === 'density') return base;
    const lane = Math.floor((y - box.y) / lay.laneHeight);
    const j = lay.index.at(lane, x);
    if (j < 0) return base;
    const item = lay.mitem[j];
    if (item < 0) return base;
    return {
      ...base,
      element: { owner: lay.owner, key: lay.keys[item] },
      extent: rangesLocation(lay.dataSystem, [lay.starts[item], lay.ends[item]]),
    };
  },

  legend(lay): LegendSpec | undefined {
    if (lay.mode === 'density')
      return { kind: 'ramp', colors: [lay.tokens.colors['seq-density-low'], lay.tokens.colors['seq-density-high']] };
    if (!lay.categories.size) return undefined;
    const labels = [...lay.categories.keys()].sort();
    return { kind: 'categories', labels, colors: labels.map((l) => lay.categories.get(l) as number) };
  },
};

// Each pass draws in a function of its own, called through a table (`drawPass`: nothing allocated once warm).

/** The data pass: placeholders, the density strip, the marks. */
function drawMarks(lay: RangesLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  const gap = lay.gap;
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
  // items and merged runs alike are discrete marks: squircles where big enough, round only where a run really ends,
  // so the base mark has the shape its emphasis outlines (_spec/README.md "Render", Shapes)
  for (let j = 0; j < lay.m; j++) {
    const item = lay.mitem[j];
    a[0] = lay.mx0[j];
    a[1] = box.y + lay.mlane[j] * lay.laneHeight + gap / 2;
    a[2] = lay.mx1[j] - lay.mx0[j];
    a[3] = lay.laneHeight - gap;
    a[4] = item >= 0 ? lay.fills[item] : tokens.colors['seq-range'];
    a[5] = NoColor;
    a[6] = 0;
    a[7] = lay.mradius[j];
    a[8] = lay.mround[j];
    out.rectArgs();
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
    g[0] = lay.mx0[j] + 3;
    g[1] = y + h / 2 + lay.fontSize * 0.35;
    g[2] = lay.mlabel[j];
    g[3] = lay.mx1[j] - lay.mx0[j] - 6;
    out.textArgs(lay.labels[lay.mitem[j]] as string, Anchor.start, lay.font, lay.mink[j]);
  }
}

/** The hover pass: an outline around each mark of the hovered item, concentric. */
function drawHover(lay: RangesLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  const gap = lay.gap;
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || !h.element) return;
  // outside the mark by the offset, concentric: it never touches an outline drawn on the mark
  const half = lay.hoverWidth / 2 + lay.hoverOffset;
  const a = out.args;
  for (let j = 0; j < lay.m; j++) {
    const item = lay.mitem[j];
    if (item < 0 || lay.keys[item] !== h.element.key) continue;
    const r = lay.mradius[j];
    // concentric: the outline's radius is the mark's plus its offset
    a[0] = lay.mx0[j] - half;
    a[1] = box.y + lay.mlane[j] * lay.laneHeight + gap / 2 - half;
    a[2] = lay.mx1[j] - lay.mx0[j] + 2 * half;
    a[3] = lay.laneHeight - gap + 2 * half;
    a[4] = NoColor;
    a[5] = tokens.colors['seq-hover'];
    a[6] = lay.hoverWidth;
    a[7] = r > 0 ? r + half : 0;
    a[8] = lay.mround[j];
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
function addToStrip(strip: RangesLay['strip'], src: LayerSource, axis: AxisTransform, s: number, e: number) {
  let pieces = itemPiece;
  if (src.path) pieces = toSection(src, s, e, mappedPieces);
  else {
    itemPiece[0] = s;
    itemPiece[1] = e;
  }
  const diff = strip.diff;
  const v = strip.v;
  const nk = axis.offsets.length;
  const dpr = axis.dpr;
  const last = axis.columns - 1;
  for (let k = 0; k + 1 < pieces.length; k += 2)
    for (let a = 0; a < nk; a++) {
      if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
      const f0 = Math.floor((axis.sx0 - axis.x0) * dpr);
      const f1 = Math.floor((axis.sx1 - 1e-6 - axis.x0) * dpr);
      const c0 = f0 < 0 ? 0 : f0 >= last ? last : f0;
      const c1 = Math.max(c0, f1 < 0 ? 0 : f1 >= last ? last : f1);
      diff[c0] += v;
      diff[c1 + 1] -= v;
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
  }
  const order = s.order;
  for (let j = 0; j < m; j++) order[j] = j;
  const { mlane, mx0 } = lay;
  sortArrayRange(order, 0, m, (o, i, j) => mlane[o[i]] - mlane[o[j]] || mx0[o[i]] - mx0[o[j]] || o[i] - o[j]);
  const { x0, x1, lane, item, round } = s;
  for (let k = 0; k < m; k++) {
    const j = order[k];
    x0[k] = lay.mx0[j];
    x1[k] = lay.mx1[j];
    lane[k] = lay.mlane[j];
    item[k] = lay.mitem[j];
    round[k] = lay.mround[j];
  }
  lay.m = 0;
  for (let k = 0; k < m; ) {
    let end = x1[k];
    let it = item[k];
    let r = round[k];
    let e = k + 1;
    while (e < m && lane[e] === lane[k] && x0[e] <= end + px) {
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
) {
  const cols = axis.columns;
  const strip = lay.strip;
  if (strip.diff.length !== cols + 1) {
    strip.diff = new Float64Array(cols + 1);
    strip.counts = new Float64Array(cols);
  } else strip.diff.fill(0);
  if (!strip.image || strip.image.width !== cols)
    strip.image = { width: cols, height: 1, data: new Uint8Array(4 * cols) };
  else strip.image.data.fill(0);
  const { diff, counts } = strip;
  const data = strip.image.data;
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
    const color = mixColors(low, high, max > 0 ? counts[c] / max : 0);
    data[4 * c] = color >>> 24;
    data[4 * c + 1] = (color >>> 16) & 0xff;
    data[4 * c + 2] = (color >>> 8) & 0xff;
    data[4 * c + 3] = color & 0xff;
  }
  lay.density = strip.image;
}
