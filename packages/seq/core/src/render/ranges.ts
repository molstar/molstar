/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { RangesDensity, RangesList } from '@molstar/seqx-core/data';
import { rangesLocation } from '@molstar/seqx-core/query';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { type AxisTransform, Cut } from './axis.js';
import { type DensityStep, densityStep, fontSize } from './config.js';
import { Anchor, type ImageData, NoColor, type PassName, type PassWriter, Round } from './draw-lists.js';
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
import { drawSpans } from './marks.js';
import { squircleRadius } from './squircle.js';
import { categoryToken, inkOn, mixColors, resolveColor, type SeqTokens, withAlpha } from './tokens.js';

// The `ranges` layer type (_spec/README.md "Render", Layer types): items in sticky lanes. Each item is carried onto
// the section along the layer's path (a domain on its gene: a piece an exon) and drawn as a squircle where it is big
// enough, square below, with square ends where the view or a piece's edge cuts it (the frame rule). Below a device
// pixel a position, the marks of a lane merge into one bar a run; at the density step the layer is one heat strip, a
// device column a pixel. Items past the lanes collapse into a density row at the bottom (`ranges.overflow`:
// `extra-lanes`), or the whole layer does (`layer`). An item's key is its `id`, else its `row`, else its tile and
// index; a key seen in an earlier tile is not drawn again (a range repeated in the tiles it touches), and within one
// tile a repeated key is another piece of one item (a wrapped row).

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
  /** The items in view: element key, lane id, extent on the data's system, label, fill. */
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
  index: LaneIndex;
  /** The heat strip and its band: the whole box, or the bottom row past the lanes. */
  density?: ImageData;
  densityY: number;
  densityH: number;
  /** The strip's buffers, kept while the box's columns do not change. */
  strip: { image?: ImageData; diff: Float64Array; counts: Float64Array };
  pending: number[];
  failed: number[];
  sticky: StickyLanes;
  font: string;
  fontSize: number;
  hoverWidth: number;
  radius: { ratio: number; max: number; min: number };
  categories: Map<string, number>;
}

function grow<T extends Int32Array | Uint8Array | Float32Array>(a: T, n: number): T {
  return n <= a.length ? a : new (a.constructor as new (n: number) => T)(Math.max(n, 2 * a.length));
}

function fieldAt(fields: RangesList['fields'], name: string, i: number): string | undefined {
  const f = fields?.find((x) => x.name === name);
  if (!f) return undefined;
  const c: Column = f.column;
  const v = c instanceof CategoryColumn ? c.get(i) : c[i];
  return v === null || v === undefined ? undefined : String(v);
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
    index: new LaneIndex(),
    densityY: 0,
    densityH: 0,
    strip: { diff: new Float64Array(0), counts: new Float64Array(0) },
    pending: [],
    failed: [],
    sticky: new StickyLanes(),
    font: ctx.config.itemFont,
    fontSize: fontSize(ctx.config.itemFont),
    hoverWidth: ctx.config.outline.hover,
    radius: { ratio: ctx.config.squircleRatio, max: ctx.config.squircleMax, min: ctx.config.squircleMin },
    categories: new Map(),
  };
}

function pushMark(lay: RangesLay, x0: number, x1: number, lane: number, item: number, round: number) {
  const j = lay.m++;
  if (j >= lay.mx0.length) {
    lay.mx0 = grow(lay.mx0, j + 1);
    lay.mx1 = grow(lay.mx1, j + 1);
    lay.mlane = grow(lay.mlane, j + 1);
    lay.mitem = grow(lay.mitem, j + 1);
    lay.mround = grow(lay.mround, j + 1);
    lay.mlabel = grow(lay.mlabel, j + 1);
  }
  lay.mx0[j] = x0;
  lay.mx1[j] = x1;
  lay.mlane[j] = lane;
  lay.mitem[j] = item;
  lay.mround[j] = round;
  lay.mlabel[j] = 0;
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
    lay.keys.length = 0;
    lay.ids.length = 0;
    lay.starts.length = 0;
    lay.ends.length = 0;
    lay.labels.length = 0;
    lay.fills.length = 0;
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

    // the items: a key seen in an earlier tile is that range again; repeated within a tile, another piece of it
    const earlier = new Set<string | number>();
    const occurrences = new Map<string | number, number>();
    const densities: { table: RangesDensity; level: number }[] = [];
    for (const e of read.entries) {
      const t = e.tables.ranges;
      if (!t) continue;
      if (t.form === 'density') {
        densities.push({ table: t, level: e.region.level });
        continue;
      }
      const here: (string | number)[] = [];
      for (let i = 0; i < t.start.length; i++) {
        const id = t.id?.get(i);
        const key: string | number =
          id !== null && id !== undefined ? id : t.row ? t.row[i] : `${e.region.level}/${e.region.tile}:${i}`;
        if (earlier.has(key)) continue;
        here.push(key);
        const k = occurrences.get(key) ?? 0;
        occurrences.set(key, k + 1);
        const label = fieldAt(t.fields, 'label', i) ?? id ?? undefined;
        const category = fieldAt(t.fields, 'category', i);
        lay.keys.push(key);
        lay.ids.push(`${typeof key === 'number' ? 'n' : 's'}${key}#${k}`);
        lay.starts.push(t.start[i]);
        lay.ends.push(t.end[i]);
        lay.labels.push(label);
        const fill = fixed ?? tokens.colors[category === undefined ? 'seq-range' : categoryToken(category)];
        if (category !== undefined) lay.categories.set(category, fill);
        lay.fills.push(fill);
      }
      for (const k of here) earlier.add(k);
    }
    const n = lay.keys.length;
    const maxLanes = Math.max(1, Math.floor(box.height / config.laneHeight));
    lay.laneHeight = box.height / maxLanes;
    let row = false;
    if (step !== 'density' && !densities.length) {
      const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => lay.starts[a] - lay.starts[b] || a - b);
      lay.lanes = grow(lay.lanes, n);
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

    // marks: each item's pieces on the section, square where cut, round where the item really ends
    const pieces: number[] = [];
    const spans: number[] = [];
    const px = 1 / axis.dpr;
    for (let i = 0; i < n; i++) {
      const lane = lay.lanes[i];
      if (lane < 0) continue;
      toSection(src, lay.starts[i], lay.ends[i], pieces);
      let lo = Number.POSITIVE_INFINITY;
      let hi = Number.NEGATIVE_INFINITY;
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        lo = Math.min(lo, pieces[k]);
        hi = Math.max(hi, pieces[k + 1]);
      }
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        spans.length = 0;
        axis.spans(pieces[k], pieces[k + 1], spans);
        for (let s = 0; s + 2 < spans.length; s += 3) {
          let round = 0;
          if (!(spans[s + 2] & Cut.left) && pieces[k] === lo) round |= Round.left;
          if (!(spans[s + 2] & Cut.right) && pieces[k + 1] === hi) round |= Round.right;
          const x0 = axis.snap(spans[s]);
          const x1 = Math.max(axis.snap(spans[s + 1]), x0 + px);
          pushMark(lay, x0, x1, lane, i, round);
        }
      }
    }
    if (lay.mode === 'bars') mergeRuns(lay, px);
    // labels inside marks wide enough, faded past their width (drawn in the labels pass)
    if (lay.mode === 'items' && lay.laneHeight >= 10)
      for (let j = 0; j < lay.m; j++) {
        const label = lay.mitem[j] >= 0 ? lay.labels[lay.mitem[j]] : undefined;
        if (label && lay.mx1[j] - lay.mx0[j] >= 24) lay.mlabel[j] = ctx.text.measure(lay.font, label);
      }
    lay.index.build(lay.m, lay.mlane, lay.mx0, lay.mx1, lay.laneCount);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    const { box, tokens } = lay;
    const gap = lay.laneHeight >= 6 ? 2 : 0;
    if (pass === 'data') {
      drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
      drawSpans(out, lay.failed, box.y, box.height, withAlpha(tokens.colors['seq-conflict'], 0.25));
      if (lay.density) out.image(box.x, lay.densityY, box.width, lay.densityH, lay.density);
      if (lay.mode === 'density') return;
      for (let j = 0; j < lay.m; j++) {
        const y = box.y + lay.mlane[j] * lay.laneHeight + gap / 2;
        const h = lay.laneHeight - gap;
        const w = lay.mx1[j] - lay.mx0[j];
        const item = lay.mitem[j];
        const fill = item >= 0 ? lay.fills[item] : tokens.colors['seq-range'];
        // a mark cut at both ends is square
        const r =
          lay.mode === 'items' && lay.mround[j]
            ? squircleRadius(w, h, lay.radius.ratio, lay.radius.max, lay.radius.min)
            : 0;
        out.rect(lay.mx0[j], y, w, h, fill, NoColor, 0, r, lay.mround[j]);
      }
      return;
    }
    if (pass === 'labels') {
      for (let j = 0; j < lay.m; j++) {
        if (!(lay.mlabel[j] > 0)) continue;
        const item = lay.mitem[j];
        const y = box.y + lay.mlane[j] * lay.laneHeight + gap / 2;
        const h = lay.laneHeight - gap;
        const w = lay.mx1[j] - lay.mx0[j];
        out.text(
          lay.labels[item] as string,
          lay.mx0[j] + 3,
          y + h / 2 + lay.fontSize * 0.35,
          Anchor.start,
          lay.font,
          inkOn(lay.fills[item], tokens),
          lay.mlabel[j],
          w - 6,
        );
      }
      return;
    }
    if (pass === 'hover') {
      const h = emphasis?.hover;
      if (!h || h.layer !== lay.key || !h.element) return;
      const half = lay.hoverWidth / 2;
      for (let j = 0; j < lay.m; j++) {
        const item = lay.mitem[j];
        if (item < 0 || lay.keys[item] !== h.element.key) continue;
        const y = box.y + lay.mlane[j] * lay.laneHeight + gap / 2;
        const hh = lay.laneHeight - gap;
        const w = lay.mx1[j] - lay.mx0[j];
        const r = lay.mround[j] ? squircleRadius(w, hh, lay.radius.ratio, lay.radius.max, lay.radius.min) : 0;
        // concentric: the outline's radius is the mark's plus its offset
        out.rect(
          lay.mx0[j] - half,
          y - half,
          w + 2 * half,
          hh + 2 * half,
          NoColor,
          tokens.colors['seq-hover'],
          lay.hoverWidth,
          r > 0 ? r + half : 0,
          lay.mround[j],
        );
      }
    }
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

/** Merges a lane's marks that touch within a device pixel into one bar a run: at most one mark a device column. */
function mergeRuns(lay: RangesLay, px: number) {
  const order = Array.from({ length: lay.m }, (_, j) => j).sort(
    (a, b) => lay.mlane[a] - lay.mlane[b] || lay.mx0[a] - lay.mx0[b],
  );
  const x0 = Array.from(order, (j) => lay.mx0[j]);
  const x1 = Array.from(order, (j) => lay.mx1[j]);
  const lane = Array.from(order, (j) => lay.mlane[j]);
  const item = Array.from(order, (j) => lay.mitem[j]);
  const round = Array.from(order, (j) => lay.mround[j]);
  lay.m = 0;
  for (let k = 0; k < order.length; ) {
    let end = x1[k];
    let it = item[k];
    let r = round[k];
    let e = k + 1;
    while (e < order.length && lane[e] === lane[k] && x0[e] <= end + px) {
      end = Math.max(end, x1[e]);
      if (item[e] !== it) it = -1;
      r = (r & Round.left) | (round[e] & Round.right);
      e++;
    }
    pushMark(lay, x0[k], end, lane[k], it, r);
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
  const pieces: number[] = [];
  const spans: number[] = [];
  const add = (s: number, e: number, v: number) => {
    toSection(src, s, e, pieces);
    for (let k = 0; k + 1 < pieces.length; k += 2) {
      spans.length = 0;
      axis.spans(pieces[k], pieces[k + 1], spans);
      for (let i = 0; i + 2 < spans.length; i += 3) {
        const c0 = axis.column(spans[i]);
        const c1 = Math.max(c0, axis.column(spans[i + 1] - 1e-6));
        diff[c0] += v;
        diff[c1 + 1] -= v;
      }
    }
  };
  for (let i = 0; i < lay.keys.length; i++)
    if (!overflowOnly || lay.lanes[i] < 0) add(lay.starts[i], Math.max(lay.ends[i], lay.starts[i] + 1), 1);
  for (const { table, level } of densities)
    for (let i = 0; i < table.bin.length; i++) add(table.bin[i] * level, (table.bin[i] + 1) * level, table.count[i]);
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
