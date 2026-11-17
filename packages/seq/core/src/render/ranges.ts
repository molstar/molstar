/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { RangesDensity, RangesList } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import { type AxisTransform, Cut } from './axis.js';
import { type DensityStep, densityStep } from './config.js';
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
import type { Box, Emphasis, Hit, LayerType, LayoutContext, LegendSpec, Overlay, StoreReader } from './layer-type.js';
import { reportColor } from './letters.js';
import { drawSpans } from './marks.js';
import { squircleRadius } from './squircle.js';
import { categoryToken, inkOn, mixColors, resolveColor, type SeqTokens, withAlpha } from './tokens.js';

// The `ranges` layer type (_spec/README.md "Render", Layer types): items in sticky lanes. Each item is carried onto
// the section along the layer's path (a domain on its gene: a piece an exon) and drawn as a squircle where it is big
// enough, square below, with square ends where the view or a piece's edge cuts it (the frame rule). Below a device
// pixel a position, the marks of a lane merge into one bar a run; at the density step, or when the items overflow the
// lanes, the layer is one heat strip, a device column a pixel. Ranges cut across tiles are deduped by `id` or `row`.

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
  owner: string;
  box: Box;
  axis: AxisTransform;
  step: DensityStep;
  mode: Mode;
  tokens: SeqTokens;
  /** The items in view, deduped: element key, extent on the data's system, label, fill. */
  keys: (string | number)[];
  ids: string[];
  starts: number[];
  ends: number[];
  labels: (string | undefined)[];
  fills: number[];
  lanes: Int32Array;
  laneCount: number;
  laneHeight: number;
  /** Items that did not fit the lanes (the layer then draws its density). */
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
  density?: ImageData;
  pending: number[];
  failed: number[];
  sticky: StickyLanes;
  font: string;
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
    owner: prep.owner,
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
    pending: [],
    failed: [],
    sticky: new StickyLanes(),
    font: ctx.config.itemFont,
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
    Object.assign(lay, { key: ctx.key, system: axis.system, owner: prep.owner, box, axis, tokens });
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
    const step = densityStep(axis.devicePx, config.density);
    lay.step = step;
    const regions = dataRegions(src, axis);
    const level = layerLevel(src, data, axis);
    const read = readTiles(ctx, data, src, regions, level);
    placeholderSpans(src, axis, read.pending, lay.pending);
    placeholderSpans(src, axis, read.failed, lay.failed);

    // the items, each once: a range repeated in the tiles it touches has one key, start and end
    const seen = new Set<string>();
    const densities: { table: RangesDensity; level: number }[] = [];
    for (const e of read.entries) {
      const t = e.tables.ranges;
      if (!t) continue;
      if (t.form === 'density') {
        densities.push({ table: t, level: e.region.level });
        continue;
      }
      for (let i = 0; i < t.start.length; i++) {
        const key: string | number = t.id ? (t.id.get(i) ?? i) : t.row ? t.row[i] : `${t.start[i]}:${t.end[i]}`;
        const id = `${typeof key === 'number' ? 'n' : 's'}${key}@${t.start[i]}:${t.end[i]}`;
        if (seen.has(id)) continue;
        seen.add(id);
        const label = fieldAt(t.fields, 'label', i) ?? (t.id ? (t.id.get(i) ?? undefined) : undefined);
        const category = fieldAt(t.fields, 'category', i);
        lay.keys.push(key);
        lay.ids.push(id);
        lay.starts.push(t.start[i]);
        lay.ends.push(t.end[i]);
        lay.labels.push(label);
        const fill = fixed ?? tokens.colors[category === undefined ? 'seq-range' : categoryToken(category)];
        if (category !== undefined) lay.categories.set(category, fill);
        lay.fills.push(fill);
      }
    }
    const n = lay.keys.length;
    const maxLanes = Math.max(1, Math.floor(box.height / config.laneHeight));
    lay.laneHeight = box.height / maxLanes;
    if (step !== 'density' && !densities.length) {
      const order = Array.from({ length: n }, (_, i) => i).sort((a, b) => lay.starts[a] - lay.starts[b] || a - b);
      lay.lanes = grow(lay.lanes, n);
      const lanes = lay.lanes.subarray(0, n);
      lay.laneCount = lay.sticky.assign(lay.ids, lay.starts, lay.ends, order, maxLanes, ctx.settled, lanes);
      for (let i = 0; i < n; i++) if (lanes[i] < 0) lay.overflow++;
    }
    lay.mode = step === 'density' || densities.length || lay.overflow ? 'density' : step === 'bars' ? 'bars' : 'items';

    if (lay.mode === 'density') {
      densityStrip(lay, src, axis, densities, tokens);
      return lay;
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
    // labels inside marks wide enough, faded past their width
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
    if (pass === 'data') {
      drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
      drawSpans(out, lay.failed, box.y, box.height, withAlpha(tokens.colors['seq-conflict'], 0.25));
      if (lay.density) {
        out.image(box.x, box.y, box.width, box.height, lay.density);
        return;
      }
      const gap = lay.laneHeight >= 6 ? 2 : 0;
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
        if (lay.mlabel[j] > 0) {
          const label = lay.labels[item] as string;
          const size = Number(/(\d+(?:\.\d+)?)px/.exec(lay.font)?.[1] ?? 10);
          out.text(
            label,
            lay.mx0[j] + 3,
            y + h / 2 + size * 0.35,
            Anchor.start,
            lay.font,
            inkOn(fill, tokens),
            lay.mlabel[j],
            w - 6,
          );
        }
      }
      return;
    }
    if (pass === 'hover') {
      const h = emphasis?.hover;
      if (!h || h.layer !== lay.key || !h.element) return;
      const gap = lay.laneHeight >= 6 ? 2 : 0;
      for (let j = 0; j < lay.m; j++) {
        const item = lay.mitem[j];
        if (item < 0 || lay.keys[item] !== h.element.key) continue;
        const y = box.y + lay.mlane[j] * lay.laneHeight + gap / 2;
        const hh = lay.laneHeight - gap;
        const w = lay.mx1[j] - lay.mx0[j];
        const r = lay.mround[j] ? squircleRadius(w, hh, lay.radius.ratio, lay.radius.max, lay.radius.min) : 0;
        // concentric: the outline's radius is the mark's plus its offset
        out.rect(
          lay.mx0[j] - 1,
          y - 1,
          w + 2,
          hh + 2,
          NoColor,
          tokens.colors['seq-hover'],
          1.5,
          r > 0 ? r + 1 : 0,
          lay.mround[j],
        );
      }
    }
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height || x < box.x || x >= box.x + box.width) return undefined;
    const position = lay.axis.positionAt(x);
    const base = { layer: lay.key, system: lay.system, ...(position >= 0 ? { position } : {}) };
    if (lay.mode === 'density') return base;
    const lane = Math.floor((y - box.y) / lay.laneHeight);
    const j = lay.index.at(lane, x);
    if (j < 0) return base;
    const item = lay.mitem[j];
    return item < 0 ? base : { ...base, element: { owner: lay.owner, key: lay.keys[item] } };
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

/** One heat strip: items (or a source's density bins) counted a device column, colored low to high. */
function densityStrip(
  lay: RangesLay,
  src: LayerSource,
  axis: AxisTransform,
  densities: readonly { table: RangesDensity; level: number }[],
  tokens: SeqTokens,
) {
  const cols = axis.columns;
  const diff = new Float64Array(cols + 1);
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
  for (let i = 0; i < lay.keys.length; i++) add(lay.starts[i], Math.max(lay.ends[i], lay.starts[i] + 1), 1);
  for (const { table, level } of densities)
    for (let i = 0; i < table.bin.length; i++) add(table.bin[i] * level, (table.bin[i] + 1) * level, table.count[i]);
  let max = 0;
  let run = 0;
  const counts = new Float64Array(cols);
  for (let c = 0; c < cols; c++) {
    run += diff[c];
    counts[c] = run;
    if (run > max) max = run;
  }
  const data = new Uint8Array(4 * cols);
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
  lay.density = { width: cols, height: 1, data };
}
