/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ValuesTable } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { type DensityStep, densityStep, perPosition, type RenderConfig } from './config.js';
import { ValueColorFallback, ValueEnv } from './layer-color.js';
import { drawPass, grown, NoColor, type PassName, type PassWriter, passDraws, Round } from './draw-lists.js';
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
  type Overlay,
  type StoreReader,
  sectionId,
} from './layer-type.js';
import { drawFailed, drawSpans, hoverCell } from './marks.js';
import type { SeqTokens } from './tokens.js';
import { sortArrayRange } from '@molstar/core/data/util/sort';
import type { Channel } from '@molstar/seqx-core/theme/index';
import { channelForms, channelLegends, Extents, type LayerChannels, PlainDatum } from './channels.js';

// The `values` layer type (_spec/README.md "Render", Layer types): bars of continuous data. A position of 4 CSS px or
// more is a bar (a 1 px gap, rounded tops); below, one stepped profile over the device columns' extents (a faint area
// under a thin top line), from per-position values or a coarser level's bins (`values.aggregate`: `minmax`, so negative
// extremes show, or `max`). The top tenth of the well is headroom, a hairline marks the baseline, and a missing value
// (NaN) is a gap, never a stub or a zero.
// The domain is the layer's `domain` param, else open: computed from the values in view when the view
// settles, held while it does not, so heights do not shimmer during a pan. Its `scale` names a scale node: its y domain
// is the scale's (an open end the union over the layers naming it, held from the last settle). Its color is its
// `color` channel: uniform (one fill), or by `value` (a bar's color its value's, `col`; a line's color changes with it).
//
// `encoding` `line` and `area` (EC-57) draw the same bars as a line through the cells' centers: held flat over the half
// cells at its ends and at gaps (a hole in the data, a break of the axis), broken by a missing value (which keeps its
// mark at the baseline), its color changing halfway between two positions; `area` fills under it at `FillAlpha`. At
// the bars and density steps the line follows each device column's maximum, the column's extent still drawn.

/** How a values layer draws: bars (default), a line, or a line over the area under it. */
export type ValuesEncoding = 'bars' | 'line' | 'area';

/** The line's thickness, CSS px, and the area's opacity (the channels' defaults). */
const LineWidth = 1.5;

export interface ValuesPrep {
  readonly encoding: ValuesEncoding;
  readonly src?: LayerSource;
  readonly channels: LayerChannels;
  readonly domain?: readonly [number, number];
  /** The scale node its y reads, by address. */
  readonly scale?: string;
}

export interface ValuesLay {
  key: string;
  system: string;
  section: string;
  box: Box;
  axis: AxisTransform;
  step: DensityStep;
  tokens: SeqTokens;
  fill: number;
  encoding: ValuesEncoding;
  /** A line's points in x order (`ord`: indices of the bars) and the color of each bar (the fill, or its value's). */
  ord: Int32Array;
  col: Uint32Array;
  /** The color channel; `col` holds a bar's color where it is not uniform. */
  color?: Channel;
  perBar: boolean;
  datum: PlainDatum;
  extents: Extents;
  /** Bars followed on the screen by a break of the axis (a line does not cross it: it is held, as at a gap). */
  cut: Uint8Array;
  domain: [number, number] | undefined;
  /** Bars: snapped x extents and the low and high value of each (equal for a position; NaN: missing). */
  n: number;
  x0: Float32Array;
  x1: Float32Array;
  lo: Float32Array;
  v: Float32Array;
  pending: number[];
  failed: number[];
  /** A device column's low and high, kept while the box's columns do not change. */
  colLo: Float32Array;
  colHi: Float32Array;
  hoverWidth: number;
  hoverOffset: number;
  /** What the hover's squircle reads (its rule). */
  config: RenderConfig;
}

/** The next bar's index, the bars' arrays grown to hold it (the caller writes its numbers: nothing passed is boxed). */
function nextBar(lay: ValuesLay): number {
  const j = lay.n++;
  if (j >= lay.x0.length) {
    lay.x0 = grown(lay.x0, j + 1);
    lay.x1 = grown(lay.x1, j + 1);
    lay.lo = grown(lay.lo, j + 1);
    lay.v = grown(lay.v, j + 1);
  }
  return j;
}

export const ValuesType: LayerType<ValuesPrep, ValuesLay> = {
  name: 'values',

  prepare(layer: ResolvedLayer, model: StaticModel): ValuesPrep {
    const src = layerSource(model, layer);
    const d = layer.params.domain;
    const domain =
      Array.isArray(d) && d.length === 2 && typeof d[0] === 'number' && typeof d[1] === 'number' && d[0] < d[1]
        ? ([d[0], d[1]] as const)
        : undefined;
    const e = layer.params.encoding;
    const scale = typeof layer.params.scale === 'string' ? layer.params.scale : undefined;
    return {
      encoding: e === 'line' || e === 'area' ? e : 'bars',
      ...(src ? { src } : {}),
      channels: channelForms(layer.theme),
      ...(domain ? { domain } : {}),
      ...(scale ? { scale } : {}),
    };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: ValuesLay): ValuesLay {
    const { axis, box, config, tokens } = ctx;
    const lay: ValuesLay = prev ?? {
      key: ctx.key,
      system: axis.system,
      section: sectionId(ctx.section),
      box,
      axis,
      step: 'cells',
      tokens,
      fill: 0,
      encoding: prep.encoding,
      ord: new Int32Array(0),
      col: new Uint32Array(0),
      perBar: false,
      datum: new PlainDatum(),
      extents: new Extents(),
      cut: new Uint8Array(0),
      domain: undefined,
      n: 0,
      x0: new Float32Array(64),
      x1: new Float32Array(64),
      lo: new Float32Array(64),
      v: new Float32Array(64),
      pending: [],
      failed: [],
      colLo: new Float32Array(0),
      colHi: new Float32Array(0),
      hoverWidth: config.outline.hover,
      hoverOffset: config.emphasis.hoverOffset,
      config,
    };
    lay.key = ctx.key;
    lay.system = axis.system;
    lay.section = sectionId(ctx.section);
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.hoverWidth = config.outline.hover;
    lay.hoverOffset = config.emphasis.hoverOffset;
    lay.config = config;
    lay.n = 0;
    lay.pending.length = 0;
    lay.failed.length = 0;
    const color = prep.channels.get('color', ValueColorFallback, ctx, ValueEnv);
    lay.color = color;
    lay.perBar = color.constant === undefined;
    lay.fill = typeof color.constant === 'number' ? color.constant : tokens.colors['seq-value'];
    lay.encoding = prep.encoding;
    const src = prep.src;
    if (!src) return lay;
    const step = densityStep(axis.scale, axis.dpr, config);
    lay.step = step;
    const regions = dataRegions(src, axis);
    const level = perPosition(step)
      ? layerLevel(src, data, axis, true)
      : layerLevel(src, data, axis, false, config.values.coarser);
    const read = readTiles(ctx, data, src, regions, level);
    placeholderSpans(src, axis, read.pending, lay.pending);
    placeholderSpans(src, axis, read.failed, lay.failed);
    const columns = !perPosition(step);
    if (lay.colHi.length !== axis.columns) {
      lay.colHi = new Float32Array(axis.columns);
      lay.colLo = new Float32Array(axis.columns);
    }
    const colHi = lay.colHi.fill(Number.NaN);
    const colLo = lay.colLo.fill(Number.NaN);
    const px = 1 / axis.dpr;
    const minmax = config.values.aggregate === 'minmax';
    const pl = placer;
    pl.lay = lay;
    pl.axis = axis;
    pl.src = src;
    pl.direct = config.values.direct && !src.path;
    pl.cols = columns;
    pl.minmax = minmax;
    pl.px = px;
    for (const e of read.entries) {
      const t: ValuesTable | undefined = e.tables.values;
      if (!t) continue;
      if (t.form === 'positions') {
        // the table's positions in each region, in order (the regions are sorted and disjoint)
        const end = t.start + t.value.length;
        for (let r = 0; r + 1 < regions.length; r += 2) {
          const a = Math.max(t.start, regions[r]);
          const b = Math.min(end, regions[r + 1]);
          for (let p = a; p < b; p++) {
            const v = t.value[p - t.start];
            pl.lo = v;
            pl.hi = v;
            place(pl, p, p + 1);
          }
        }
      } else {
        const w = e.region.level;
        for (let i = 0; i < t.max.length; i++) {
          const a = t.start + i * w;
          if (!overlapsRegions(regions, a, a + w)) continue;
          const none = !(t.count[i] > 0);
          pl.hi = none ? Number.NaN : t.max[i];
          pl.lo = minmax ? (none ? Number.NaN : t.min[i]) : pl.hi;
          place(pl, a, a + w);
        }
      }
    }
    pl.lay = undefined;
    pl.axis = undefined;
    pl.src = undefined;
    if (columns)
      for (let c = 0; c < colHi.length; c++) {
        const j = nextBar(lay);
        lay.x0[j] = box.x + c * px;
        lay.x1[j] = box.x + (c + 1) * px;
        lay.lo[j] = minmax ? colLo[c] : colHi[c];
        lay.v[j] = colHi[c];
      }
    if (lay.perBar) colorBars(lay, color, ctx);
    lay.domain = yDomain(prep, ctx, lay);
    orderLine(lay);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(ValuesDraws, lay, out, emphasis);
  },

  hit(lay, at): Hit | undefined {
    const x = at[0];
    const y = at[1];
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    return { layer: lay.key, section: lay.section, position, system: lay.system };
  },

  legend(lay, options) {
    return channelLegends([lay.color], { ...options, layer: lay.key });
  },

  domain(lay) {
    return lay.domain;
  },
};

// Each pass draws in a function of its own, called through a table (`drawPass`: nothing allocated once warm).

/** The data pass: bars, or the line (and area) through the same marks. */
function drawData(lay: ValuesLay, out: PassWriter) {
  if (lay.encoding !== 'bars') drawLine(lay, out);
  else drawBars(lay, out);
}

/**
 * The look of the bars (_spec/README.md "Render", Layer types, `values`): a stepped profile below `BarMinPx` CSS px a
 * value (an area at `FillAlpha` under a top line of `ProfileLine`), bars from it (a gap, rounded tops, at `BarAlpha`).
 * An alpha multiplies the series color's own (`(c & 0xffffff00) | (c & 255) * alpha / 255`, written out where a color
 * is stored: a color past a small integer is boxed in a local or an argument, an allocation a frame).
 */
const FillAlpha = 89;
const BarAlpha = 153;
/** The top line's width and the bars' gap and corner radius, CSS px. */
const ProfileLine = 0.75;
const BarGap = 1;
const BarRound = 2;
/** A profile's run goes on while its colors differ by at most this (the sum over red, green and blue). */
const RunColorStep = 48;

/** A value takes bars from this many CSS px; below, the profile. */
const BarMinPx = 4;
/** The share of the height the domain's top takes: the rest is headroom (10%). */
const Fill = 0.9;

/** The hairline at the baseline (zero, else the domain's low), a device pixel thick. No number is passed to a call. */
function drawBaseline(lay: ValuesLay, out: PassWriter) {
  const { box } = lay;
  const lo = lay.domain ? lay.domain[0] : 0;
  const hi = lay.domain ? lay.domain[1] : 1;
  const base = lo <= 0 && hi >= 0 ? 0 : lo;
  const span = (hi - lo || 1) / Fill;
  const yBase = box.y + box.height - ((Math.max(lo, Math.min(hi, base)) - lo) / span) * box.height;
  const px = 1 / lay.axis.dpr;
  const a = out.args;
  a[0] = box.x;
  a[1] = yBase - px;
  a[2] = box.width;
  a[3] = px;
  a[4] = lay.tokens.colors['seq-sequence'];
  a[5] = NoColor;
  a[6] = 0;
  a[7] = 0;
  a[8] = 0;
  out.rectArgs();
}

/**
 * The data pass of bars: placeholders, a hairline at the baseline, then a profile (runs of consecutive values: the area
 * under the steps, the top line over it; a run goes on across small steps of a per-value color) or bars. A missing
 * value is a gap. Nothing allocated once warm: the numbers go through `out.args`, the points through `linePoints`.
 */
function drawBars(lay: ValuesLay, out: PassWriter) {
  const { box, tokens } = lay;
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  const lo = lay.domain ? lay.domain[0] : 0;
  const hi = lay.domain ? lay.domain[1] : 1;
  const base = lo <= 0 && hi >= 0 ? 0 : lo;
  const span = (hi - lo || 1) / Fill;
  const y1 = box.y + box.height;
  const bh = box.height;
  const yBase = y1 - ((Math.max(lo, Math.min(hi, base)) - lo) / span) * bh;
  const px = 1 / lay.axis.dpr;
  drawBaseline(lay, out);
  const a = out.args;
  const n = lay.n;
  const bars = perPosition(lay.step) && lay.axis.scale >= BarMinPx;
  if (bars) {
    for (let j = 0; j < n; j++) {
      const v = lay.v[j];
      if (Number.isNaN(v)) continue;
      // from the baseline to the value (a column's extent past it too): the free end is the round one
      const yv = y1 - ((Math.max(lo, Math.min(hi, v)) - lo) / span) * bh;
      const yl = y1 - ((Math.max(lo, Math.min(hi, Math.min(lay.lo[j], base))) - lo) / span) * bh;
      const top = Math.min(yv, yBase);
      const w = lay.x1[j] - lay.x0[j];
      // a bar narrower than the gap and a device pixel (a cut at the window's edge) takes no gap
      const gap = w > BarGap + px ? BarGap : 0;
      a[0] = lay.x0[j] + gap / 2;
      a[1] = top;
      a[2] = w - gap;
      a[3] = Math.max(yv, yl, yBase) - top;
      if (lay.perBar) a[4] = ((lay.col[j] & 0xffffff00) | ((((lay.col[j] & 255) * BarAlpha + 127) / 255) | 0)) >>> 0;
      else a[4] = ((lay.fill & 0xffffff00) | ((((lay.fill & 255) * BarAlpha + 127) / 255) | 0)) >>> 0;
      a[5] = NoColor;
      a[6] = 0;
      a[7] = Math.max(0, Math.min(BarRound, a[2] / 2, a[3] / 2));
      a[8] = Round.both | (v < base ? Round.cutTop : Round.cutBottom);
      out.rectArgs();
    }
    return;
  }
  // a column's extent past the baseline on the other side of its value (`minmax`: a negative extreme under a positive
  // maximum, the maximum of a column below zero up to it): a block of the area's alpha, from the value's end
  if (!perPosition(lay.step)) {
    for (let j = 0; j < n; j++) {
      if (Number.isNaN(lay.v[j]) || !(lay.lo[j] < base)) continue;
      const yv = y1 - ((Math.max(lo, Math.min(hi, lay.v[j])) - lo) / span) * bh;
      const yl = y1 - ((Math.max(lo, Math.min(hi, lay.lo[j])) - lo) / span) * bh;
      // the profile's area covers the baseline to the maximum: this block is what lies beyond it
      const from = Math.max(yv, yBase);
      a[0] = lay.x0[j];
      a[1] = from;
      a[2] = lay.x1[j] - lay.x0[j];
      a[3] = yl - from;
      if (lay.perBar) a[4] = ((lay.col[j] & 0xffffff00) | ((((lay.col[j] & 255) * FillAlpha + 127) / 255) | 0)) >>> 0;
      else a[4] = ((lay.fill & 0xffffff00) | ((((lay.fill & 255) * FillAlpha + 127) / 255) | 0)) >>> 0;
      a[5] = NoColor;
      a[6] = 0;
      a[7] = 0;
      a[8] = 0;
      out.rectArgs();
    }
  }
  const pts = linePoints;
  const ord = lay.ord;
  let count = 0;
  let prev = -1;
  for (let k = 0; k <= n; k++) {
    const j = k < n ? ord[k] : -1;
    const valid = j >= 0 && !Number.isNaN(lay.v[j]);
    const next =
      valid &&
      count > 0 &&
      prev >= 0 &&
      !lay.cut[prev] &&
      lay.x0[j] - lay.x1[prev] <= px * 1.01 &&
      lay.x0[j] >= lay.x1[prev] - px * 1.01;
    // a per-value color: a new run (the area and line in the new color) where it has moved by more than a step
    const recolor =
      next &&
      Math.abs((lay.col[j] >>> 24) - (a[0] >>> 24)) +
        Math.abs(((lay.col[j] >>> 16) & 255) - ((a[0] >>> 16) & 255)) +
        Math.abs(((lay.col[j] >>> 8) & 255) - ((a[0] >>> 8) & 255)) >
        RunColorStep;
    if (count > 0 && (!next || recolor)) {
      // the area under the steps down to the baseline, then the line over it
      pts[2 * count] = pts[2 * count - 2];
      pts[2 * count + 1] = yBase;
      pts[2 * count + 2] = pts[0];
      pts[2 * count + 3] = yBase;
      a[4] = ((a[0] & 0xffffff00) | ((((a[0] & 255) * FillAlpha + 127) / 255) | 0)) >>> 0;
      a[5] = NoColor;
      a[6] = 0;
      a[7] = 0;
      out.pathArgs(pts, count + 2, true);
      a[4] = NoColor;
      a[5] = a[0];
      a[6] = ProfileLine;
      a[7] = 0;
      out.pathArgs(pts, count, false);
      count = 0;
      if (recolor) {
        // the next run's line starts where this one ended: the top line stays one line
        pts[0] = lay.x1[prev];
        pts[1] = y1 - ((Math.max(lo, Math.min(hi, lay.v[prev])) - lo) / span) * bh;
        count = 1;
      }
    }
    if (!valid) {
      prev = -1;
      continue;
    }
    const y = y1 - ((Math.max(lo, Math.min(hi, lay.v[j])) - lo) / span) * bh;
    if (count < 2) a[0] = lay.col[j];
    pts[2 * count] = lay.x0[j];
    pts[2 * count + 1] = y;
    pts[2 * count + 2] = lay.x1[j];
    pts[2 * count + 3] = y;
    count += 2;
    prev = j;
  }
}

/**
 * The hover's outline on the value hovered: its bar's or its step's own extent, from the baseline to the value's end
 * (a position that has no value in view, the cell's column as every layer's).
 */
function hoverValue(lay: ValuesLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || h.position === undefined) return;
  const axis = lay.axis;
  if (h.system !== axis.system || !axis.cell(h.position)) return;
  const cx = (axis.sx0 + axis.sx1) / 2;
  let j = -1;
  for (let i = 0; i < lay.n; i++)
    if (lay.x0[i] <= cx && cx < lay.x1[i] && !Number.isNaN(lay.v[i])) {
      j = i;
      break;
    }
  if (j < 0) {
    hoverCell(lay, out, emphasis);
    return;
  }
  const { box } = lay;
  const lo = lay.domain ? lay.domain[0] : 0;
  const hi = lay.domain ? lay.domain[1] : 1;
  const base = lo <= 0 && hi >= 0 ? 0 : lo;
  const span = (hi - lo || 1) / Fill;
  const y1 = box.y + box.height;
  const bh = box.height;
  const yBase = y1 - ((Math.max(lo, Math.min(hi, base)) - lo) / span) * bh;
  const yv = y1 - ((Math.max(lo, Math.min(hi, lay.v[j])) - lo) / span) * bh;
  const yl = y1 - ((Math.max(lo, Math.min(hi, Math.min(lay.lo[j], base))) - lo) / span) * bh;
  const top = Math.min(yv, yBase);
  const bottom = Math.max(yv, yl, yBase);
  // outside the extent by the offset, as every hover outline; a sliver of a bar still gets a visible frame
  const d = lay.hoverWidth / 2 + lay.hoverOffset;
  const c = lay.config;
  const cw = lay.x1[j] - lay.x0[j];
  const eh = Math.max(bottom - top, 1 / lay.axis.dpr);
  // the extent's squircle by the frame rule (as `hoverCell`), the outline's its plus `d`
  const r =
    cw < c.squircleMin || eh < c.squircleMin ? 0 : Math.min(eh * c.squircleRatio, c.squircleMax, cw / 2, eh / 2);
  const g = out.args;
  g[0] = lay.x0[j] - d;
  g[1] = top - d;
  g[2] = cw + 2 * d;
  g[3] = eh + 2 * d;
  g[4] = NoColor;
  g[5] = lay.tokens.colors['seq-hover'];
  g[6] = lay.hoverWidth;
  g[7] = r > 0 ? r + d : 0;
  g[8] = Round.both;
  out.rectArgs();
}

const ValuesDraws = passDraws<ValuesLay, Emphasis | undefined>({ data: drawData, hover: hoverValue });

/** A layout's placement state: one, reused (a layout runs to its end), so placing a value allocates nothing. */
interface Placer {
  lay: ValuesLay | undefined;
  axis: AxisTransform | undefined;
  src: LayerSource | undefined;
  /** The data is on the section's own system: placed without the mapping. */
  direct: boolean;
  /** A bar a device column (else a bar a position). */
  cols: boolean;
  minmax: boolean;
  px: number;
  /** The value being placed, its low (its high with `max`) and high: set by the layout, read by `mark`. */
  lo: number;
  hi: number;
}

const placer: Placer = {
  lay: undefined,
  axis: undefined,
  src: undefined,
  direct: true,
  cols: false,
  minmax: true,
  px: 0.5,
  lo: 0.5,
  hi: 0.5,
};

/** A placement's pieces on the section, through the map (kept). */
const pieces: number[] = [];

// Placing passes no number made per value to a call (a call that is not inlined boxes each: an allocation a value):
// positions are integers, the value is in the placer, the extent in the axis's `sx0`, `sx1`, and `snap`, `column` and
// the bar's store are inlined on purpose.

/**
 * The value in the placer at data positions `[a, b)`: every piece of the axis it lands on in view. On the section's own
 * system the positions are the section's (no mapping, no arrays): the same marks as through `toSection`.
 */
function place(pl: Placer, a: number, b: number) {
  const axis = pl.axis!;
  if (pl.direct) {
    axis.piecesIn(a, b);
    const ka = axis.ka;
    const kb = axis.kb;
    for (let k = ka; k <= kb; k++) if (axis.spanIn(k, a, b)) mark(pl, axis);
    return;
  }
  toSection(pl.src!, a, b, pieces);
  for (let i = 0; i + 1 < pieces.length; i += 2) {
    axis.piecesIn(pieces[i], pieces[i + 1]);
    const ka = axis.ka;
    const kb = axis.kb;
    for (let k = ka; k <= kb; k++) if (axis.spanIn(k, pieces[i], pieces[i + 1])) mark(pl, axis);
  }
}

/** The extent `spanIn` left in the axis (unsnapped): a bar a position, or its device columns' extents widened. */
function mark(pl: Placer, axis: AxisTransform) {
  const lay = pl.lay!;
  const lo = pl.lo;
  const hi = pl.hi;
  const dpr = axis.dpr;
  if (!pl.cols) {
    const x0 = Math.round(axis.sx0 * dpr) / dpr;
    const j = nextBar(lay);
    lay.x0[j] = x0;
    lay.x1[j] = Math.max(Math.round(axis.sx1 * dpr) / dpr, x0 + pl.px);
    lay.lo[j] = lo;
    lay.v[j] = hi;
    return;
  }
  if (Number.isNaN(hi)) return;
  const colHi = lay.colHi;
  const colLo = lay.colLo;
  const last = axis.columns - 1;
  const f0 = Math.floor((axis.sx0 - axis.x0) * dpr);
  const f1 = Math.floor((axis.sx1 - 1e-6 - axis.x0) * dpr);
  const c0 = f0 < 0 ? 0 : f0 >= last ? last : f0;
  const c1 = Math.max(c0, f1 < 0 ? 0 : f1 >= last ? last : f1);
  for (let c = c0; c <= c1; c++) {
    if (Number.isNaN(colHi[c]) || hi > colHi[c]) colHi[c] = hi;
    if (Number.isNaN(colLo[c]) || lo < colLo[c]) colLo[c] = lo;
  }
}

function overlapsRegions(regions: readonly number[], a: number, b: number): boolean {
  for (let i = 0; i + 1 < regions.length; i += 2) if (a < regions[i + 1] && b > regions[i]) return true;
  return false;
}

/** The values in view, min to max; `[0, 1]` when there are none, widened around a single value. */
function openDomain(lay: ValuesLay): [number, number] {
  let lo = Number.POSITIVE_INFINITY;
  let hi = Number.NEGATIVE_INFINITY;
  for (let j = 0; j < lay.n; j++) {
    const v = lay.v[j];
    if (Number.isNaN(v)) continue;
    if (lay.lo[j] < lo) lo = lay.lo[j];
    if (v > hi) hi = v;
  }
  if (lo > hi) return [0, 1];
  if (lo === hi) return lo >= 0 ? [0, hi || 1] : [lo, 0];
  return [Math.min(lo, 0), hi];
}

/** A line's scratch points (x, y pairs), kept and grown at layout. */
let linePoints = new Float64Array(256);

/** The bars in x order and each one's color: the line's points. */
function orderLine(lay: ValuesLay) {
  const n = lay.n;
  if (lay.ord.length < n) {
    lay.ord = new Int32Array(Math.max(n, 2 * lay.ord.length));
    lay.cut = new Uint8Array(lay.ord.length);
  }
  if (!lay.perBar) {
    if (lay.col.length < n) lay.col = new Uint32Array(lay.ord.length);
    lay.col.fill(lay.fill, 0, n);
  }
  lay.cut.fill(0, 0, n);
  const ord = lay.ord;
  const x0 = lay.x0;
  // the bars come in x order (columns, a section's own positions): the sort, and its closure, only where they do not
  let sorted = true;
  for (let i = 0; i < n; i++) {
    ord[i] = i;
    if (i > 0 && x0[i] < x0[i - 1]) sorted = false;
  }
  if (!sorted) {
    const o = ord.subarray(0, n);
    sortArrayRange(o, 0, n, (q, a, b) => x0[q[a]] - x0[q[b]] || q[a] - q[b]);
  }
  // the breaks of the axis in view, each after the bar whose right edge is at it (a position a bar only)
  if (perPosition(lay.step)) {
    const breaks = lay.axis.breaks([]);
    const px = 1 / lay.axis.dpr;
    for (let i = 0; i + 2 < breaks.length; i += 3) {
      const x = lay.axis.snap(breaks[i]);
      let lo = 0;
      let hi = n;
      while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (lay.x1[ord[mid]] < x - px / 2) lo = mid + 1;
        else hi = mid;
      }
      if (lo < n && Math.abs(lay.x1[ord[lo]] - x) <= px / 2) lay.cut[ord[lo]] = 1;
    }
  }
  linePoints = grown(linePoints, 4 * n + 16);
}

/**
 * The line's data pass: placeholders, the missing marks and the columns' extents, then each run of consecutive values as
 * paths, one a color: the first cell held flat over its first half, the centers, the last held over its second half;
 * a color changes at the boundary between two cells, where the line is halfway between their centers. The area is a
 * polygon under each path, down to the baseline. No closure, no destructuring: nothing allocated once warm.
 */
function drawLine(lay: ValuesLay, out: PassWriter) {
  const { box, tokens } = lay;
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  const lo = lay.domain ? lay.domain[0] : 0;
  const hi = lay.domain ? lay.domain[1] : 1;
  const base = lo <= 0 && hi >= 0 ? 0 : lo;
  const span = (hi - lo || 1) / Fill;
  const y1 = box.y + box.height;
  const bh = box.height;
  const yBase = y1 - ((Math.max(lo, Math.min(hi, base)) - lo) / span) * bh;
  const each = perPosition(lay.step);
  const px = 1 / lay.axis.dpr;
  drawBaseline(lay, out);
  const a = out.args;
  const pts = linePoints;
  const area = lay.encoding === 'area';
  const n = lay.n;
  const ord = lay.ord;
  let count = 0;
  let color = 0;
  for (let k = 0; k < n; k++) {
    const j = ord[k];
    const v = lay.v[j];
    // a missing value is a gap: it breaks the line, nothing is drawn
    if (Number.isNaN(v)) continue;
    // the extent a column holds below its maximum, as a bar over it
    if (!each && lay.lo[j] < v) {
      const top = y1 - ((Math.max(lo, Math.min(hi, Math.max(v, base))) - lo) / span) * bh;
      const h = y1 - ((Math.max(lo, Math.min(hi, Math.min(lay.lo[j], base))) - lo) / span) * bh - top;
      out.fill(lay.x0[j], top, lay.x1[j] - lay.x0[j], h, lay.col[j]);
    }
  }
  // the runs: walk the ordered bars once more, drawing a path where the run or its color ends
  let prev = -1;
  for (let k = 0; k <= n; k++) {
    const j = k < n ? ord[k] : -1;
    const valid = j >= 0 && !Number.isNaN(lay.v[j]);
    const joins =
      valid &&
      prev >= 0 &&
      !lay.cut[prev] &&
      lay.x0[j] - lay.x1[prev] <= px * 1.01 &&
      lay.x0[j] >= lay.x1[prev] - px * 1.01;
    if (count > 0 && (!joins || (valid && lay.col[j] !== color))) {
      const last = prev;
      const yl = y1 - ((Math.max(lo, Math.min(hi, lay.v[last])) - lo) / span) * bh;
      if (joins) {
        // the color changes at the boundary: the line goes on to it, halfway to the next center
        const yn = y1 - ((Math.max(lo, Math.min(hi, lay.v[j])) - lo) / span) * bh;
        pts[2 * count] = lay.x1[last];
        pts[2 * count + 1] = (yl + yn) / 2;
        count++;
      } else {
        pts[2 * count] = lay.x1[last];
        pts[2 * count + 1] = yl;
        count++;
      }
      flush(out, a, pts, count, color, area, yBase);
      count = 0;
    }
    if (!valid) {
      prev = -1;
      continue;
    }
    const v = lay.v[j];
    const y = y1 - ((Math.max(lo, Math.min(hi, v)) - lo) / span) * bh;
    if (count === 0) {
      color = lay.col[j];
      if (joins && prev >= 0) {
        // a run starting at a color change: from the boundary, halfway between the centers
        const yp = y1 - ((Math.max(lo, Math.min(hi, lay.v[prev])) - lo) / span) * bh;
        pts[0] = lay.x0[j];
        pts[1] = (yp + y) / 2;
        count = 1;
      } else {
        // the first cell held flat over its first half
        pts[0] = lay.x0[j];
        pts[1] = y;
        count = 1;
      }
    }
    pts[2 * count] = (lay.x0[j] + lay.x1[j]) / 2;
    pts[2 * count + 1] = y;
    count++;
    prev = j;
  }
}

/** One path of the line (and the polygon under it, behind): its points `pts` (`count` of them) in `color`. */
function flush(
  out: PassWriter,
  a: Float64Array,
  pts: Float64Array,
  count: number,
  color: number,
  area: boolean,
  yBase: number,
) {
  if (area) {
    // the polygon: the line, down to the baseline at its last x, back along it to its first x
    pts[2 * count] = pts[2 * count - 2];
    pts[2 * count + 1] = yBase;
    pts[2 * count + 2] = pts[0];
    pts[2 * count + 3] = yBase;
    a[4] = ((color & 0xffffff00) | ((((color & 255) * FillAlpha + 127) / 255) | 0)) >>> 0;
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    out.pathArgs(pts, count + 2, true);
  }
  a[4] = NoColor;
  a[5] = color;
  a[6] = LineWidth;
  a[7] = 0;
  out.pathArgs(pts, count, false);
}

/** Each bar's color by its value (missing: the missing token) into `col`, the open extents of the color's domain folded. */
function colorBars(lay: ValuesLay, color: Channel, ctx: LayoutContext) {
  const d = lay.datum;
  const missing = lay.tokens.colors['seq-missing'];
  if (lay.col.length < lay.n) lay.col = new Uint32Array(Math.max(lay.n, 2 * lay.col.length));
  lay.extents.reset([color]);
  for (let j = 0; j < lay.n; j++) {
    d.value = lay.v[j];
    const c = color.at(d);
    lay.col[j] = typeof c === 'number' ? c : missing;
    if (lay.extents.any) lay.extents.add(d);
  }
  // held, so a ramp's open domain is the data's extent (a menu's Viridis on 0 to 100 data), not [0, 1]
  if (lay.extents.any) lay.extents.flush(ctx.theme);
}

/**
 * The y domain: the layer's own `domain`, else its scale's (an explicit end wins; an open one the union held over the
 * layers naming the scale, the layer's own values in view until the view first settles), else open over its values.
 */
function yDomain(prep: ValuesPrep, ctx: LayoutContext, lay: ValuesLay): [number, number] {
  if (prep.domain) return [prep.domain[0], prep.domain[1]];
  const own = ctx.settled || !lay.domain ? openDomain(lay) : lay.domain;
  if (!prep.scale) return own;
  if (lay.n) ctx.theme.extent(prep.scale, own[0], own[1]);
  const d = ctx.theme.scale(prep.scale)?.domain;
  const held = ctx.theme.domain(prep.scale) ?? own;
  const lo = typeof d?.[0] === 'number' ? d[0] : held[0];
  const hi = typeof d?.[1] === 'number' ? d[1] : held[1];
  return lo < hi ? [lo, hi] : own;
}
