/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ValuesTable } from '@molstar/seqx-core/data';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { type DensityStep, densityStep, perPosition } from './config.js';
import { drawPass, grown, type PassName, type PassWriter, passDraws } from './draw-lists.js';
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
import { drawFailed, drawSpans, hoverCell } from './marks.js';
import { resolveColor, type SeqTokens } from './tokens.js';

// The `values` layer type (_spec/README.md "Render", Layer types): a bar a value, square (continuous data). Where a
// position is a device pixel or more, a bar a position; below, a bar a device column over the column's extent from
// the baseline (`values.aggregate`: `minmax`, so negative extremes show, or `max`), from per-position values or a
// coarser level's bins. Missing (NaN) is a thin neutral mark at the baseline,
// never a zero. The domain is the layer's `domain` param, else open: computed from the values in view when the view
// settles, held while it does not, so heights do not shimmer during a pan.

export interface ValuesPrep {
  readonly src?: LayerSource;
  readonly color?: unknown;
  readonly domain?: readonly [number, number];
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
    return {
      ...(src ? { src } : {}),
      ...(layer.theme.color !== undefined ? { color: layer.theme.color } : {}),
      ...(domain ? { domain } : {}),
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
    };
    lay.key = ctx.key;
    lay.system = axis.system;
    lay.section = sectionId(ctx.section);
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.hoverWidth = config.outline.hover;
    lay.n = 0;
    lay.pending.length = 0;
    lay.failed.length = 0;
    let fill = tokens.colors['seq-value'];
    if (prep.color !== undefined) {
      const c = resolveColor(prep.color, tokens);
      if (c === undefined) reportColor(ctx);
      else fill = c;
    }
    lay.fill = fill;
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
    if (prep.domain) lay.domain = [prep.domain[0], prep.domain[1]];
    else if (ctx.settled || !lay.domain) lay.domain = openDomain(lay);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(ValuesDraws, lay, out, emphasis);
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    return { layer: lay.key, section: lay.section, position, system: lay.system };
  },

  legend(lay): LegendSpec | undefined {
    return lay.domain ? { kind: 'ramp', domain: lay.domain, colors: [lay.fill] } : undefined;
  },
};

// Each pass draws in a function of its own, called through a table (`drawPass`: nothing allocated once warm).

/** The data pass: placeholders, and a bar a value from the baseline (missing: a mark at it). */
function drawBars(lay: ValuesLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  // no closure, no destructuring, no call with a number made here (`yOf` inlined on purpose) and one call site a
  // primitive (a site first run after its function was compiled stays a call): nothing allocated once warm
  const lo = lay.domain ? lay.domain[0] : 0;
  const hi = lay.domain ? lay.domain[1] : 1;
  const base = lo <= 0 && hi >= 0 ? 0 : lo;
  const span = hi - lo || 1;
  const y1 = box.y + box.height;
  const bh = box.height;
  const yBase = y1 - ((Math.max(lo, Math.min(hi, base)) - lo) / span) * bh;
  const missing = tokens.colors['seq-missing'];
  const each = perPosition(lay.step);
  for (let j = 0; j < lay.n; j++) {
    const v = lay.v[j];
    let top = yBase - 2;
    let h = 2;
    let fill = missing;
    if (Number.isNaN(v)) {
      if (!each) continue;
    } else {
      // the extent from the baseline: a column's low below it shows (negative extremes)
      top = y1 - ((Math.max(lo, Math.min(hi, Math.max(v, base))) - lo) / span) * bh;
      h = y1 - ((Math.max(lo, Math.min(hi, Math.min(lay.lo[j], base))) - lo) / span) * bh - top;
      fill = lay.fill;
    }
    out.fill(lay.x0[j], top, lay.x1[j] - lay.x0[j], h, fill);
  }
}

const ValuesDraws = passDraws<ValuesLay, Emphasis | undefined>({ data: drawBars, hover: hoverCell });

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
  const n = axis.offsets.length;
  if (pl.direct) {
    for (let k = 0; k < n; k++) if (axis.spanIn(k, a, b)) mark(pl, axis);
    return;
  }
  toSection(pl.src!, a, b, pieces);
  for (let i = 0; i + 1 < pieces.length; i += 2)
    for (let k = 0; k < n; k++) if (axis.spanIn(k, pieces[i], pieces[i + 1])) mark(pl, axis);
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
