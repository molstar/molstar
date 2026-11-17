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
import type { PassName, PassWriter } from './draw-lists.js';
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
import { drawSpans, outlineCell } from './marks.js';
import { resolveColor, type SeqTokens, withAlpha } from './tokens.js';

// The `values` layer type (_spec/README.md "Render", Layer types): a bar a value, square (continuous data). Where a
// position is a device pixel or more, a bar a position; below, a bar a device column, its height the column's
// maximum (from per-position values or a coarser level's bins). Missing (NaN) is a thin neutral mark at the baseline,
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
  box: Box;
  axis: AxisTransform;
  step: DensityStep;
  tokens: SeqTokens;
  fill: number;
  domain: [number, number] | undefined;
  /** Bars: snapped x extents and values (NaN: missing). */
  n: number;
  x0: Float32Array;
  x1: Float32Array;
  v: Float32Array;
  pending: number[];
  failed: number[];
}

function grow(a: Float32Array, n: number): Float32Array {
  return n <= a.length ? a : new Float32Array(Math.max(n, 2 * a.length));
}

function push(lay: ValuesLay, x0: number, x1: number, v: number) {
  const j = lay.n++;
  if (j >= lay.x0.length) {
    lay.x0 = grow(lay.x0, j + 1);
    lay.x1 = grow(lay.x1, j + 1);
    lay.v = grow(lay.v, j + 1);
  }
  lay.x0[j] = x0;
  lay.x1[j] = x1;
  lay.v[j] = v;
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
      box,
      axis,
      step: 'cells',
      tokens,
      fill: 0,
      domain: undefined,
      n: 0,
      x0: new Float32Array(64),
      x1: new Float32Array(64),
      v: new Float32Array(64),
      pending: [],
      failed: [],
    };
    Object.assign(lay, { key: ctx.key, system: axis.system, box, axis, tokens });
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
    const step = densityStep(axis.devicePx, config.density);
    lay.step = step;
    const regions = dataRegions(src, axis);
    const level = perPosition(step) ? layerLevel(src, data, axis, true) : layerLevel(src, data, axis);
    const read = readTiles(ctx, data, src, regions, level);
    placeholderSpans(src, axis, read.pending, lay.pending);
    placeholderSpans(src, axis, read.failed, lay.failed);
    const pieces: number[] = [];
    const spans: number[] = [];
    const px = 1 / axis.dpr;
    // the value of data positions [a, b): every piece it lands on in view
    const place = (a: number, b: number, v: number, cols: Float32Array | undefined) => {
      toSection(src, a, b, pieces);
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        spans.length = 0;
        axis.spans(pieces[k], pieces[k + 1], spans);
        for (let s = 0; s + 2 < spans.length; s += 3) {
          if (!cols) {
            const x0 = axis.snap(spans[s]);
            push(lay, x0, Math.max(axis.snap(spans[s + 1]), x0 + px), v);
            continue;
          }
          const c0 = axis.column(spans[s]);
          const c1 = Math.max(c0, axis.column(spans[s + 1] - 1e-6));
          for (let c = c0; c <= c1; c++) if (Number.isNaN(cols[c]) || v > cols[c]) cols[c] = v;
        }
      }
    };
    const columns = perPosition(step) ? undefined : new Float32Array(axis.columns).fill(Number.NaN);
    for (const e of read.entries) {
      const t: ValuesTable | undefined = e.tables.values;
      if (!t) continue;
      if (t.form === 'positions') {
        for (let i = 0; i < t.value.length; i++) {
          const p = t.start + i;
          if (!inRegions(regions, p)) continue;
          place(p, p + 1, t.value[i], columns);
        }
      } else {
        const w = e.region.level;
        for (let i = 0; i < t.max.length; i++) {
          const a = t.start + i * w;
          if (!overlapsRegions(regions, a, a + w)) continue;
          place(a, a + w, t.count[i] > 0 ? t.max[i] : Number.NaN, columns);
        }
      }
    }
    if (columns) for (let c = 0; c < columns.length; c++) push(lay, box.x + c * px, box.x + (c + 1) * px, columns[c]);
    if (prep.domain) lay.domain = [prep.domain[0], prep.domain[1]];
    else if (ctx.settled || !lay.domain) lay.domain = openDomain(lay);
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    const { box, tokens } = lay;
    if (pass === 'data') {
      drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
      drawSpans(out, lay.failed, box.y, box.height, withAlpha(tokens.colors['seq-conflict'], 0.25));
      const [lo, hi] = lay.domain ?? [0, 1];
      const base = lo <= 0 && hi >= 0 ? 0 : lo;
      const span = hi - lo || 1;
      const yOf = (v: number) => box.y + box.height - ((Math.max(lo, Math.min(hi, v)) - lo) / span) * box.height;
      const yBase = yOf(base);
      const missing = tokens.colors['seq-missing'];
      for (let j = 0; j < lay.n; j++) {
        const v = lay.v[j];
        const w = lay.x1[j] - lay.x0[j];
        if (Number.isNaN(v)) {
          if (perPosition(lay.step)) out.rect(lay.x0[j], yBase - 2, w, 2, missing);
          continue;
        }
        const y = yOf(v);
        out.rect(lay.x0[j], Math.min(y, yBase), w, Math.abs(yBase - y), lay.fill);
      }
      return;
    }
    if (pass === 'hover') {
      const h = emphasis?.hover;
      if (!h || h.layer !== lay.key || h.position === undefined) return;
      outlineCell(out, lay.axis, h.position, box, tokens.colors['seq-hover'], 1.5);
    }
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height) return undefined;
    const position = lay.axis.positionAt(x);
    if (position < 0) return undefined;
    return { layer: lay.key, position, system: lay.system };
  },

  legend(lay): LegendSpec | undefined {
    return lay.domain ? { kind: 'ramp', domain: lay.domain, colors: [lay.fill] } : undefined;
  },
};

function inRegions(regions: readonly number[], p: number): boolean {
  for (let i = 0; i + 1 < regions.length; i += 2) if (p >= regions[i] && p < regions[i + 1]) return true;
  return false;
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
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  if (lo > hi) return [0, 1];
  if (lo === hi) return lo >= 0 ? [0, hi || 1] : [lo, 0];
  return [Math.min(lo, 0), hi];
}
