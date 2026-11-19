/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArrayRange } from '@molstar/core/data/util/sort';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { type DensityStep, densityStep, perPosition } from './config.js';
import { drawPass, grown, NoColor, type PassName, type PassWriter, passDraws } from './draw-lists.js';
import { collectItems, countStrip, itemHit, type ItemRows, newItemRows } from './items.js';
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
import { reportColor } from './letters.js';
import { drawFailed, drawSpans } from './marks.js';
import type { StripBuffers } from './ranges.js';
import { categoryToken, resolveColor, type SeqTokens } from './tokens.js';

// The `variants` layer type (_spec/README.md "Render", Layer types): a site a stem from the track's floor and a round
// head, the heads of one site stacked up the stem (squeezed together when they do not fit). At the `bars` and `density`
// steps (and from a source's density bins) a strip counts the variants a device column, in the density tokens. A hit is
// a head, or its stem.

/** The head's radius and the stem's least height under it, CSS px. */
const Radius = 4;
const Stem = 6;
const Sides = 16;

export interface VariantsPrep {
  readonly src?: LayerSource;
  readonly owner: string;
  readonly color?: unknown;
}

export interface VariantsLay {
  key: string;
  system: string;
  section: string;
  owner: string;
  dataSystem: string;
  src?: LayerSource;
  box: Box;
  axis: AxisTransform;
  step: DensityStep;
  tokens: SeqTokens;
  rows: ItemRows;
  fixed: number | undefined;
  /** Marks in x order: the head's center x and y, the row, the fill. */
  m: number;
  mx: Float32Array;
  my: Float32Array;
  mrow: Int32Array;
  mfill: Uint32Array;
  /** Marks by row in the order they were made (sorted into the above). */
  px: Float64Array;
  prow: Int32Array;
  pn: number;
  order: Int32Array;
  density?: StripBuffers['image'];
  strip: StripBuffers;
  pending: number[];
  failed: number[];
  hoverWidth: number;
  hoverOffset: number;
}

const itemPiece: number[] = [0, 0];
const mapped: number[] = [];
const ring = new Float64Array(2 * Sides);
const Cos = Float64Array.from({ length: Sides }, (_, i) => Math.cos((2 * Math.PI * i) / Sides));
const Sin = Float64Array.from({ length: Sides }, (_, i) => Math.sin((2 * Math.PI * i) / Sides));

function nextPoint(lay: VariantsLay): number {
  const j = lay.pn++;
  if (j >= lay.px.length) {
    lay.px = grown(lay.px, j + 1);
    lay.prow = grown(lay.prow, j + 1);
  }
  return j;
}

export const VariantsType: LayerType<VariantsPrep, VariantsLay> = {
  name: 'variants',
  fallback: 'ranges',

  prepare(layer: ResolvedLayer, model: StaticModel): VariantsPrep {
    const src = layerSource(model, layer);
    return {
      ...(src ? { src } : {}),
      owner: layer.data?.address ?? '',
      ...(layer.theme.color !== undefined ? { color: layer.theme.color } : {}),
    };
  },

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay: Overlay, prev?: VariantsLay): VariantsLay {
    const { axis, box, config, tokens } = ctx;
    const lay: VariantsLay = prev ?? {
      key: ctx.key,
      system: axis.system,
      section: sectionId(ctx.section),
      owner: prep.owner,
      dataSystem: '',
      box,
      axis,
      step: 'cells',
      tokens,
      rows: newItemRows(),
      fixed: undefined,
      m: 0,
      mx: new Float32Array(64),
      my: new Float32Array(64),
      mrow: new Int32Array(64),
      mfill: new Uint32Array(64),
      px: new Float64Array(64),
      prow: new Int32Array(64),
      pn: 0,
      order: new Int32Array(64),
      strip: { diff: new Float64Array(0), counts: new Float64Array(0), v: 1 },
      pending: [],
      failed: [],
      hoverWidth: config.outline.hover,
      hoverOffset: config.emphasis.hoverOffset,
    };
    lay.key = ctx.key;
    lay.system = axis.system;
    lay.section = sectionId(ctx.section);
    lay.owner = prep.owner;
    lay.dataSystem = prep.src?.dataSystem ?? '';
    lay.src = prep.src;
    lay.box = box;
    lay.axis = axis;
    lay.tokens = tokens;
    lay.hoverWidth = config.outline.hover;
    lay.hoverOffset = config.emphasis.hoverOffset;
    lay.m = 0;
    lay.pn = 0;
    lay.rows.n = 0;
    lay.pending.length = 0;
    lay.failed.length = 0;
    lay.density = undefined;
    const src = prep.src;
    if (!src) return lay;
    lay.fixed = undefined;
    if (prep.color !== undefined) {
      lay.fixed = resolveColor(prep.color, tokens);
      if (lay.fixed === undefined) reportColor(ctx);
    }
    const step = densityStep(axis.scale, axis.dpr, config);
    lay.step = step;
    const regions = dataRegions(src, axis);
    const read = readTiles(ctx, data, src, regions, layerLevel(src, data, axis));
    placeholderSpans(src, axis, read.pending, lay.pending);
    placeholderSpans(src, axis, read.failed, lay.failed);
    const rows = lay.rows;
    const densities = collectItems(rows, read.entries, ['consequence', 'category']);
    if (!perPosition(step) || densities.length) {
      lay.density = countStrip(lay.strip, src, axis, rows, densities, tokens);
      return lay;
    }
    // a head at the center of each piece of a site in view
    for (let i = 0; i < rows.n; i++) {
      let pieces = itemPiece;
      if (src.path) pieces = toSection(src, rows.starts[i], rows.ends[i], mapped);
      else {
        itemPiece[0] = rows.starts[i];
        itemPiece[1] = Math.max(rows.ends[i], rows.starts[i] + 1);
      }
      for (let k = 0; k + 1 < pieces.length; k += 2) {
        axis.piecesIn(pieces[k], pieces[k + 1]);
        const ka = axis.ka;
        const kb = axis.kb;
        for (let a = ka; a <= kb; a++) {
          if (!axis.spanIn(a, pieces[k], pieces[k + 1])) continue;
          const j = nextPoint(lay);
          lay.px[j] = axis.snap((axis.sx0 + axis.sx1) / 2);
          lay.prow[j] = i;
        }
      }
    }
    stack(lay);
    return lay;
  },

  draw(lay, _pass: PassName, out: PassWriter, emphasis?: Emphasis) {
    drawPass(VariantsDraws, lay, out, emphasis);
  },

  hit(lay, x, y): Hit | undefined {
    const { box } = lay;
    if (y < box.y || y >= box.y + box.height || x < box.x || x >= box.x + box.width) return undefined;
    const position = lay.axis.positionAt(x);
    const base = { layer: lay.key, section: lay.section, system: lay.system, ...(position >= 0 ? { position } : {}) };
    if (lay.density || !lay.m) return base;
    // the marks by x: those within a head's width of the pointer; the head under it, else a stem it is on
    let lo = 0;
    let hi = lay.m;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (lay.mx[mid] < x - Radius - 1) lo = mid + 1;
      else hi = mid;
    }
    let best = -1;
    let dist = Number.POSITIVE_INFINITY;
    const floor = box.y + box.height;
    for (let j = lo; j < lay.m && lay.mx[j] <= x + Radius + 1; j++) {
      const d = Math.hypot(x - lay.mx[j], y - lay.my[j]);
      const onHead = d <= Radius;
      const onStem = Math.abs(x - lay.mx[j]) <= 1.5 && y >= lay.my[j] && y <= floor;
      if ((onHead || onStem) && d < dist) {
        best = j;
        dist = d;
      }
    }
    if (best < 0) return base;
    const i = lay.mrow[best];
    return itemHit(base, lay.owner, lay.dataSystem, lay.src, lay.rows.keys[i], lay.rows.starts[i], lay.rows.ends[i]);
  },
};

/**
 * The marks in x order, the heads of one site stacked up its stem: a head's center y, the spacing squeezed so the
 * stack fits the track.
 */
function stack(lay: VariantsLay) {
  const n = lay.pn;
  lay.order = grown(lay.order, n);
  const order = lay.order.subarray(0, n);
  for (let i = 0; i < n; i++) order[i] = i;
  const { px, prow } = lay;
  sortArrayRange(order, 0, n, (o, a, b) => px[o[a]] - px[o[b]] || prow[o[a]] - prow[o[b]]);
  if (lay.mx.length < n) {
    const size = Math.max(n, 2 * lay.mx.length);
    lay.mx = new Float32Array(size);
    lay.my = new Float32Array(size);
    lay.mrow = new Int32Array(size);
    lay.mfill = new Uint32Array(size);
  }
  const { box, tokens } = lay;
  const floor = box.y + box.height - 1;
  const room = box.height - 1 - Stem - 2 * Radius;
  for (let s = 0; s < n; ) {
    let e = s + 1;
    while (e < n && px[order[e]] === px[order[s]]) e++;
    const gap = e - s > 1 ? Math.min(2 * Radius + 1, Math.max(1, room) / (e - s - 1)) : 0;
    for (let k = s; k < e; k++) {
      const o = order[k];
      lay.mx[k] = px[o];
      lay.my[k] = floor - Stem - Radius - (k - s) * gap;
      lay.mrow[k] = prow[o];
      const v = lay.rows.values[prow[o]];
      lay.mfill[k] = lay.fixed ?? tokens.colors[v === undefined ? 'seq-conflict' : categoryToken(v)];
    }
    s = e;
  }
  lay.m = n;
}

/** The head: a polygon of `Sides` corners around (x, y) at radius r, written to `ring` (r grown by the caller's offset). */
function head(x: number, y: number, r: number) {
  for (let i = 0; i < Sides; i++) {
    ring[2 * i] = x + r * Cos[i];
    ring[2 * i + 1] = y + r * Sin[i];
  }
}

// Each pass draws in a function of its own, called through a table (`drawPass`: nothing allocated once warm).

/** The data pass: placeholders, the strip, or the stems and heads (a stem under every head, heads over stems). */
function drawVariants(lay: VariantsLay, out: PassWriter, _emphasis: Emphasis | undefined) {
  const { box, tokens } = lay;
  drawSpans(out, lay.pending, box.y, box.height, tokens.colors['seq-pending']);
  drawFailed(out, lay.failed, box.y, box.height, tokens);
  const a = out.args;
  if (lay.density) {
    a[0] = box.x;
    a[1] = box.y;
    a[2] = box.width;
    a[3] = box.height;
    out.imageArgs(lay.density);
    return;
  }
  const px = 1 / lay.axis.dpr;
  const floor = box.y + box.height;
  for (let j = 0; j < lay.m; j++) {
    a[0] = lay.mx[j] - px / 2;
    a[1] = lay.my[j];
    a[2] = px;
    a[3] = floor - lay.my[j];
    a[4] = lay.mfill[j];
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    a[8] = 0;
    out.rectArgs();
  }
  for (let j = 0; j < lay.m; j++) {
    head(lay.mx[j], lay.my[j], Radius);
    a[4] = lay.mfill[j];
    a[5] = NoColor;
    a[6] = 0;
    a[7] = 0;
    out.pathArgs(ring, Sides, true);
  }
}

/** The hover pass: an outline around the hovered site's heads, outside them by the offset. */
function drawHover(lay: VariantsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key || !h.element) return;
  const r = Radius + lay.hoverWidth / 2 + lay.hoverOffset;
  const a = out.args;
  for (let j = 0; j < lay.m; j++) {
    if (lay.rows.keys[lay.mrow[j]] !== h.element.key) continue;
    head(lay.mx[j], lay.my[j], r);
    a[4] = NoColor;
    a[5] = lay.tokens.colors['seq-hover'];
    a[6] = lay.hoverWidth;
    a[7] = 0;
    out.pathArgs(ring, Sides, true);
  }
}

const VariantsDraws = passDraws<VariantsLay, Emphasis | undefined>({ data: drawVariants, hover: drawHover });
