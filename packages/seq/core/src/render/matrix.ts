/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The matrix layer (_spec/README.md "Render", The matrix layer): items × states, states whole (EC-220 to EC-222). The
// states are rows in the band, scrolled by the band's `bandTop` and zoomed by its `rowScale` (the rows layer's band
// state; fit to the band when absent); the items' columns over the strip's x. At the source's finest level each item
// is a rect a state row at its own extent, a long one solid across the columns it covers (drawn first: the short ones
// over it); at a coarser level the bins are one image tile a store tile (bins × states, never smoothed) in the shared
// image cache. Values are coloured by the layer's `color` channel by value (its open domain the range of the tiles
// read), codes by the genotype classes' tokens; a missing value transparent, an uncalled code `seq-missing`. A tile not
// loaded is pending across the band, never empty (EC-169); the layer's wants are bounded together (`readTiles`'s key).

import type { Node, Value } from '@molstar/core/repr/tree';
import { type MatrixBins, type MatrixItems, type MatrixTable, tableRows } from '@molstar/seqx-core/data';
import { theme as themeForm } from '@molstar/seqx-core/lang';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { Channel, ChannelName } from '@molstar/seqx-core/theme/index';
import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { channelForms, channelLegends, Extents, type LayerChannels, PlainDatum } from './channels.js';
import { grown, type ImageData, NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { colorOf } from './grid.js';
import { ImageCache } from './images.js';
import {
  dataRegions,
  fittedLevel,
  type LayerSource,
  layerLevel,
  layerSource,
  placeholderSpans,
  readTiles,
  toSection,
} from './layer-data.js';
import { type Box, type Hit, type LayerType, type LayoutContext, sectionId, type StoreReader } from './layer-type.js';
import type { SeqTokens } from './tokens.js';

const ValueEnv = { numeric: (n: string) => n === 'value' };
const DefaultColor = themeForm.field('value', { list: false });

/** A genotype class's token (`0` no copy, `1` some, `2` every copy, `3` partly called); `255` uncalled. */
const CodeTokens = ['seq-ramp-0', 'seq-ramp-2', 'seq-ramp-4', 'seq-ramp-1'] as const;
const Uncalled = 255;

/** A row's height when the band does not say and the states do not fit it: a state at least this many px. */
export const MatrixMinRowPx = 1;

export interface MatrixPrep {
  readonly owner: string;
  readonly src?: LayerSource;
  /** The states system's address (`Hit.ySystem`). */
  readonly statesSystem: string;
  /** The states' count when the spec says it (inline `states`, or the source's coverage at resolve). */
  readonly states?: number;
  /** The states' names given inline. */
  readonly names?: readonly string[];
  readonly channels: LayerChannels;
}

type ScratchHit = { -readonly [K in keyof Hit]: Hit[K] };

export class MatrixLay {
  key = '';
  section = '';
  system = '';
  statesSystem = '';
  box: Box = { x: 0, y: 0, width: 0, height: 0 };
  axis: AxisTransform | undefined;
  tokens: SeqTokens | undefined;
  color: Channel | undefined;
  colorVersion = 0;
  readonly datum = new PlainDatum();
  readonly extents = new Extents();
  /** The states' count, px a state row, the scroll in px, the first and past the last state drawn. */
  states = 0;
  rowScale = 0;
  bandTop = 0;
  first = 0;
  last = 0;
  /** The level read (`Hit.bin` past the finest), and whether it read items (the finest) or bins. */
  level = 1;
  items = false;
  /** Rects (x, y, w, h each) and their fills. */
  n = 0;
  rx = new Float64Array(64);
  rc = new Uint32Array(16);
  /** Bin tiles: their boxes (x, y, w, h each) and images. */
  readonly images: ImageData[] = [];
  readonly ib: number[] = [];
  cache: ImageCache | undefined;
  /** The tables read (a hover's cell: `matrixCellAt`) and the states' names (inline, else the coverage's). */
  readonly tables: MatrixTable[] = [];
  names: readonly string[] = NoNames;
  /** Pending spans on x (x0, x1, cut triples): the pending fill across the band. */
  readonly pend: number[] = [];
  readonly scratch: ScratchHit = { layer: '', section: '', system: '' };
}

const OwnImages = new WeakMap<object, ImageCache>();
const NoNames: readonly string[] = [];

function prepare(layer: ResolvedLayer, model: StaticModel): MatrixPrep {
  const owner = layer.data?.address ?? '';
  const a = model.annotations.get(owner);
  const statesSystem = a?.shape === 'matrix' ? (a.states ?? '') : '';
  const states = model.systems.get(statesSystem)?.length;
  const src = a?.shape === 'matrix' ? layerSource(model, layer) : undefined;
  const given = names(nodeAt(model.spec, owner)?.params?.states);
  return {
    owner,
    ...(src ? { src } : {}),
    statesSystem,
    ...(states !== undefined ? { states } : {}),
    ...(given ? { names: given } : {}),
    channels: channelForms(layer.theme),
  };
}

/** A node of the spec by its exact address (`/data/gt`). */
function nodeAt(spec: Node, address: string): Node | undefined {
  let n: Node | undefined = spec;
  for (const id of address.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === id);
  return n;
}

function names(v: Value | undefined): readonly string[] | undefined {
  return Array.isArray(v) && v.every((x) => typeof x === 'string') ? (v as string[]) : undefined;
}

/** The value range of a table's cells (values or bins' means), kept by the table. */
const Ranges = new WeakMap<MatrixTable, readonly [number, number]>();

function rangeOf(t: MatrixTable): readonly [number, number] {
  let r = Ranges.get(t);
  if (r) return r;
  const v = t.form === 'items' ? t.value : t.mean;
  let lo = Number.NaN;
  let hi = Number.NaN;
  if (v)
    for (let i = 0; i < v.length; i++) {
      const x = v[i];
      if (Number.isNaN(x)) continue;
      if (!(x >= lo)) lo = x;
      if (!(x <= hi)) hi = x;
    }
  Ranges.set(t, (r = [lo, hi]));
  return r;
}

function codeColor(lay: MatrixLay, code: number): number {
  if (code === Uncalled) return lay.tokens!.colors['seq-missing'];
  const t = CodeTokens[code];
  return t ? lay.tokens!.colors[t] : NoColor;
}

/** A cell's fill: state `s` of column `i` (an item or a bin). */
function cellColor(lay: MatrixLay, t: MatrixTable, i: number, s: number): number {
  const k = i * t.states + s;
  if (t.form === 'items') return t.codes ? codeColor(lay, t.codes[k]) : t.value ? colorOf(lay, t.value[k]) : NoColor;
  if (t.classes) {
    // a bin's majority class, its share the colour's strength (design: "its majority with its share")
    const c = codeColor(lay, t.classes[k]);
    if (c === NoColor || !t.shares) return c;
    return ((c & 0xffffff00) | Math.round(((c & 0xff) * t.shares[k]) / 255)) >>> 0;
  }
  return t.mean ? colorOf(lay, t.mean[k]) : NoColor;
}

function pushRect(lay: MatrixLay, x: number, y: number, w: number, h: number, c: number) {
  const k = lay.n++;
  lay.rx = grown(lay.rx, 4 * k + 4);
  lay.rc = grown(lay.rc, k + 1);
  lay.rx[4 * k] = x;
  lay.rx[4 * k + 1] = y;
  lay.rx[4 * k + 2] = w;
  lay.rx[4 * k + 3] = h;
  lay.rc[k] = c;
}

const Pieces: number[] = [];
const Spans: number[] = [];

/**
 * Items' rects over the state rows drawn: the long ones (more than a position) first, the short ones over them. An items
 * tile holds every item with an end in it, uncut: one starting before `from` (an earlier loaded tile's end) is that
 * tile's, drawn once.
 */
function layoutItems(lay: MatrixLay, src: LayerSource, t: MatrixItems, a: AxisTransform, from: number) {
  const minW = 1 / a.dpr;
  for (let pass = 0; pass < 2; pass++)
    for (let i = 0; i < t.start.length; i++) {
      const s0 = t.start[i];
      const e0 = Math.max(t.end[i], s0 + 1);
      if (s0 < from || e0 - s0 > 1 !== (pass === 0)) continue;
      Pieces.length = 0;
      toSection(src, s0, e0, Pieces);
      for (let p = 0; p + 1 < Pieces.length; p += 2) {
        Spans.length = 0;
        a.spans(Pieces[p], Pieces[p + 1], Spans);
        for (let q = 0; q + 2 < Spans.length; q += 3) {
          const x0 = Math.min(Spans[q], Spans[q + 1]);
          const w = Math.max(minW, Math.abs(Spans[q + 1] - Spans[q]));
          for (let s = lay.first; s < lay.last; s++) {
            const c = cellColor(lay, t, i, s);
            if (c === NoColor) continue;
            pushRect(lay, x0, rowY(lay, s), w, lay.rowScale, c);
          }
        }
      }
    }
}

/** A bins tile's image: a column a bin, a row a state, in display order (a reversed axis's bins right to left). */
function buildBins(lay: MatrixLay, t: MatrixBins, n: number): ImageData | undefined {
  if (n <= 0 || t.states <= 0) return undefined;
  const img = lay.cache!.slot(n, t.states);
  const d = img.data;
  const flip = lay.axis!.dir === -1;
  for (let s = 0; s < t.states; s++)
    for (let i = 0; i < n; i++) {
      const q = 4 * (s * n + (flip ? n - 1 - i : i));
      const fill = cellColor(lay, t, i, s);
      if (fill === NoColor) {
        d[q] = d[q + 1] = d[q + 2] = d[q + 3] = 0;
        continue;
      }
      d[q] = fill >>> 24;
      d[q + 1] = (fill >>> 16) & 0xff;
      d[q + 2] = (fill >>> 8) & 0xff;
      d[q + 3] = fill & 0xff;
    }
  return img;
}

/** A state row's top in the strip's frame. */
function rowY(lay: MatrixLay, s: number): number {
  return lay.box.y + s * lay.rowScale - lay.bandTop;
}

export const MatrixType: LayerType<MatrixPrep, MatrixLay> = {
  name: 'matrix',
  shape: 'matrix',
  fallback: 'values',
  axes: 1,
  prepare,

  layout(prep, ctx: LayoutContext, data: StoreReader, _overlay, prev?: MatrixLay): MatrixLay {
    const { axis, box, config } = ctx;
    const lay = prev ?? new MatrixLay();
    lay.key = ctx.key;
    lay.section = sectionId(ctx.section);
    lay.system = axis.system;
    lay.statesSystem = prep.statesSystem;
    lay.box = box;
    lay.axis = axis;
    lay.tokens = ctx.tokens;
    lay.n = 0;
    lay.images.length = 0;
    lay.ib.length = 0;
    lay.pend.length = 0;
    lay.tables.length = 0;
    if (ctx.images) lay.cache = ctx.images;
    else if (!lay.cache) {
      let own = OwnImages.get(config);
      if (!own) OwnImages.set(config, (own = new ImageCache(config.imageBytes)));
      lay.cache = own;
    }
    const src = prep.src;
    if (!src) return lay;
    const color = prep.channels.get('color' as ChannelName, DefaultColor, ctx, ValueEnv);
    if (color !== lay.color) {
      lay.color = color;
      lay.colorVersion++;
    }
    lay.extents.reset([color]);
    const levels = data.coverage(src.source)?.levels ?? [1];
    const regions = dataRegions(src, axis);
    // past the store's bounds on the layer's wants together: coarser, then zoom in
    lay.level = fittedLevel(ctx, data, src, regions, layerLevel(src, data, axis));
    lay.items = lay.level === levels[0];
    const read = readTiles(ctx, data, src, regions, lay.level, { layer: ctx.key });
    // the states: the spec's count, else the tables' (a source's coverage not known at resolve)
    let states = prep.states ?? 0;
    for (const e of read.entries)
      if (e.tables.matrix && !prep.states) states = Math.max(states, e.tables.matrix.states);
    lay.states = states;
    lay.names =
      prep.names ??
      data.coverage(src.source)?.shapes.find((x) => x.shape === 'matrix' && x.system === src.system)?.states ??
      NoNames;
    // the band's rows: its zoom (px a state) when it says, else the states fit to the band; its scroll clamped
    const rs = ctx.rowScale ?? Math.max(MatrixMinRowPx, box.height / Math.max(1, states));
    lay.rowScale = rs;
    lay.bandTop = Math.max(0, Math.min(ctx.bandTop ?? 0, states * rs - box.height));
    lay.first = Math.max(0, Math.floor(lay.bandTop / rs));
    lay.last = Math.min(states, Math.ceil((lay.bandTop + box.height) / rs));
    // the color's open domain: the range of the tiles read
    if (lay.extents.any) {
      let lo = Number.NaN;
      let hi = Number.NaN;
      for (const e of read.entries) {
        const t = e.tables.matrix;
        if (!t) continue;
        const [l, h] = rangeOf(t);
        if (!(l >= lo)) lo = l;
        if (!(h <= hi)) hi = h;
      }
      if (!Number.isNaN(lo)) {
        lay.datum.value = lo;
        lay.extents.add(lay.datum);
        lay.datum.value = hi;
        lay.extents.add(lay.datum);
        lay.extents.flush(ctx.theme);
      }
    }
    if (read.pending.length) placeholderSpans(src, axis, read.pending, lay.pend);
    const tokens = ctx.tokens.version;
    for (const e of read.entries) {
      const t = e.tables.matrix;
      if (!t) continue;
      lay.tables.push(t);
      if (t.form === 'items') {
        let from = Number.NEGATIVE_INFINITY;
        for (const o of read.entries) if (o.region.start < e.region.start && o.region.end > from) from = o.region.end;
        layoutItems(lay, src, t, axis, from);
        continue;
      }
      const n = tableRows(t);
      const key = `${lay.key}|m|${e.region.source}|${e.region.system}|${e.region.level}|${e.region.tile}|${e.version}|${lay.colorVersion}.${tokens}|${axis.dir}`;
      let img = lay.cache.get(key);
      if (!img) {
        img = buildBins(lay, t, n);
        if (!img) continue;
        lay.cache.put(key, img);
      }
      Pieces.length = 0;
      toSection(src, t.start, t.start + n * lay.level, Pieces);
      if (Pieces.length !== 2) continue;
      const d = axis.dpr;
      const x0 = Math.round(axis.x(Pieces[0]) * d) / d;
      const x1 = Math.round(axis.x(Pieces[1]) * d) / d;
      if (Number.isNaN(x0) || Number.isNaN(x1)) continue;
      lay.images.push(img);
      lay.ib.push(Math.min(x0, x1), rowY(lay, 0), Math.abs(x1 - x0), states * rs);
    }
    return lay;
  },

  draw(lay, pass: PassName, out: PassWriter) {
    if (pass !== 'data' || !lay.tokens || !lay.axis) return;
    const b = lay.box;
    // the states scrolled: what is past the band is cut
    const g = out.args;
    g[0] = b.x;
    g[1] = b.y;
    g[2] = b.width;
    g[3] = b.height;
    // through the args: numbers passed to a call that is not inlined are boxed (on a ring's warp writer)
    out.clipArgs();
    if (lay.pend.length) {
      const pending = lay.tokens.colors['seq-pending'];
      for (let k = 0; k + 2 < lay.pend.length; k += 3) {
        const x0 = Math.min(lay.pend[k], lay.pend[k + 1]);
        out.fill(x0, b.y, Math.abs(lay.pend[k + 1] - lay.pend[k]), b.height, pending);
      }
    }
    for (let k = 0; k < lay.images.length; k++) {
      g[0] = lay.ib[4 * k];
      g[1] = lay.ib[4 * k + 1];
      g[2] = lay.ib[4 * k + 2];
      g[3] = lay.ib[4 * k + 3];
      out.imageArgs(lay.images[k]);
    }
    for (let k = 0; k < lay.n; k++)
      out.fill(lay.rx[4 * k], lay.rx[4 * k + 1], lay.rx[4 * k + 2], lay.rx[4 * k + 3], lay.rc[k]);
    out.unclip();
  },

  hit(lay, at) {
    const a = lay.axis;
    if (!a || !(lay.rowScale > 0) || !lay.states) return undefined;
    const b = lay.box;
    if (at[1] < b.y || at[1] >= b.y + b.height) return undefined;
    const s = Math.floor((at[1] - b.y + lay.bandTop) / lay.rowScale);
    const i = a.positionOf(at, 0);
    if (i < 0 || s < 0 || s >= lay.states) return undefined;
    const h = lay.scratch;
    h.layer = lay.key;
    h.section = lay.section;
    h.system = lay.system;
    h.position = i;
    h.ySystem = lay.statesSystem;
    h.yPosition = s;
    h.target = 'cell';
    h.bin = lay.items || lay.level <= 1 ? undefined : lay.level;
    return h;
  },

  legend(lay, options) {
    return channelLegends([lay.color], { ...options, layer: lay.key });
  },
};
