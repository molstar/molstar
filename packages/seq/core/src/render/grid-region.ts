/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GridTable } from '@molstar/seqx-core/data';
import type { AxisTransform } from './axis.js';
import { type ImageData, NoColor } from './draw-lists.js';
import { colorOf, type GridLay, type GridPrep, keepsCell } from './grid.js';
import { dataRegions, fittedLevel, gridLevel, type LayerSource, readGridTiles } from './layer-data.js';
import type { LayoutContext, StoreReader } from './layer-type.js';

// A grid read by region (_spec/README.md "Render", The grid layer: by region): a `grid` node with a `source` is read
// through the store in 2D tiles at the finest level whose bin is at least `grid.binPx` device pixels, a ring of
// `grid.marginTiles` round the view. Each store tile is one image tile (its cells coloured by the value theme, a sparse
// tile's `fill` honoured, NaN transparent, a symmetric grid's mirror tile transposed), kept in the shared image cache
// by (layer, the entry's region and version, the colours' version, the triangle, the mirror). A tile not loaded draws
// the coarser held tile under it; with none held, a pending rect (the pending fill), never empty (EC-169).

/** Where a region grid reads: its source on x, and its y system's name as the source knows it. */
export interface GridRegion {
  readonly src: LayerSource;
  readonly yName: string;
}

/** The value range of a table, kept by the table (a tile's colours are held to the range of the tiles read). */
const Ranges = new WeakMap<GridTable, readonly [number, number]>();

function rangeOf(t: GridTable): readonly [number, number] {
  let r = Ranges.get(t);
  if (r) return r;
  let lo = Number.NaN;
  let hi = Number.NaN;
  const v = t.value;
  for (let i = 0; i < v.length; i++) {
    const x = v[i];
    if (Number.isNaN(x)) continue;
    if (!(x >= lo)) lo = x;
    if (!(x <= hi)) hi = x;
  }
  if (t.form === 'sparse' && !Number.isNaN(t.fill)) {
    if (!(t.fill >= lo)) lo = t.fill;
    if (!(t.fill <= hi)) hi = t.fill;
  }
  Ranges.set(t, (r = [lo, hi]));
  return r;
}

/** A tile read: its table, its bins' box (x, y, width, height in bins, as drawn: a mirror's swapped), its key. */
export interface Placed {
  readonly table: GridTable;
  readonly transposed: boolean;
  readonly x: number;
  readonly y: number;
  readonly w: number;
  readonly h: number;
  readonly key: string;
}

/** The view's x range as regions (`fittedLevel`). */
const XRange = [0, 0];

/** A reader that asks for nothing (what a coarser level holds is read, never wanted). */
const NoWants = { want: () => {} };

/** The tiles loaded over the view at `level`, as drawn (no want: what a coarser level holds under a pending one). */
function loadedAt(
  ctx: LayoutContext,
  data: StoreReader,
  region: GridRegion,
  x: readonly [number, number],
  y: { readonly system: string; readonly start: number; readonly end: number },
  level: number,
  want: boolean,
  marginTiles: number,
) {
  const read = readGridTiles(want ? ctx : NoWants, data, region.src, x, y, level, want ? marginTiles : 0);
  const out: Placed[] = [];
  for (const { entry, transposed } of read.entries) {
    const t = entry.tables.grid;
    const r = entry.region;
    if (!t || !r.y) continue;
    // the served tile's bins: its region's at the level
    const bx = Math.floor(r.start / level);
    const by = Math.floor(r.y.start / level);
    const bw = Math.ceil(r.end / level) - bx;
    const bh = Math.ceil(r.y.end / level) - by;
    const key = `${r.system}|${r.tile}|${r.y.system}|${r.y.tile}|${level}|${entry.version}`;
    out.push(
      transposed
        ? { table: t, transposed, x: by, y: bx, w: bh, h: bw, key }
        : { table: t, transposed, x: bx, y: by, w: bw, h: bh, key },
    );
  }
  return { tiles: out, pending: read.pending, tooMany: read.tooMany };
}

/** A tile's cells as an image, in display order (a reversed axis's right to left), outside the triangle transparent. */
function buildTile(lay: GridLay, p: Placed, level: number): ImageData | undefined {
  if (p.w <= 0 || p.h <= 0) return undefined;
  const img = lay.cache!.slot(p.w, p.h);
  const d = img.data;
  const t = p.table;
  // the served table's cells by its own (x, y) bins: a mirror reads them swapped
  const value = (i: number, j: number) => {
    const sx = p.transposed ? p.y + j : p.x + i;
    const sy = p.transposed ? p.x + i : p.y + j;
    if (t.form === 'dense') {
      const a = sx - t.x;
      const b = sy - t.y;
      return a < 0 || b < 0 || a >= t.nx || b >= t.ny ? Number.NaN : t.value[a * t.ny + b];
    }
    return Number.NaN;
  };
  // a sparse tile: its fill, then its listed cells
  let sparse: Float32Array | undefined;
  if (t.form === 'sparse') {
    sparse = new Float32Array(p.w * p.h).fill(t.fill);
    for (let k = 0; k < t.value.length; k++) {
      const i = (p.transposed ? t.yBin[k] : t.xBin[k]) - p.x;
      const j = (p.transposed ? t.xBin[k] : t.yBin[k]) - p.y;
      if (i >= 0 && j >= 0 && i < p.w && j < p.h) sparse[i * p.h + j] = t.value[k];
    }
  }
  const flipX = lay.axis!.dir === -1;
  const flipY = lay.axisY!.dir === -1;
  for (let j = 0; j < p.h; j++) {
    const row = flipY ? p.h - 1 - j : j;
    for (let i = 0; i < p.w; i++) {
      const q = 4 * (row * p.w + (flipX ? p.w - 1 - i : i));
      const v = keepsCell(lay, p.x + i, p.y + j) ? (sparse ? sparse[i * p.h + j] : value(i, j)) : Number.NaN;
      const fill = colorOf(lay, v);
      if (fill === NoColor) {
        d[q] = d[q + 1] = d[q + 2] = d[q + 3] = 0;
        continue;
      }
      d[q] = fill >>> 24;
      d[q + 1] = (fill >>> 16) & 0xff;
      d[q + 2] = (fill >>> 8) & 0xff;
      d[q + 3] = fill & 0xff;
    }
  }
  return img;
}

/** A box of positions on the plot, snapped to the device grid (tiles meeting at a fraction would seam), into `out`. */
function boxOf(a: AxisTransform, ay: AxisTransform, xs: number, xe: number, ys: number, ye: number, out: number[]) {
  const d = a.dpr;
  const x0 = Math.round(a.x(xs) * d) / d;
  const x1 = Math.round(a.x(xe) * d) / d;
  const y0 = Math.round(ay.x(ys) * d) / d;
  const y1 = Math.round(ay.x(ye) * d) / d;
  if (Number.isNaN(x0) || Number.isNaN(x1) || Number.isNaN(y0) || Number.isNaN(y1)) return false;
  out.push(Math.min(x0, x1), Math.min(y0, y1), Math.abs(x1 - x0), Math.abs(y1 - y0));
  return true;
}

/**
 * Lays a region grid out: its level, its tiles' images (the coarser held ones first, under the pending), pending rects
 * where nothing is held; its colours held to the range of the tiles read. `lay.level` is the level read (`Hit.bin`).
 */
export function layoutRegionGrid(
  lay: GridLay,
  prep: GridPrep,
  region: GridRegion,
  ctx: LayoutContext,
  data: StoreReader,
  a: AxisTransform,
  ay: AxisTransform,
): void {
  const cov = data.coverage(region.src.source);
  const levels = cov?.levels ?? [1];
  const xr = dataRegions(region.src, a);
  const yr = dataRegions({ ...region.src, system: region.yName, path: undefined, back: undefined }, ay);
  if (!xr.length || !yr.length) return;
  const x: [number, number] = [xr[0], xr[xr.length - 1]];
  const y = { system: region.yName, start: yr[0], end: yr[yr.length - 1] };
  const cfg = ctx.config.grid;
  const fit = gridLevel(
    levels,
    { x: 1 / a.devicePx, y: 1 / ay.devicePx },
    { x: x[1] - x[0], y: y.end - y.start },
    data,
    cfg.binPx,
  );
  // past the store's bounds on the layer's wants together: coarser, then zoom in
  XRange[0] = x[0];
  XRange[1] = x[1];
  const L = fittedLevel(ctx, data, region.src, XRange, fit, y);
  lay.binLevel = L;
  const now = loadedAt(ctx, data, region, x, y, L, true, cfg.marginTiles);
  let under: Placed[] = [];
  let underLevel = L;
  if (now.pending.length)
    for (const c of levels) {
      if (c <= L) continue;
      const held = loadedAt(ctx, data, region, x, y, c, false, 0).tiles;
      if (held.length) {
        under = held;
        underLevel = c;
        break;
      }
    }
  // the colours' domain: the range of the tiles drawn (open domains held to it)
  let lo = Number.NaN;
  let hi = Number.NaN;
  for (const p of [...under, ...now.tiles]) {
    const [l, h] = rangeOf(p.table);
    if (!(l >= lo)) lo = l;
    if (!(h <= hi)) hi = h;
  }
  if (lay.extents.any && !Number.isNaN(lo)) {
    lay.datum.value = lo;
    lay.extents.add(lay.datum);
    lay.datum.value = hi;
    lay.extents.add(lay.datum);
    lay.extents.flush(ctx.theme);
  }
  const tokens = ctx.tokens.version;
  const xLen = region.src.length ?? Number.MAX_SAFE_INTEGER;
  // a tile's last bin may run past its system's end on either axis (a level wider than one position): its box ends at
  // the end, never past it (a position past the axis has no place, and the tile would be dropped)
  const yLen = cov?.systems.find((s) => s.name === region.yName)?.length ?? Number.MAX_SAFE_INTEGER;
  const place = (p: Placed, level: number) => {
    const key = `${lay.key}|r|${p.key}|${p.transposed ? 't' : ''}|${lay.colorVersion}.${tokens}|${lay.orderKey}`;
    let img = lay.cache!.get(key);
    if (!img) {
      img = buildTile(lay, p, level);
      if (!img) return;
      lay.cache!.put(key, img);
    }
    const xs = p.x * level;
    const ys = p.y * level;
    if (boxOf(a, ay, xs, Math.min(xs + p.w * level, xLen), ys, Math.min(ys + p.h * level, yLen), lay.ib))
      lay.images.push(img);
  };
  // nothing held under a pending rect: the pending fill there (EC-169: loading, never empty)
  if (!under.length)
    for (let k = 0; k + 3 < now.pending.length; k += 4)
      boxOf(a, ay, now.pending[k], now.pending[k + 1], now.pending[k + 2], now.pending[k + 3], lay.pend);
  for (const p of under) place(p, underLevel);
  for (const p of now.tiles) place(p, L);
  // what a hover reads (`regionCellAt`): the tiles drawn at the level read
  lay.tiles = now.tiles;
}

/**
 * A region grid's cell value at positions (i, j) from the tiles drawn at the level read: NaN where the tile holds no
 * value (a sparse tile's unlisted cell is its `fill`), undefined where no tile is loaded (pending).
 */
export function regionCellAt(lay: GridLay, i: number, j: number): number | undefined {
  const L = lay.binLevel;
  const bx = Math.floor(i / L);
  const by = Math.floor(j / L);
  for (const p of lay.tiles) {
    if (bx < p.x || by < p.y || bx >= p.x + p.w || by >= p.y + p.h) continue;
    const t = p.table;
    // the served table's own bins: a mirror's swapped
    const sx = p.transposed ? by : bx;
    const sy = p.transposed ? bx : by;
    if (t.form === 'dense') {
      const a = sx - t.x;
      const b = sy - t.y;
      return a < 0 || b < 0 || a >= t.nx || b >= t.ny ? Number.NaN : t.value[a * t.ny + b];
    }
    for (let k = 0; k < t.value.length; k++) if (t.xBin[k] === sx && t.yBin[k] === sy) return t.value[k];
    return t.fill;
  }
  return undefined;
}
