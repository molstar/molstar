/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot (seqspec:format#layer-types `dot_plot`, seqspec:interaction#hover), part of the renderer API: two axes
// over one section's layout coordinate system, x across as every track's and y down the track showing the same window
// scaled to the track's height (`dotPlotY`); on a stacked section, its lane's block both ways.
// Prepare places each grid row's two ends as a pair's (their placements' paths) into one float32 grid per block pair,
// sized to the indices its rows place within the blocks' pieces (the row's own Float32Array where it is the whole grid
// and its paths the identity), and, where both blocks are shown whole in one direction, a mean pyramid in shown order
// (2^L × 2^L cells a bin, NaN left out, down to one bin); what it allocates is at most `DotPlotLimits.cells` floats
// (past it the pair is left out and reported), and the same data placed the same way reuses both (a params patch).
// Pairs are drawn by their first and last cells: a pair of single positions a dot at every zoom (a circle in its cell,
// its rect below a device px), a matched stretch a diagonal segment. A grid is drawn as image tiles of `DOT_PLOT_TILE`
// samples a side over its own extent, at the level whose sample is a device pixel or more, built through the color
// theme on first use and kept per prepared state, theme, axis and backend. A progressive backend builds a budget of
// samples a composer frame (`samplesPerFrame`, every lane's draw sharing it), drawing a coarser level meanwhile (the
// one with its blocks about 32 samples a side, built whatever the budget) and asking for another frame (`v.more`: the
// composer draws the plot's entries again alone); past `DotPlotLimits.tiles` the least recently drawn tiles go, never
// one drawn in the last two frames. The overlay draws every channel's lit runs as rows (the composer draws them as
// columns) and the rectangle being dragged (the track's first dot plot), and the hovered cell's outline (the first
// with a grid there). Allocation-free once warm.

import { endsOf } from "@molstar/seq-core/format";
import { type ComposedPath, Issue, type SeqSet } from "@molstar/seq-core/model";
import { AxisGeometry, type AxisView } from "../axis.ts";
import { MarkBit } from "../marks.ts";
import { partOf } from "../prepare.ts";
import type { LayerInput, PrepareContext } from "../registry.ts";
import { ThemeLocation } from "../themes.ts";
import { inkOf } from "./arcs.ts";
import {
  BackendCaps,
  type BackendImages,
  GlyphBatch,
  GlyphShape,
  NO_COLOR,
  PathBatch,
  PathMode,
  RectBatch,
  SkinToken,
} from "./backend.ts";
import { visibleRange } from "./cells.ts";
import { keyAlpha, packAlpha } from "./density.ts";
import { hoverEmphasis, outlineItem } from "./hover.ts";
import { DOT_R } from "./sizes.ts";
import type { DrawView } from "./view.ts";

/**
 * The bounds: the floats a block pair's grid allocates (its own copy, unless one Float32Array row is the whole grid on
 * the identity, which it references, plus its mean pyramid, about a third of its cells; past it the pair is left out
 * and reported), the samples a progressive backend builds a frame, and the tiles a backend keeps (soft: a tile drawn in
 * the last two frames stays; past it the least recently drawn go).
 */
export const DotPlotLimits = { cells: 2 ** 24, samplesPerFrame: 2 ** 17, tiles: 96 };
/** What frames did: tiles built and the samples they read (the tests' and benches' counters). */
export const DotPlotVisits = { tiles: 0, samples: 0 };

/** A level of a grid's pyramid: bins of 2^L × 2^L shown cells from bin (bx0, by0), w × h, x-major; NaN: none. */
export interface DotPlotLevel {
  readonly bx0: number;
  readonly by0: number;
  readonly w: number;
  readonly h: number;
  readonly mean: Float32Array;
}
/** A block pair's values: block bx's indices [x0, x0 + nx) across by's [y0, y0 + ny), x-major; NaN: none. */
export interface DotPlotGrid {
  /** Its place in `DotPlotState.grids`. */
  readonly index: number;
  readonly bx: number;
  readonly by: number;
  readonly x0: number;
  readonly nx: number;
  readonly y0: number;
  readonly ny: number;
  /** The first annotation whose rows fill it (the theme's location). */
  readonly annotation: number;
  readonly values: Float32Array;
  /**
   * Both blocks shown whole in one direction (shown k is `u`, `first + dir·u` the block index): levels 1.. of the mean
   * pyramid in shown order, down to one bin (a 1 × 1 grid has that one), so a tile's sample is one bin; null otherwise
   * (a sample is its cells, through `indexOf`).
   */
  readonly levels: readonly DotPlotLevel[] | null;
  readonly firstX: number;
  readonly dirX: number;
  readonly firstY: number;
  readonly dirY: number;
}
/**
 * Grids per block pair (`gridAt` finds one), and pairs as segments: per segment its x block and first and last block
 * indices, its y block's, its annotation and row.
 */
export interface DotPlotState {
  readonly set: SeqSet;
  readonly grids: readonly DotPlotGrid[];
  /** Grid by block pair, `bx · blocks + by`. */
  readonly byPair: ReadonlyMap<number, number>;
  readonly blocks: number;
  readonly count: number;
  readonly aBlock: Int32Array;
  readonly a0: Int32Array;
  readonly a1: Int32Array;
  readonly bBlock: Int32Array;
  readonly b0: Int32Array;
  readonly b1: Int32Array;
  readonly annotation: Int32Array;
  readonly row: Int32Array;
  /** The layer is its track's first dot plot (`firstDotPlot`): it draws the overlay, once for the track. */
  readonly first: boolean;
  /** The layer is its track's last dot plot: it outlines the hovered cell where no dot plot before it had its grid. */
  readonly last: boolean;
  /** The renderer's hysteresis (whole px): the height given, the last three targets, the target held at; 0: none. */
  readonly sized: Int32Array;
}

/**
 * The first `dot_plot` layer of layer k's track, among `n` layers in order (their track and declared type by index):
 * the one that draws a track's dot plot overlay and that the drag's rectangle names. Edge rate.
 */
export function firstDotPlot(n: number, track: (j: number) => number, type: (j: number) => string, k: number): number {
  const t = track(k);
  for (let j = 0; j < n; j++) if (track(j) === t && type(j) === "dot_plot") return j;
  return k;
}

const part = new Int32Array(2);
/** Block index of source position i through a placement's path (its larger part's first); -1: none. */
const blockIndex = (path: ComposedPath | null, i: number) => (partOf(path, i, part) ? part[0] : -1);

interface Placed {
  block: number;
  path: ComposedPath | null;
}
/** A grid row on a block pair, and the block indices it places (lo inclusive, hi exclusive). */
interface Part {
  values: Float32Array | Float64Array;
  A: Placed;
  B: Placed;
  s: number;
  t: number;
  u: number;
  w: number;
  annotation: number;
}

/** Headless and pure (seqspec:viewer#extension-points): grids, their pyramids and segments in block indices. */
export function prepareDotPlot(input: LayerInput, rc?: PrepareContext): DotPlotState {
  const { set, layer, section } = input;
  const at = new Map<string, Placed[]>();
  for (const pl of layer.placements) {
    const k = `${pl.annotation}:${pl.slot}`;
    const list = at.get(k) ?? [];
    list.push({ block: pl.block, path: pl.path >= 0 ? layer.paths[pl.path] : null });
    at.set(k, list);
  }
  const nb = section.blocks.length;
  const parts = new Map<number, Part[]>();
  const seg: number[][] = [];
  for (const id of layer.data) {
    const a = set.annotations[id];
    if (!a?.pairs) continue;
    const kept = layer.kept[layer.data.indexOf(id)];
    const bSlot = new Int32Array(a.count).fill(-1);
    for (const e of endsOf(a)) if (e.side === 1) for (const r of e.rows) bSlot[r] = e.slot;
    for (let r = 0; r < a.count; r++) {
      if (kept && !kept[r]) continue;
      const [s, t, u, w] = [a.start[r], a.stop[r], a.pairs.bStart[r], a.pairs.bStop[r]];
      for (const A of at.get(`${id}:${a.slot[r]}`) ?? [])
        for (const B of at.get(`${id}:${bSlot[r]}`) ?? []) {
          // A ligand block is one position for its atoms: no grid or segment across it.
          if (section.blocks[A.block]?.ligand || section.blocks[B.block]?.ligand) continue;
          if (!a.grid) {
            const ends = [blockIndex(A.path, s), blockIndex(A.path, Math.max(s, t - 1))];
            const bEnds = [blockIndex(B.path, u), blockIndex(B.path, Math.max(u, w - 1))];
            if (ends[0] >= 0 && ends[1] >= 0 && bEnds[0] >= 0 && bEnds[1] >= 0)
              seg.push([A.block, ends[0], ends[1], B.block, bEnds[0], bEnds[1], id, r]);
            continue;
          }
          const key = A.block * nb + B.block;
          const list = parts.get(key) ?? [];
          list.push({ values: a.grid.rows[r], A, B, s, t, u, w, annotation: id });
          parts.set(key, list);
        }
    }
  }
  const grids: DotPlotGrid[] = [];
  const byPair = new Map<number, number>();
  for (const [key, list] of parts) {
    const g = gridOf(input, list, grids.length, rc);
    if (!g) continue;
    byPair.set(key, grids.length);
    grids.push(g);
  }
  const col = (c: number) => Int32Array.from(seg, (x) => x[c]);
  const ls = section.layers;
  const own = ls.findIndex((l) => l.ref === layer.ref);
  const first =
    own >= 0 &&
    firstDotPlot(
      ls.length,
      (j) => ls[j].track,
      (j) => ls[j].type,
      own,
    ) === own;
  let last = own >= 0;
  for (let j = own + 1; last && j < ls.length; j++)
    if (ls[j].track === ls[own].track && ls[j].type === "dot_plot") last = false;
  return {
    set,
    grids,
    byPair,
    blocks: nb,
    count: seg.length,
    aBlock: col(0),
    a0: col(1),
    a1: col(2),
    bBlock: col(3),
    b0: col(4),
    b1: col(5),
    annotation: col(6),
    row: col(7),
    first,
    last,
    sized: new Int32Array(5),
  };
}

/** A block's piece in its own indices, [lo, hi): its `range`, else (a whole unit, a cyclic one rotated) its domain. */
function pieceOf(input: LayerInput, b: number, out: Int32Array) {
  const r = input.section.blocks[b];
  const ranged = r.stop > r.start;
  out[0] = ranged ? r.start : r.base;
  out[1] = ranged ? r.stop : r.base + r.length;
}
const pieceX = new Int32Array(2);
const pieceY = new Int32Array(2);
/**
 * One block pair's grid from its rows (prepare-time): sized to the block indices they place within the pieces, the
 * single row's own Float32Array where it is the whole grid on the identity, its pyramid; null (reported) past the cap.
 */
function gridOf(input: LayerInput, list: Part[], index: number, rc?: PrepareContext): DotPlotGrid | null {
  const { A: A0, B: B0 } = list[0];
  const [bx, by] = [A0.block, B0.block];
  pieceOf(input, bx, pieceX);
  pieceOf(input, by, pieceY);
  let [x0, x1, y0, y1] = [Number.MAX_SAFE_INTEGER, -1, Number.MAX_SAFE_INTEGER, -1];
  /** The block indices rows [s, t) place within the piece, into `out`'s hull: on the identity the clip, else a walk. */
  const range = (p: ComposedPath | null, s: number, t: number, piece: Int32Array, out: number[]) => {
    if (!p) {
      const [lo, hi] = [Math.max(s, piece[0]), Math.min(t, piece[1])];
      if (hi > lo) [out[0], out[1]] = [Math.min(out[0], lo), Math.max(out[1], hi)];
      return;
    }
    for (let i = s; i < t; i++) {
      const j = blockIndex(p, i);
      if (j < piece[0] || j >= piece[1]) continue;
      if (j < out[0]) out[0] = j;
      if (j + 1 > out[1]) out[1] = j + 1;
    }
  };
  const ex = [x0, x1];
  const ey = [y0, y1];
  for (const q of list) {
    range(q.A.path, q.s, q.t, pieceX, ex);
    range(q.B.path, q.u, q.w, pieceY, ey);
  }
  [x0, x1, y0, y1] = [ex[0], ex[1], ey[0], ey[1]];
  if (x1 <= x0 || y1 <= y0) return null;
  const [nx, ny] = [x1 - x0, y1 - y0];
  const label = (b: number) => {
    const row = input.section.blocks[b];
    const s = input.set;
    return row.label >= 0 ? s.strings.get(row.label) : s.strings.get(s.units.key[row.unit] ?? -1);
  };
  const one = list.length === 1 ? list[0] : null;
  const whole =
    one !== null &&
    !one.A.path &&
    !one.B.path &&
    one.values instanceof Float32Array &&
    one.s === x0 &&
    one.t === x1 &&
    one.u === y0 &&
    one.w === y1;
  const sx = input.section.blocks[bx];
  const sy = input.section.blocks[by];
  const simple = (r: typeof sx) => r.positions === null && r.layout.count === 1;
  const g = {
    index,
    bx,
    by,
    x0,
    nx,
    y0,
    ny,
    annotation: list[0].annotation,
    values: whole ? (one.values as Float32Array) : EMPTY,
    levels: null as DotPlotLevel[] | null,
    firstX: sx.first,
    dirX: sx.dir < 0 ? -1 : 1,
    firstY: sy.first,
    dirY: sy.dir < 0 ? -1 : 1,
  };
  const pyramided = simple(sx) && simple(sy);
  // What it allocates: its own copy unless it references the row, and its pyramid.
  const own = whole ? 0 : nx * ny;
  const floats = own + (pyramided ? pyramidSize(g) : 0);
  if (floats > DotPlotLimits.cells) {
    const what = whole ? "its mean pyramid" : "a copy of its cells and its mean pyramid";
    const message =
      `"${input.layer.ref}": the grid on ${label(bx)} × ${label(by)} is ${nx} × ${ny} = ${nx * ny} cells, and ${what} ` +
      `would take ${floats} floats, past the ${DotPlotLimits.cells} a block pair may allocate: it is left out (narrow ` +
      `the blocks' range, or give the block pair one Float32Array row covering exactly the cells it places, which ` +
      `is referenced, not copied)`;
    rc?.report?.(Issue.warning("svs.render.omitted", "", message));
    return null;
  }
  // The same data placed the same way (a layer's color or other params patched) reuses its cells and pyramid.
  const key = `${list.map(partKey).join(";")}|${x0},${nx},${y0},${ny},${g.firstX},${g.dirX},${g.firstY},${g.dirY},${pyramided}`;
  const byData = built.get(list[0].values) ?? new Map<string, Built>();
  built.set(list[0].values, byData);
  const was = byData.get(key);
  if (was) {
    g.values = was.values;
    g.levels = was.levels;
    return g;
  }
  if (!whole) {
    const values = new Float32Array(nx * ny).fill(Number.NaN);
    for (const q of list) {
      const bw = q.w - q.u;
      const ys = Int32Array.from({ length: bw }, (_, j) => (q.B.path ? blockIndex(q.B.path, q.u + j) : q.u + j) - y0);
      for (let i = q.s; i < q.t; i++) {
        const x = (q.A.path ? blockIndex(q.A.path, i) : i) - x0;
        if (x < 0 || x >= nx) continue;
        const o = (i - q.s) * bw;
        for (let j = 0; j < bw; j++) if (ys[j] >= 0 && ys[j] < ny) values[x * ny + ys[j]] = q.values[o + j];
      }
    }
    g.values = values;
  }
  if (pyramided) g.levels = pyramid(g);
  // A few placements per data (a range or layout change makes a new one): the oldest go.
  if (byData.size >= 8) byData.clear();
  byData.set(key, { values: g.values, levels: g.levels });
  return g;
}
const EMPTY = new Float32Array(0);
interface Built {
  values: Float32Array;
  levels: DotPlotLevel[] | null;
}
/** Per data (a grid row's values), its grids by how they were placed: re-prepares at new params reuse them. */
const built = new WeakMap<object, Map<string, Built>>();
/** An object's id for the placement key (rows and paths by identity); 0: none. */
const ids = new WeakMap<object, number>();
let lastId = 0;
function idOf(o: object | null): number {
  if (!o) return 0;
  const id = ids.get(o);
  if (id !== undefined) return id;
  ids.set(o, ++lastId);
  return lastId;
}
const partKey = (q: Part) => `${idOf(q.values)}:${q.s}:${q.t}:${q.u}:${q.w}:${idOf(q.A.path)}:${idOf(q.B.path)}`;

/** `u` of block index i (shown order of a block shown whole). */
const uOf = (first: number, dir: number, i: number) => (dir > 0 ? i - first : first - i);
type Placing = Omit<DotPlotGrid, "levels" | "values">;
/** The grid's hull in shown order, [ux0, ux1] × [uy0, uy1] (inclusive), into `hull`. */
const hull = new Int32Array(4);
function hullOf(g: Placing) {
  const a = uOf(g.firstX, g.dirX, g.x0);
  const b = uOf(g.firstX, g.dirX, g.x0 + g.nx - 1);
  const c = uOf(g.firstY, g.dirY, g.y0);
  const d = uOf(g.firstY, g.dirY, g.y0 + g.ny - 1);
  hull[0] = Math.min(a, b);
  hull[1] = Math.max(a, b);
  hull[2] = Math.min(c, d);
  hull[3] = Math.max(c, d);
}
/** The bins of levels 1.. down to one bin (what `pyramid` allocates). */
function pyramidSize(g: Placing): number {
  hullOf(g);
  let n = 0;
  for (let L = 1; L < LEVELS; L++) {
    const w = (hull[1] >> L) - (hull[0] >> L) + 1;
    const h = (hull[3] >> L) - (hull[2] >> L) + 1;
    n += w * h;
    if (w === 1 && h === 1) break;
  }
  return n;
}
/**
 * Levels 1.. of a grid's mean pyramid in shown order, down to one bin: bin j of level L holds u ∈ [j·2^L, (j + 1)·2^L).
 * Each level from the one below, weighted by its cells with a value (counts kept only while building).
 */
function pyramid(g: Omit<DotPlotGrid, "levels">): DotPlotLevel[] {
  hullOf(g);
  const [ux0, ux1, uy0, uy1] = [hull[0], hull[1], hull[2], hull[3]];
  const levels: DotPlotLevel[] = [];
  let prevCount: Float64Array | null = null;
  for (let L = 1; L < LEVELS; L++) {
    const [bx0, by0] = [ux0 >> L, uy0 >> L];
    const [w, h] = [(ux1 >> L) - bx0 + 1, (uy1 >> L) - by0 + 1];
    const sum = new Float64Array(w * h);
    const count = new Float64Array(w * h);
    const prev = levels[levels.length - 1];
    if (!prev || !prevCount) {
      for (let x = 0; x < g.nx; x++) {
        const bi = (uOf(g.firstX, g.dirX, g.x0 + x) >> 1) - bx0;
        for (let y = 0; y < g.ny; y++) {
          const val = g.values[x * g.ny + y];
          if (Number.isNaN(val)) continue;
          const o = bi * h + (uOf(g.firstY, g.dirY, g.y0 + y) >> 1) - by0;
          sum[o] += val;
          count[o]++;
        }
      }
    } else
      for (let x = 0; x < prev.w; x++)
        for (let y = 0; y < prev.h; y++) {
          const c = prevCount[x * prev.h + y];
          if (c === 0) continue;
          const o = (((prev.bx0 + x) >> 1) - bx0) * h + ((prev.by0 + y) >> 1) - by0;
          sum[o] += prev.mean[x * prev.h + y] * c;
          count[o] += c;
        }
    const mean = new Float32Array(w * h);
    for (let k = 0; k < mean.length; k++) mean[k] = count[k] > 0 ? sum[k] / count[k] : Number.NaN;
    levels.push({ bx0, by0, w, h, mean });
    prevCount = count;
    if (w === 1 && h === 1) break;
  }
  return levels;
}

/** The grid's sample at level L over shown [kx0, kx1) × [ky0, ky1) into `mean[0]` (NaN: no value), its bin's mean. */
const mean = new Float64Array(1);
/** A sample's y cells' grid rows (the per-cell path), grown only. */
let rowsOf = new Int32Array(256);
function sampleAt(ax: NonNullable<AxisView["axis"]>, g: DotPlotGrid, L: number, sx: number, sy: number) {
  const levels = g.levels;
  if (levels !== null) {
    if (L === 0) {
      const x = (g.dirX > 0 ? g.firstX + sx : g.firstX - sx) - g.x0;
      const y = (g.dirY > 0 ? g.firstY + sy : g.firstY - sy) - g.y0;
      mean[0] = x >= 0 && x < g.nx && y >= 0 && y < g.ny ? g.values[x * g.ny + y] : Number.NaN;
      return;
    }
    if (L <= levels.length) {
      const lv = levels[L - 1];
      const x = sx - lv.bx0;
      const y = sy - lv.by0;
      mean[0] = x >= 0 && x < lv.w && y >= 0 && y < lv.h ? lv.mean[x * lv.h + y] : Number.NaN;
      return;
    }
    // Coarser than the pyramid's last level, its one bin: the sample holding it.
    const lv = levels[levels.length - 1];
    const shift = L - levels.length;
    mean[0] = sx === lv.bx0 >> shift && sy === lv.by0 >> shift ? lv.mean[0] : Number.NaN;
    return;
  }
  // Blocks not shown whole in one direction: the sample's cells, through the layout map.
  const step = 2 ** L;
  const ky0 = sy * step;
  const ky1 = Math.min(ax.shown[g.by], ky0 + step);
  if (rowsOf.length < step) rowsOf = new Int32Array(step);
  for (let k = ky0; k < ky1; k++) {
    const j = ax.indexOf(g.by, k as never) - g.y0;
    rowsOf[k - ky0] = j >= 0 && j < g.ny ? j : -1;
  }
  let sum = 0;
  let n = 0;
  const kx1 = Math.min(ax.shown[g.bx], (sx + 1) * step);
  for (let k = sx * step; k < kx1; k++) {
    const i = ax.indexOf(g.bx, k as never) - g.x0;
    if (i < 0 || i >= g.nx) continue;
    for (let r = 0; r < ky1 - ky0; r++) {
      const val = rowsOf[r] >= 0 ? g.values[i * g.ny + rowsOf[r]] : Number.NaN;
      if (Number.isNaN(val)) continue;
      sum += val;
      n++;
    }
  }
  mean[0] = n > 0 ? sum / n : Number.NaN;
}

/** A sample's mean at level L (tests; `sampleAt` is the draw's, through a module slot): NaN where none has a value. */
export function dotPlotSample(ax: NonNullable<AxisView["axis"]>, g: DotPlotGrid, L: number, sx: number, sy: number) {
  sampleAt(ax, g, L, sx, sy);
  return mean[0];
}

// ---- The y axis --------------------------------------------------------------------------------------------------------

/**
 * The device y of block b's shown k's top edge: its x (device-snapped, `xd`) across the data area mapped down the
 * entry, so y shows the section's window scaled to the track's height. Ints in and out.
 */
export function dotPlotY(v: DrawView, b: number, k: number): number {
  const w = v.right - v.left;
  return w > 0 ? Math.round((v.xd(b, k) - v.left * v.dpr) * (v.height / w) + v.top * v.dpr) | 0 : 0;
}

/**
 * The block and shown position at content x `io[0]` (CSS px, as `blockX`) into `out` (block, shown): on a stacked
 * section `only`'s (its lane's block), else the block whose start is at or left of it (the search hit.ts does inline);
 * outside a block, with `clamp`, the nearer block's end, else false.
 */
export function shownAtX(v: AxisView, io: Float64Array, clamp: boolean, only: number, out: Int32Array): boolean {
  const ax = v.axis;
  if (!ax || ax.count === 0) return false;
  const x = io[0];
  let b = ax.stacked ? Math.max(0, only) : 0;
  for (let hi = ax.stacked ? b : ax.count - 1; b < hi; ) {
    const mid = (b + hi + 1) >> 1;
    if (v.blockX[mid] <= x) b = mid;
    else hi = mid - 1;
  }
  let k = Math.floor((x - v.blockX[b]) / (ax.span[b] * v.cell));
  const n = ax.shown[b];
  if (k < 0 || k >= n) {
    if (!clamp || n === 0) return false;
    // In the gap after block b: the nearer of its last and the next block's first.
    const next =
      !ax.stacked && b + 1 < ax.count && k >= n
        ? v.blockX[b + 1] - x < x - (v.blockX[b] + n * ax.span[b] * v.cell)
        : false;
    if (next) {
      b++;
      k = 0;
    } else k = k < 0 ? 0 : n - 1;
  }
  out[0] = b;
  out[1] = k;
  return true;
}

/**
 * The y axis's block and shown position at surface y `io[1]` of an entry at [top, top + height) (CSS px): its content
 * x, `left + (y − top) · width / height`, into `io[0]`, then `shownAtX`.
 */
export function dotPlotAt(
  v: AxisView,
  top: number,
  height: number,
  io: Float64Array,
  clamp: boolean,
  only: number,
  out: Int32Array,
): boolean {
  const y = clamp ? Math.max(top, Math.min(top + height - 0.01, io[1])) : io[1];
  if (height <= 0 || y < top || y >= top + height) return false;
  io[0] = v.left + ((y - top) * (v.right - v.left)) / height;
  return shownAtX(v, io, clamp, only, out);
}

// ---- Tiles -----------------------------------------------------------------------------------------------------------

/** Samples across a tile's side. */
export const DOT_PLOT_TILE = 256;
const LEVELS = 32;
/**
 * A grid's tiles at one level: over the grid's own samples (tiles [tx0, tx0 + nx) × [ty0, ty0 + ny)), per tile its
 * image id (-1 until built, -2 empty) and the frame it was last drawn in.
 */
export interface TileLevel {
  readonly tx0: number;
  readonly ty0: number;
  readonly nx: number;
  readonly ny: number;
  readonly ids: Int32Array;
  readonly used: Int32Array;
}
/** One backend's tiles: per grid and level (`index · LEVELS + L`). */
export interface DotPlotTileSet {
  images: BackendImages;
  levels: (TileLevel | undefined)[];
  /** Images held. */
  count: number;
}
/**
 * Per prepared state (the engine's `cache`): the theme and axis the tiles were made under, a set per backend, colors,
 * each grid's shown extent, and the frame count with its sample budget.
 */
export class DotPlotTiles {
  state: DotPlotState | null = null;
  theme: unknown = null;
  axis: unknown = null;
  readonly sets: DotPlotTileSet[] = [];
  /** Per segment its color key (alpha folded in, as density's), made under `theme`. */
  keys: Int32Array | null = null;
  /** Per grid its shown extent under `axis`: [kx0, kx1) × [ky0, ky1), 4 ints a grid; null until drawn. */
  extent: Int32Array | null = null;
  /** The composer frame last drawn in (`DrawView.frameId`), and this cache's own count of frames. */
  at = -1;
  frame = 0;
  /** Samples a progressive backend may still build this frame (every lane and grid of the layer share it). */
  budget = 0;
  dispose(): void {
    for (const t of this.sets) for (let k = 0; k < t.levels.length; k++) drop(t, k);
    this.sets.length = 0;
    this.keys = this.extent = null;
  }
}
function drop(t: DotPlotTileSet, k: number) {
  const lv = t.levels[k];
  for (let i = 0; lv && i < lv.ids.length; i++)
    if (lv.ids[i] >= 0) {
      t.images.dispose(lv.ids[i]);
      t.count--;
    }
  t.levels[k] = undefined;
}
/**
 * Where the engine gives no cache (an SVG export, a test's frame): one per prepared state, its images the backend's,
 * dropped with it (an export's backend lives for one export).
 */
const own = new WeakMap<DotPlotState, DotPlotTiles>();
/** Edge-rate, apart from the draw. */
function newTiles(s: DotPlotState): DotPlotTiles {
  const c = new DotPlotTiles();
  own.set(s, c);
  return c;
}
function newSet(c: DotPlotTiles, images: BackendImages): DotPlotTileSet {
  // Packed: a slot per grid and level, none made yet.
  const levels = new Array<TileLevel | undefined>(LEVELS * Math.max(1, c.state?.grids.length ?? 1)).fill(undefined);
  const t = { images, levels, count: 0 };
  c.sets.push(t);
  return t;
}
function tilesOf(v: DrawView, s: DotPlotState): DotPlotTiles {
  const c = v.cache instanceof DotPlotTiles ? v.cache : (own.get(s) ?? newTiles(s));
  if (c.state !== s || c.theme !== v.color || c.axis !== v.axis) {
    c.dispose();
    c.state = s;
    c.theme = v.color;
    c.axis = v.axis;
  }
  if (!c.extent) c.extent = extentsOf(v, s);
  return c;
}
function tileSet(c: DotPlotTiles, images: BackendImages): DotPlotTileSet {
  for (let i = 0; i < c.sets.length; i++) if (c.sets[i].images === images) return c.sets[i];
  return newSet(c, images);
}
/** Each grid's shown extent (edge rate, once per axis): its indices' shown positions' hull per block; empty: none. */
function extentsOf(v: DrawView, s: DotPlotState): Int32Array {
  const out = new Int32Array(4 * s.grids.length);
  const side = (b: number, i0: number, n: number, o: number) => {
    let [lo, hi] = [0x7fffffff, -1];
    for (let i = i0; i < i0 + n; i++) {
      const k = v.shownOf(b, i);
      if (k < 0) continue;
      if (k < lo) lo = k;
      if (k > hi) hi = k;
    }
    out[o] = hi >= lo ? lo : 0;
    out[o + 1] = hi >= lo ? hi + 1 : 0;
  };
  for (const g of s.grids) {
    side(g.bx, g.x0, g.nx, 4 * g.index);
    side(g.by, g.y0, g.ny, 4 * g.index + 2);
  }
  return out;
}

/** The quiet gray drawn without a theme (a color key; as a pixel, `pixelOf`). */
const QUIET = 0x8a8f98;
/** A cell's color key at its value under the theme (alpha folded in, as density's); `NO_COLOR` where it gives none. */
function cellKey(v: DrawView, g: DotPlotGrid, k: number, value: number): number {
  const theme = v.color;
  const ax = v.axis;
  if (!theme || !ax) return QUIET;
  // A clean location: another draw's row or item would make the theme read that row's value, not the cell's.
  ThemeLocation.reset(v.loc);
  const loc = ThemeLocation.at(v.loc, ax, v.layer, g.bx, ax.indexOf(g.bx, k as never));
  loc.annotation = g.annotation;
  loc.value = value;
  const c = theme.value(loc);
  const a = c >= 0 && theme.alpha ? theme.alpha(loc) : 255;
  ThemeLocation.reset(loc);
  return c < 0 ? NO_COLOR : packAlpha(c, a);
}
/** A color key as a pixel's bits (ImageData's byte order, R first in memory); 0: none. */
function pixelOf(key: number): number {
  if (key === NO_COLOR || key < 0) return 0;
  const a = 255 - ((key >>> 24) & 255);
  return ((a << 24) | ((key & 255) << 16) | (key & 0xff00) | ((key >> 16) & 255)) >>> 0;
}
const px = new Uint32Array(DOT_PLOT_TILE * DOT_PLOT_TILE);
/** The samples the last `build` read (its cost against the budget). */
const built1 = new Int32Array(1);
/**
 * Grid g's tile (tx, ty) at level L (edge rate): per sample in the grid's extent [s0, s1) × [t0, t1) its bin's mean
 * through the color theme; -2: empty.
 */
function build(v: DrawView, images: BackendImages, g: DotPlotGrid, L: number, tx: number, ty: number, e: Int32Array) {
  const ax = v.axis;
  built1[0] = 0;
  if (!ax) return -2;
  const step = 2 ** L;
  const T = DOT_PLOT_TILE;
  const W = Math.min(T, Math.ceil(ax.shown[g.bx] / step) - tx * T);
  const H = Math.min(T, Math.ceil(ax.shown[g.by] / step) - ty * T);
  if (W <= 0 || H <= 0) return -2;
  const x0 = Math.max(0, e[0] - tx * T);
  const x1 = Math.min(W, e[1] - tx * T);
  const y0 = Math.max(0, e[2] - ty * T);
  const y1 = Math.min(H, e[3] - ty * T);
  if (x1 <= x0 || y1 <= y0) return -2;
  DotPlotVisits.tiles++;
  built1[0] = (x1 - x0) * (y1 - y0);
  DotPlotVisits.samples += built1[0];
  px.fill(0, 0, W * H);
  let any = false;
  for (let x = x0; x < x1; x++)
    for (let y = y0; y < y1; y++) {
      sampleAt(ax, g, L, tx * T + x, ty * T + y);
      if (Number.isNaN(mean[0])) continue;
      const p = pixelOf(cellKey(v, g, (tx * T + x) * step, mean[0]));
      px[y * W + x] = p;
      if (p !== 0) any = true;
    }
  return any ? images.fromPixels(W, H, px) : -2;
}

/** The level whose sample (2^L shown cells a side) is at least a device pixel both ways. */
function levelOf(v: DrawView, g: DotPlotGrid): number {
  const w = v.right - v.left;
  const cx = v.cell * v.span[g.bx] * v.dpr;
  const cy = (v.cell * v.span[g.by] * v.dpr * v.height) / (w > 0 ? w : 1);
  const m = cx < cy ? cx : cy;
  return m >= 1 || m <= 0 ? 0 : Math.min(LEVELS - 1, Math.ceil(Math.log2(1 / m)));
}
/** Samples a side of the level a draw always has (built past the budget): its blocks about 16 to 32 samples a side. */
const COARSE = 32;
/** The coarsest level a draw falls back to, built whatever the budget: the larger block `COARSE` samples or fewer. */
function coarseOf(ax: NonNullable<AxisView["axis"]>, g: DotPlotGrid): number {
  let L = 0;
  while (L < LEVELS - 1 && Math.max(ax.shown[g.bx], ax.shown[g.by]) > COARSE * 2 ** L) L++;
  return L;
}

const r = new Float64Array(8);
const xs = new Int32Array(2);
const ys = new Int32Array(2);
/** Level L's tiles of grid g over its own samples (edge rate: made on first use). */
function levelTiles(c: DotPlotTiles, set: DotPlotTileSet, g: DotPlotGrid, L: number): TileLevel {
  const key = g.index * LEVELS + L;
  const lv = set.levels[key];
  if (lv) return lv;
  const e = c.extent as Int32Array;
  const o = 4 * g.index;
  const step = 2 ** L;
  const T = DOT_PLOT_TILE;
  const tx0 = Math.floor(Math.floor(e[o] / step) / T);
  const ty0 = Math.floor(Math.floor(e[o + 2] / step) / T);
  const nx = e[o + 1] > e[o] ? Math.floor((Math.ceil(e[o + 1] / step) - 1) / T) + 1 - tx0 : 0;
  const ny = e[o + 3] > e[o + 2] ? Math.floor((Math.ceil(e[o + 3] / step) - 1) / T) + 1 - ty0 : 0;
  const made = {
    tx0,
    ty0,
    nx,
    ny,
    ids: new Int32Array(nx * ny).fill(-1),
    used: new Int32Array(nx * ny),
  };
  set.levels[key] = made;
  return made;
}
/** Grid g's extent at level L in samples, [s0, s1) × [t0, t1), into `ext`. */
const ext = new Int32Array(4);
function extentAt(c: DotPlotTiles, g: DotPlotGrid, L: number) {
  const e = c.extent as Int32Array;
  const o = 4 * g.index;
  const step = 2 ** L;
  ext[0] = Math.floor(e[o] / step);
  ext[1] = Math.ceil(e[o + 1] / step);
  ext[2] = Math.floor(e[o + 2] / step);
  ext[3] = Math.ceil(e[o + 3] / step);
}
/**
 * Tile (tx, ty) of level L, built if the frame's budget allows (or `force`); its id, -2 empty (or outside the grid),
 * -1 not built (the frame asks for another).
 */
function tileAt(
  v: DrawView,
  c: DotPlotTiles,
  set: DotPlotTileSet,
  g: DotPlotGrid,
  L: number,
  tx: number,
  ty: number,
  force: boolean,
) {
  const lv = levelTiles(c, set, g, L);
  const x = tx - lv.tx0;
  const y = ty - lv.ty0;
  if (x < 0 || x >= lv.nx || y < 0 || y >= lv.ny) return -2;
  const t = y * lv.nx + x;
  lv.used[t] = c.frame;
  if (lv.ids[t] !== -1) return lv.ids[t];
  if (!force && c.budget <= 0) return -1;
  extentAt(c, g, L);
  lv.ids[t] = build(v, set.images, g, L, tx, ty, ext);
  c.budget -= built1[0];
  if (lv.ids[t] >= 0) set.count++;
  return lv.ids[t];
}
/**
 * Level L's samples [a0, a1) × [b0, b1) of tile (tx, ty) drawn over their shown cells; a partial last sample takes its
 * fraction of a source pixel, so every cell lands on its own pixels.
 */
function drawPart(
  v: DrawView,
  g: DotPlotGrid,
  id: number,
  L: number,
  tx: number,
  ty: number,
  a0: number,
  a1: number,
  b0: number,
  b1: number,
) {
  const ax = v.axis as NonNullable<AxisView["axis"]>;
  const step = 2 ** L;
  const T = DOT_PLOT_TILE;
  const kx1 = Math.min(ax.shown[g.bx], a1 * step);
  const ky1 = Math.min(ax.shown[g.by], b1 * step);
  const x0 = v.xd(g.bx, a0 * step);
  const y0 = dotPlotY(v, g.by, b0 * step);
  r[0] = a0 - tx * T;
  r[1] = b0 - ty * T;
  r[2] = (kx1 - a0 * step) / step;
  r[3] = (ky1 - b0 * step) / step;
  r[4] = x0;
  r[5] = y0;
  r[6] = v.xd(g.bx, kx1) - x0;
  r[7] = dotPlotY(v, g.by, ky1) - y0;
  if (r[2] > 0 && r[3] > 0 && r[6] > 0 && r[7] > 0) v.backend.image(id, r, false);
}
/**
 * Level L's tiles over its samples [a0, a1) × [b0, b1): drawn when every one is built or buildable (`force`, or the
 * budget); false (nothing drawn) otherwise.
 */
function cover(
  v: DrawView,
  c: DotPlotTiles,
  set: DotPlotTileSet,
  g: DotPlotGrid,
  L: number,
  a0: number,
  a1: number,
  b0: number,
  b1: number,
  force: boolean,
): boolean {
  const T = DOT_PLOT_TILE;
  const tx0 = Math.floor(a0 / T);
  const tx1 = Math.ceil(a1 / T);
  const ty0 = Math.floor(b0 / T);
  const ty1 = Math.ceil(b1 / T);
  for (let ty = ty0; ty < ty1; ty++)
    for (let tx = tx0; tx < tx1; tx++) if (tileAt(v, c, set, g, L, tx, ty, force) === -1) return false;
  for (let ty = ty0; ty < ty1; ty++)
    for (let tx = tx0; tx < tx1; tx++) {
      const id = tileAt(v, c, set, g, L, tx, ty, false);
      if (id >= 0)
        drawPart(
          v,
          g,
          id,
          L,
          tx,
          ty,
          Math.max(a0, tx * T),
          Math.min(a1, (tx + 1) * T),
          Math.max(b0, ty * T),
          Math.min(b1, (ty + 1) * T),
        );
    }
  return true;
}
/**
 * Grid g over the window at its level, within its extent: each tile in view, built within the budget; one not built
 * yet is drawn from the finest coarser level built (the coarse one, `coarseOf`, is always built), and the frame asks
 * for another.
 */
function drawGrid(v: DrawView, c: DotPlotTiles, set: DotPlotTileSet, g: DotPlotGrid) {
  visibleRange(v, g.bx, xs);
  visibleRange(v, g.by, ys);
  if (xs[0] >= xs[1] || ys[0] >= ys[1]) return;
  const L = levelOf(v, g);
  const top = coarseOf(v.axis as NonNullable<AxisView["axis"]>, g);
  const T = DOT_PLOT_TILE;
  const step = 2 ** L;
  extentAt(c, g, L);
  const u0 = Math.max(ext[0], Math.floor(xs[0] / step));
  const u1 = Math.min(ext[1], Math.ceil(xs[1] / step));
  const w0 = Math.max(ext[2], Math.floor(ys[0] / step));
  const w1 = Math.min(ext[3], Math.ceil(ys[1] / step));
  for (let ty = Math.floor(w0 / T); ty * T < w1; ty++)
    for (let tx = Math.floor(u0 / T); tx * T < u1; tx++) {
      const a0 = Math.max(u0, tx * T);
      const a1 = Math.min(u1, (tx + 1) * T);
      const b0 = Math.max(w0, ty * T);
      const b1 = Math.min(w1, (ty + 1) * T);
      if (a1 <= a0 || b1 <= b0) continue;
      const id = tileAt(v, c, set, g, L, tx, ty, L >= top);
      if (id >= 0) drawPart(v, g, id, L, tx, ty, a0, a1, b0, b1);
      if (id !== -1) continue;
      v.more = true;
      // Meanwhile the finest coarser level built over the same shown cells (the coarse level is built on demand).
      for (let L2 = L + 1; L2 <= top; L2++) {
        const k = 2 ** (L2 - L);
        if (
          cover(
            v,
            c,
            set,
            g,
            L2,
            Math.floor(a0 / k),
            Math.ceil(a1 / k),
            Math.floor(b0 / k),
            Math.ceil(b1 / k),
            L2 >= top,
          )
        )
          break;
      }
    }
}
/**
 * Past `DotPlotLimits.tiles`, the least recently drawn tiles go, one frame's at a time; never one drawn in the last two
 * frames (the cap is soft while every tile is in use). Edge rate.
 */
function evict(c: DotPlotTiles, set: DotPlotTileSet) {
  const levels = set.levels;
  while (set.count > DotPlotLimits.tiles) {
    let oldest = 0x7fffffff;
    for (let k = 0; k < levels.length; k++) {
      const lv = levels[k];
      for (let t = 0; lv && t < lv.ids.length; t++)
        if (lv.ids[t] >= 0 && lv.used[t] <= c.frame - 2 && lv.used[t] < oldest) oldest = lv.used[t];
    }
    if (oldest === 0x7fffffff) return;
    for (let k = 0; k < levels.length; k++) {
      const lv = levels[k];
      for (let t = 0; lv && t < lv.ids.length && set.count > DotPlotLimits.tiles; t++)
        if (lv.ids[t] >= 0 && lv.used[t] === oldest) {
          set.images.dispose(lv.ids[t]);
          lv.ids[t] = -1;
          set.count--;
        }
    }
  }
}

// ---- Segments --------------------------------------------------------------------------------------------------------

const paths = new PathBatch(64);
const dots = new RectBatch(64);
const circles = new GlyphBatch(64);
const CIRCLE = GlyphShape.circle;
/** Per segment its color key under the theme (edge rate). */
function keysOf(v: DrawView, s: DotPlotState, c: DotPlotTiles): Int32Array {
  if (c.keys) return c.keys;
  const keys = new Int32Array(s.count);
  const t = v.color;
  const ax = v.axis;
  for (let k = 0; k < s.count; k++) {
    ThemeLocation.reset(v.loc);
    if (ax) ThemeLocation.at(v.loc, ax, v.layer, s.aBlock[k], s.a0[k]);
    v.loc.annotation = s.annotation[k];
    v.loc.row = s.row[k];
    const col = t ? t.value(v.loc) : NO_COLOR;
    const a = col >= 0 && t?.alpha ? t.alpha(v.loc) : 255;
    keys[k] = packAlpha(col, a);
  }
  ThemeLocation.reset(v.loc);
  c.keys = keys;
  return keys;
}
/**
 * Segment k's device box into `seg`: its first cells' near corner, its last's far corner (a reversed end: mirrored).
 * Doubles, computed inline as `xd` and `dotPlotY` do but never cut to int32: a position 10^8 off the view at letter
 * zoom is past 2^31 device px, which `| 0` wrapped back on screen (VL-12).
 */
const seg = new Float64Array(4);
const segRuns = new Int32Array(16);
/**
 * Block b's shown index nearest i among [i, stop] taken in direction `dir` (i itself when shown); -1: none shown. A
 * shown run is a contiguous index run, so the nearest is a run's end.
 */
function nextShown(v: DrawView, b: number, i: number, stop: number, dir: number): number {
  if (v.shownOf(b, i) >= 0) return i;
  const n = Math.min(segRuns.length >> 1, v.runToShown(b, dir > 0 ? i : stop, (dir > 0 ? stop : i) + 1, segRuns));
  let best = -1;
  for (let r = 0; r < n; r++)
    for (let e = 0; e < 2; e++) {
      const j = v.indexOf(b, (e === 0 ? segRuns[2 * r] : segRuns[2 * r + 1] - 1) as never);
      if (j >= 0 && dir * (j - i) >= 0 && (best < 0 || dir * (j - best) < 0)) best = j;
    }
  return best;
}
/** The segment's ends moved inward to the first and last of its steps shown on both axes (`ends[0]`, `ends[1]`). */
const ends = new Int32Array(2);
function shownSteps(v: DrawView, s: DotPlotState, k: number): boolean {
  const ba = s.aBlock[k];
  const bb = s.bBlock[k];
  const a0 = s.a0[k];
  const b0 = s.b0[k];
  const na = s.a1[k] - a0;
  const nb = s.b1[k] - b0;
  const n = Math.max(Math.abs(na), Math.abs(nb));
  const da = na < 0 ? -1 : 1;
  const db = nb < 0 ? -1 : 1;
  for (let side = 0; side < 2; side++) {
    // From the start forward, then from the end back; a few rounds settle a hidden stretch on either axis.
    const dir = side === 0 ? 1 : -1;
    let t = side === 0 ? 0 : n;
    for (let round = 0; round < 4 && t >= 0 && t <= n; round++) {
      const t0 = t;
      const a = a0 + Math.round((t * na) / (n || 1));
      const ea = nextShown(v, ba, a, side === 0 ? s.a1[k] : a0, dir * da);
      if (ea < 0) return false;
      t = Math.abs(na) > 0 ? Math.round((Math.abs(ea - a0) * n) / Math.abs(na)) : t;
      const b = b0 + Math.round((t * nb) / (n || 1));
      const eb = nextShown(v, bb, b, side === 0 ? s.b1[k] : b0, dir * db);
      if (eb < 0) return false;
      t = Math.abs(nb) > 0 ? Math.round((Math.abs(eb - b0) * n) / Math.abs(nb)) : t;
      if (t === t0) break;
    }
    ends[side] = t;
  }
  return ends[0] <= ends[1];
}
function segment(v: DrawView, s: DotPlotState, k: number): boolean {
  const ba = s.aBlock[k];
  const bb = s.bBlock[k];
  let p0 = v.shownOf(ba, s.a0[k]);
  let p1 = v.shownOf(ba, s.a1[k]);
  let q0 = v.shownOf(bb, s.b0[k]);
  let q1 = v.shownOf(bb, s.b1[k]);
  if (p0 < 0 || p1 < 0 || q0 < 0 || q1 < 0) {
    // An end filtered out: the stretch from its first step shown on both axes to its last (a filter's gap closed up).
    if (!shownSteps(v, s, k)) return false;
    const na = s.a1[k] - s.a0[k];
    const nb = s.b1[k] - s.b0[k];
    const n = Math.max(Math.abs(na), Math.abs(nb)) || 1;
    p0 = v.shownOf(ba, s.a0[k] + Math.round((ends[0] * na) / n));
    p1 = v.shownOf(ba, s.a0[k] + Math.round((ends[1] * na) / n));
    q0 = v.shownOf(bb, s.b0[k] + Math.round((ends[0] * nb) / n));
    q1 = v.shownOf(bb, s.b0[k] + Math.round((ends[1] * nb) / n));
    if (p0 < 0 || p1 < 0 || q0 < 0 || q1 < 0) return false;
  }
  const dpr = v.dpr;
  const wa = v.span[ba] * v.cell;
  const wb = v.span[bb] * v.cell;
  const xa = v.blockX0[ba];
  const xb = v.blockX0[bb];
  const area = v.right - v.left;
  const sy = area > 0 ? v.height / area : 0;
  const left = v.left * dpr;
  const top = v.top * dpr;
  seg[0] = Math.round((xa + (p1 >= p0 ? p0 : p0 + 1) * wa) * dpr) - v.scrollDev;
  seg[2] = Math.round((xa + (p1 >= p0 ? p1 + 1 : p1) * wa) * dpr) - v.scrollDev;
  seg[1] = Math.round((Math.round((xb + (q1 >= q0 ? q0 : q0 + 1) * wb) * dpr) - v.scrollDev - left) * sy + top);
  seg[3] = Math.round((Math.round((xb + (q1 >= q0 ? q1 + 1 : q1) * wb) * dpr) - v.scrollDev - left) * sy + top);
  return true;
}
/**
 * `seg` clipped to the box [l, r] × [t, b] (device px) along its line (Liang–Barsky); false: it misses the box. Only
 * for a segment reaching far past the view, where an end would not fit an int32.
 */
function clipSegment(l: number, r: number, t: number, b: number): boolean {
  const dx = seg[2] - seg[0];
  const dy = seg[3] - seg[1];
  let u0 = 0;
  let u1 = 1;
  for (let e = 0; e < 4; e++) {
    const p = e === 0 ? -dx : e === 1 ? dx : e === 2 ? -dy : dy;
    const q = e === 0 ? seg[0] - l : e === 1 ? r - seg[0] : e === 2 ? seg[1] - t : b - seg[1];
    if (p === 0) {
      if (q < 0) return false;
      continue;
    }
    const u = q / p;
    if (p < 0) {
      if (u > u1) return false;
      if (u > u0) u0 = u;
    } else {
      if (u < u0) return false;
      if (u < u1) u1 = u;
    }
  }
  const x0 = seg[0];
  const y0 = seg[1];
  seg[0] = Math.round(x0 + u0 * dx);
  seg[1] = Math.round(y0 + u0 * dy);
  seg[2] = Math.round(x0 + u1 * dx);
  seg[3] = Math.round(y0 + u1 * dy);
  return true;
}
/** Past this many device px off the view a segment is clipped to it before it is drawn. */
const FAR = 1 << 24;

const PROGRESSIVE = BackendCaps.Progressive;
/** The draw (seqspec:format#layer-types `dot_plot`): grids, then segments over them. */
export function drawDotPlot(v: DrawView, s: DotPlotState): void {
  const ax = v.axis;
  if (!ax || v.geometry !== AxisGeometry.linear || v.right <= v.left || v.height <= 0) return;
  const c = tilesOf(v, s);
  const set = tileSet(c, v.backend.images);
  const progressive = (v.backend.caps & PROGRESSIVE) !== 0;
  // One count and budget a composer frame, whatever lanes draw in it: a progressive backend (the screen) builds a
  // number of samples a frame; an export builds all it draws. Past the cap, tiles not drawn lately go first.
  if (c.at !== v.frameId) {
    c.at = v.frameId;
    c.frame++;
    c.budget = DotPlotLimits.samplesPerFrame;
    if (set.count > DotPlotLimits.tiles) evict(c, set);
  }
  if (!progressive) c.budget = 0x7fffffff;
  // On a stacked section a lane is one block: its own grids alone.
  const only = ax.stacked ? v.onlyBlock : -1;
  for (let i = 0; i < s.grids.length; i++) {
    const g = s.grids[i];
    if (only < 0 || (g.bx === only && g.by === only)) drawGrid(v, c, set, g);
  }
  if (s.count === 0) return;
  const keys = keysOf(v, s, c);
  const L = Math.floor(v.left * v.dpr);
  const R = Math.ceil(v.right * v.dpr);
  const T = Math.floor(v.top * v.dpr);
  const B = Math.ceil((v.top + v.height) * v.dpr);
  const width = Math.max(2, Math.round(2 * v.dpr));
  paths.reset();
  dots.reset();
  circles.reset();
  for (let k = 0; k < s.count; k++) {
    if (only >= 0 && (s.aBlock[k] !== only || s.bBlock[k] !== only)) continue;
    if (!segment(v, s, k)) continue;
    if ((seg[0] < L && seg[2] < L) || (seg[0] > R && seg[2] > R) || (seg[1] < T && seg[3] < T)) continue;
    if (seg[1] > B && seg[3] > B) continue;
    // An end far past the view: the segment cut to the view (a margin each side), along its line.
    const far =
      seg[0] < L - FAR ||
      seg[0] > R + FAR ||
      seg[2] < L - FAR ||
      seg[2] > R + FAR ||
      seg[1] < T - FAR ||
      seg[1] > B + FAR ||
      seg[3] < T - FAR ||
      seg[3] > B + FAR;
    if (far && !clipSegment(L - width, R + width, T - width, B + width)) continue;
    const x0 = seg[0] | 0;
    const y0 = seg[1] | 0;
    const x1 = seg[2] | 0;
    const y1 = seg[3] | 0;
    const key = keys[k];
    const w = Math.abs(x1 - x0);
    const h = Math.abs(y1 - y0);
    // A pair of single positions is a dot at every zoom: its cell's rect below a device px either way, else a circle
    // centered in the cell, the viewer's dot within the cell's smaller side. A stretch too small to stroke: its rect.
    const single = s.a0[k] === s.a1[k] && s.b0[k] === s.b1[k];
    if ((single && (w < 2 || h < 2)) || (!single && w <= 3 && h <= 3)) {
      dots.begin(inkOf(key), keyAlpha(key));
      dots.push(Math.min(x0, x1), Math.min(y0, y1), Math.max(1, w), Math.max(1, h));
      continue;
    }
    if (single) {
      circles.begin(inkOf(key), keyAlpha(key));
      circles.push((x0 + x1) / 2, (y0 + y1) / 2, Math.min(DOT_R * v.dpr, Math.min(w, h) / 2), CIRCLE);
      continue;
    }
    paths.begin(inkOf(key), keyAlpha(key), PathMode.stroke, width);
    paths.moveTo(x0, y0);
    paths.lineTo(x1, y1);
  }
  if (paths.n > 0) v.backend.paths(paths);
  if (dots.n > 0) v.backend.rects(dots);
  if (circles.n > 0) v.backend.glyphs(circles);
}

// ---- The overlay -----------------------------------------------------------------------------------------------------

const BITS = [MarkBit.selection, MarkBit.highlight, MarkBit.accent, MarkBit.hover];
const TOKENS = [SkinToken.selection, SkinToken.highlight, SkinToken.accent, SkinToken.hover];
const SELECTION = MarkBit.selection;
const ACCENT = MarkBit.accent;
const BORDER = SkinToken.selectionBorder;
const fills = new RectBatch(32);
const edges = new RectBatch(16);
let shown = new Int32Array(16);
const el = new Int32Array(7);
/** A row band over shown [k0, k1) of block b across the data area, in channel `bit`'s tint (a selection: its edges). */
function row(v: DrawView, b: number, k0: number, k1: number, q: number) {
  const y0 = dotPlotY(v, b, k0);
  const y1 = Math.max(y0 + 1, dotPlotY(v, b, k1));
  const x0 = Math.floor(v.left * v.dpr);
  const x1 = Math.ceil(v.right * v.dpr);
  fills.begin(TOKENS[q], BITS[q] === ACCENT ? 63 : 255);
  fills.push(x0, y0, x1 - x0, y1 - y0);
  if (BITS[q] !== SELECTION || y1 - y0 < 3) return;
  const w = Math.max(1, Math.round(v.dpr)) | 0;
  edges.begin(BORDER);
  edges.push(x0, y0, x1 - x0, w);
  edges.push(x0, y1 - w, x1 - x0, w);
}

/**
 * Every channel's lit runs as rows across the plot (the composer draws their columns), the rectangle being dragged
 * (`v.pending` on this layer), and the hovered cell outlined (seqspec:interaction#hover). The track's first dot plot
 * draws them, once for its track.
 */
export function overlayDotPlot(v: DrawView, s: DotPlotState): void {
  const ax = v.axis;
  if (!ax || v.geometry !== AxisGeometry.linear || v.right <= v.left) return;
  if (s.first) lit(v, s);
  outlineCell(v, s);
}
const span = new Int32Array(2);
/**
 * A drag's shown [min(s, e), max(s, e) + 1) cut to one position past a block's view [lo, hi) each way into `out`, both
 * ends (LC-10: one side each left a rectangle wholly off the view reversed); false: none of it in that reach.
 */
export function pendingSpan(lo: number, hi: number, s: number, e: number, out: Int32Array): boolean {
  const a = s < e ? s : e;
  const z = (s < e ? e : s) + 1;
  out[0] = a < lo - 1 ? lo - 1 : a > hi + 1 ? hi + 1 : a;
  out[1] = z < lo - 1 ? lo - 1 : z > hi + 1 ? hi + 1 : z;
  return out[1] > out[0];
}
/** The track's first dot plot's part of the overlay: the channels' rows and the rectangle being dragged. */
function lit(v: DrawView, _s: DotPlotState) {
  const ax = v.axis as NonNullable<AxisView["axis"]>;
  const only = ax.stacked ? v.onlyBlock : -1;
  fills.reset();
  edges.reset();
  for (let q = 0; q < BITS.length; q++) {
    const set = v.marks.runs(BITS[q]);
    if (!set || set.count === 0) continue;
    for (let b = 0; b < ax.count; b++) {
      if (only >= 0 && b !== only) continue;
      visibleRange(v, b, ys);
      for (let i = set.from(b); i < set.to(b) && ys[0] < ys[1]; i++) {
        let n = ax.runToShown(b, set.start[i], set.stop[i], shown);
        if (2 * n > shown.length) {
          shown = new Int32Array(4 * n);
          n = ax.runToShown(b, set.start[i], set.stop[i], shown);
        }
        for (let j = 0; j < n; j++) {
          const k0 = Math.max(ys[0], shown[2 * j]);
          const k1 = Math.min(ys[1], shown[2 * j + 1]);
          if (k1 > k0) row(v, b, k0, k1, q);
        }
      }
    }
  }
  const p = v.pending;
  if (p?.active && p.plotLayer === v.layer && p.yBlock >= 0 && p.block >= 0) {
    // Cut to a position past the view each way: `xd` holds within a cell of the view only (VL-12); wholly off, none.
    const across = pendingSpan(v.visLo[p.block], v.visHi[p.block], p.startShown, p.endShown, span);
    const a0 = span[0];
    const a1 = span[1];
    const down = pendingSpan(v.visLo[p.yBlock], v.visHi[p.yBlock], p.yStartShown, p.yEndShown, span);
    const b0 = span[0];
    const b1 = span[1];
    if (across && down) {
      const x0 = v.xd(p.block, a0);
      const y0 = dotPlotY(v, p.yBlock, b0);
      const w = Math.max(1, v.xd(p.block, a1) - x0);
      const h = Math.max(1, dotPlotY(v, p.yBlock, b1) - y0);
      if (!p.remove) {
        fills.begin(SkinToken.selection, 255);
        fills.push(x0, y0, w, h);
      }
      edges.begin(BORDER, 255, Math.max(1, Math.round(v.dpr)));
      edges.push(x0, y0, w, h);
    }
  }
  if (fills.n > 0) v.backend.rects(fills);
  if (edges.n > 0) v.backend.rects(edges);
}

/** The overlay frame and entry top the hovered cell was last outlined in: once a track entry. */
const outlined = new Float64Array(2);
/**
 * The hovered cell (the hover's own hit on this track, its partner's block and index) outlined once for the track's
 * entry, in the ink of its color in the first dot plot (in layer order) with a grid on its block pair; where none has
 * one, by the last dot plot, in the default ink.
 */
function outlineCell(v: DrawView, s: DotPlotState) {
  if (!v.axis || !hoverEmphasis(v) || !v.marks.element(v.layer, el) || el[6] < 0) return;
  if (outlined[0] === v.frameId && outlined[1] === v.top) return;
  const bx = el[0];
  const by = el[6];
  const kx = bx >= 0 ? v.shownOf(bx, el[1]) : -1;
  const ky = by >= 0 ? v.shownOf(by, el[5]) : -1;
  if (kx < 0 || ky < 0) return;
  const at = gridAt(s, bx, by);
  if (at < 0 && !s.last) return;
  outlined[0] = v.frameId;
  outlined[1] = v.top;
  const value = dotPlotValue(s, bx, el[1], by, el[5]);
  const key = Number.isNaN(value) || at < 0 ? NO_COLOR : cellKey(v, s.grids[at], kx, value);
  const x0 = v.xd(bx, kx);
  const y0 = dotPlotY(v, by, ky);
  outlineItem(v, key, x0, y0, Math.max(1, v.xd(bx, kx + 1) - x0), Math.max(1, dotPlotY(v, by, ky + 1) - y0));
}
/** The grid of block pair (bx, by); -1: none. */
export function gridAt(s: DotPlotState, bx: number, by: number): number {
  return bx < 0 || by < 0 ? -1 : (s.byPair.get(bx * s.blocks + by) ?? -1);
}
/** The grid's value at block bx's index i and by's index j; NaN: none. */
export function dotPlotValue(s: DotPlotState, bx: number, i: number, by: number, j: number): number {
  const at = gridAt(s, bx, by);
  if (at < 0) return Number.NaN;
  const g = s.grids[at];
  const x = i - g.x0;
  const y = j - g.y0;
  return x >= 0 && x < g.nx && y >= 0 && y < g.ny ? g.values[x * g.ny + y] : Number.NaN;
}
