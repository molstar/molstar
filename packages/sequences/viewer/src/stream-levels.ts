/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The levels a view asks and draws (seqspec:sources#levels-and-the-zoom-gate): the level grid's constants, an alignment's
// fine and coarse levels for a view, snapped to those its source serves (`listedLevels`), and a tile's sample rows
// computed in memory (`reduce`). A leaf: it imports core and types only, so any module may import it first.

import { TILE, type Tile, Tiles } from "@molstar/seq-core/format";
import type { AxisView } from "./axis.ts";

/** Levels per axis; label chunks are a level of their own; a 1-D system's levels from `D1`, `LV` a system (data first). */
export const L = 16;
export const LABELS = L * L;
export const D1 = 2 * LABELS;
export const LV = 32;

/**
 * Fine and coarse levels for a region of `cols` × `rows` positions, `cellPx` and `rowPx` device px each (`inp`: cellPx,
 * rowPx, cols, rows, budget; a double passed to a call V8 does not inline is boxed): into `out` [ly, lx, lyc, lxc].
 * One rounding rule: the finest level whose sample spans a device pixel; then coarser while the samples exceed `budget`
 * (the axis with more gives first; 0: no budget); coarse: the region in one tile per axis.
 */
export function tileLevels(inp: Float64Array, out: Int32Array) {
  const cols = inp[2];
  const rows = inp[3];
  let lx = inp[0] >= 1 ? 0 : Math.min(L - 1, Math.ceil(Math.log2(1 / inp[0])));
  let ly = inp[1] >= 1 ? 0 : Math.min(L - 1, Math.ceil(Math.log2(1 / inp[1])));
  for (let g = 0; g < 2 * L; g++) {
    const sx = Math.ceil(cols / 2 ** lx);
    const sy = Math.ceil(rows / 2 ** ly);
    if (inp[4] <= 0 || sx * sy <= inp[4]) break;
    if ((sx >= sy && lx < L - 1) || ly >= L - 1) lx++;
    else ly++;
  }
  out[0] = ly;
  out[1] = lx;
  out[2] = Math.min(L - 1, Math.max(ly, Tiles.levelFor(rows)));
  out[3] = Math.min(L - 1, Math.max(lx, Tiles.levelFor(cols)));
}

/**
 * The levels `lv` [ly, lx, lyc, lxc] chose, each pair snapped to the one the source serves in its place (`served`, a
 * `StreamView`'s; seqspec:sources#levels-and-the-zoom-gate): the manager asks, and the draws read, the same.
 */
export function listedLevels(served: Int16Array | null | undefined, lv: Int32Array) {
  if (!served) return;
  for (let k = 0; k < 4; k += 2) {
    const s = served[lv[k] * L + lv[k + 1]];
    if (s < 0) continue;
    lv[k] = (s / L) | 0;
    lv[k + 1] = s % L;
  }
}
/**
 * Per grid level, the level a manifest serving the levels `grid` lists (>= 0) asks in its place: the coarsest it serves
 * whose sample is no wider on either axis, else the one coarser by the fewest levels (then finer by the fewest); -1
 * where it serves none.
 */
export function servedLevels(grid: Int16Array, out: Int16Array) {
  out.fill(-1);
  for (let ly = 0; ly < L; ly++)
    for (let lx = 0; lx < L; lx++) {
      let best = -1;
      let score = Number.POSITIVE_INFINITY;
      for (let a = 0; a < L; a++)
        for (let x = 0; x < L; x++) {
          if (grid[a * L + x] < 0) continue;
          // Finer or equal on both axes first, the coarsest of them; else the fewest levels coarser, then the fewest
          // finer.
          const finer = Math.max(0, ly - a) + Math.max(0, lx - x);
          const coarser = Math.max(0, a - ly) + Math.max(0, x - lx);
          const s = coarser === 0 ? finer - 4096 : 64 * coarser + finer;
          if (s < score) [score, best] = [s, a * L + x];
        }
      out[ly * L + lx] = best;
    }
}

const li = new Float64Array(5);
/**
 * Block b's levels for the lanes in view (`lanes` from `lanesInView`; a draw's `DrawView.lanes`) into `lv`, the
 * manager's and the draws' alike: its shown positions in view at its cell, under its share of `budget` by shown
 * positions, at the levels the source serves (`served`, the layer's `StreamView.served`).
 */
export function blockLevels(
  v: AxisView,
  b: number,
  budget: number,
  lanes: Float64Array,
  lv: Int32Array,
  served?: Int16Array | null,
) {
  let shown = 0;
  for (let k = 0; k < v.visHi.length; k++) shown += Math.max(0, v.visHi[k] - v.visLo[k]);
  li[0] = v.cell * v.span[b] * v.dpr;
  li[1] = lanes[5] * v.dpr;
  li[2] = v.visHi[b] - v.visLo[b];
  li[3] = lanes[1] - lanes[0];
  li[4] = shown > 0 && budget > 0 ? Math.max(1, (budget * li[2]) / shown) : budget;
  tileLevels(li, lv);
  listedLevels(served, lv);
}

/** A miniature's overview column in a section without one (CSS px across; the layout's `OVERVIEW_WIDTH`). */
export const MAP_WIDTH = 72;
/**
 * The vertical overview's level (every row and column in a tile or two per axis) into lv[2], lv[3]: a sample per device
 * px of the section's overview column (`vnavWidth` across the scrolled region; none: `MAP_WIDTH`, for a miniature
 * naming the rows track), the manager's and the rows' map's alike, whatever box the map is drawn into.
 */
export function overviewLevels(
  v: AxisView,
  lanes: number,
  columns: number,
  budget: number,
  lv: Int32Array,
  served?: Int16Array | null,
) {
  li[0] = ((v.vnavWidth > 0 ? v.vnavWidth : MAP_WIDTH) * v.dpr) / columns;
  li[1] = ((v.regionBottom - v.regionTop) * v.dpr) / lanes;
  li[2] = columns;
  li[3] = lanes;
  li[4] = budget;
  tileLevels(li, lv);
  listedLevels(served, lv);
}

const counts = new Uint32Array(256);
const seen = new Uint8Array(256);
/**
 * Sample rows [y0, y1) of a tile of the pyramid in the wire format (seqspec:sources#the-tiles-protocol): per sample of
 * 2^ly lanes × 2^lx columns the most frequent residue class (ties to the lower code; a gap only where it has no
 * residue) and its share of all the sample's cells × 255, so gappy regions fade; 0 where no data. Each cell read once;
 * the cells read.
 */
export function reduce(
  codes: Uint8Array,
  slot: Int32Array,
  C: number,
  lanes: number,
  ly: number,
  lx: number,
  ty: number,
  tx: number,
  y0: number,
  y1: number,
  out: Tile,
): number {
  const sy = 2 ** ly;
  const sx = 2 ** lx;
  let cells = 0;
  for (let y = y0; y < y1; y++) {
    const k0 = (ty * TILE + y) * sy;
    const k1 = Math.min(k0 + sy, lanes);
    for (let x = 0; x < TILE; x++) {
      const c0 = (tx * TILE + x) * sx;
      const c1 = Math.min(c0 + sx, C);
      let n = 0;
      let m = 0;
      for (let k = k0; k < k1 && c0 < c1; k++) {
        const s = slot[k];
        if (s < 0) continue;
        for (let c = s * C + c0, e = s * C + c1; c < e; c++) {
          const code = codes[c];
          if (counts[code]++ === 0) seen[m++] = code;
        }
        n += c1 - c0;
      }
      let best = 0;
      let gap = 0;
      for (let j = 0; j < m; j++) {
        const code = seen[j];
        const k = counts[code];
        if (code === 45 || code === 46) {
          if (gap === 0 || k > counts[gap] || (k === counts[gap] && code < gap)) gap = code;
        } else if (code !== 0 && (best === 0 || k > counts[best] || (k === counts[best] && code < best))) best = code;
      }
      if (best === 0) best = gap;
      out.cls[y * TILE + x] = n > 0 ? best : 0;
      out.share[y * TILE + x] = n > 0 ? Math.round((counts[best] / n) * 255) : 0;
      for (let j = 0; j < m; j++) counts[seen[j]] = 0;
      cells += n;
    }
  }
  return cells;
}
