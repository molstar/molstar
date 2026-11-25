/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GridDense, GridTable } from './protocol.js';

// Grid levels (_spec/README.md "Data", Summaries): a grid read whole once, summarized at factor 2 a side (square
// levels), `sum` or `mean` as its coverage's `aggregate` (NaN not counted: a cell with no value below is NaN); a tile
// cut from a level; a served tile transposed for its mirror; the whole-assembly level assembled from the blocks'.

/** A whole grid at one level: `nx × ny` cells from bin 0 on each axis, x-major. */
export interface GridPlane {
  readonly nx: number;
  readonly ny: number;
  readonly value: Float32Array;
}

/** A grid table as a plane over `[0, nx) × [0, ny)` bins, cells not listed `fill` (a dense tile placed at its bins). */
export function gridPlane(t: GridTable, nx: number, ny: number): GridPlane {
  const fill = t.form === 'sparse' ? t.fill : Number.NaN;
  const value = new Float32Array(nx * ny).fill(fill);
  if (t.form === 'sparse') {
    for (let i = 0; i < t.value.length; i++) {
      const x = t.xBin[i];
      const y = t.yBin[i];
      if (x >= 0 && x < nx && y >= 0 && y < ny) value[x * ny + y] = t.value[i];
    }
  } else
    for (let i = 0; i < t.nx; i++)
      for (let j = 0; j < t.ny; j++) {
        const x = t.x + i;
        const y = t.y + j;
        if (x < nx && y < ny) value[x * ny + y] = t.value[i * t.ny + j];
      }
  return { nx, ny, value };
}

/** A plane at twice the bin width a side: each cell the `sum` or `mean` of the up to four below that have a value. */
export function halveGrid(p: GridPlane, aggregate: 'sum' | 'mean'): GridPlane {
  const nx = Math.ceil(p.nx / 2);
  const ny = Math.ceil(p.ny / 2);
  const value = new Float32Array(nx * ny);
  for (let i = 0; i < nx; i++)
    for (let j = 0; j < ny; j++) {
      let s = 0;
      let n = 0;
      for (let di = 0; di < 2; di++) {
        const x = 2 * i + di;
        if (x >= p.nx) continue;
        for (let dj = 0; dj < 2; dj++) {
          const y = 2 * j + dj;
          if (y >= p.ny) continue;
          const v = p.value[x * p.ny + y];
          if (Number.isNaN(v)) continue;
          s += v;
          n++;
        }
      }
      value[i * ny + j] = n === 0 ? Number.NaN : aggregate === 'sum' ? s : s / n;
    }
  return { nx, ny, value };
}

/** The cells of bins `[x0, x1) × [y0, y1)` cut from a plane, as a dense tile (clipped to the plane). */
export function cutGrid(p: GridPlane, x0: number, x1: number, y0: number, y1: number): GridDense {
  const a = Math.max(0, Math.min(p.nx, x0));
  const b = Math.max(a, Math.min(p.nx, x1));
  const c = Math.max(0, Math.min(p.ny, y0));
  const d = Math.max(c, Math.min(p.ny, y1));
  const ny = d - c;
  const value = new Float32Array((b - a) * ny);
  for (let i = a; i < b; i++) value.set(p.value.subarray(i * p.ny + c, i * p.ny + d), (i - a) * ny);
  return { form: 'dense', x: a, y: c, nx: b - a, ny, value };
}

/** A served tile as its mirror (x and y swapped): what a reader of a symmetric grid's other triangle draws. */
export function transposeGrid(t: GridTable): GridTable {
  if (t.form === 'sparse') return { ...t, x: t.y, y: t.x, xBin: t.yBin, yBin: t.xBin };
  const value = new Float32Array(t.value.length);
  for (let i = 0; i < t.nx; i++) for (let j = 0; j < t.ny; j++) value[j * t.nx + i] = t.value[i * t.ny + j];
  return { form: 'dense', x: t.y, y: t.x, nx: t.ny, ny: t.nx, value };
}

/** A block's place in the whole-assembly system: its sequence's plane at a level, its offset there (positions). */
export interface WholeBlock {
  readonly x: number;
  readonly y: number;
  readonly plane: GridPlane;
  /** The plane's level (its bins' width). */
  readonly level: number;
}

/**
 * The whole-assembly plane at `level` (Juicebox's "All"): each block pair's plane, at a level no wider, placed by its
 * sequences' offsets, a cell going to the whole bin holding its first position (`sum` or `mean` of those meeting it).
 */
export function wholePlane(length: number, level: number, blocks: readonly WholeBlock[], aggregate: 'sum' | 'mean') {
  const n = Math.ceil(length / level);
  const sum = new Float64Array(n * n);
  const count = new Float32Array(n * n);
  for (const b of blocks) {
    const p = b.plane;
    for (let i = 0; i < p.nx; i++) {
      const x = Math.floor((b.x + i * b.level) / level);
      if (x >= n) break;
      for (let j = 0; j < p.ny; j++) {
        const v = p.value[i * p.ny + j];
        if (Number.isNaN(v)) continue;
        const y = Math.floor((b.y + j * b.level) / level);
        if (y >= n) break;
        sum[x * n + y] += v;
        count[x * n + y]++;
      }
    }
  }
  const value = new Float32Array(n * n);
  for (let k = 0; k < value.length; k++)
    value[k] = count[k] === 0 ? Number.NaN : aggregate === 'sum' ? sum[k] : sum[k] / count[k];
  return { nx: n, ny: n, value } satisfies GridPlane;
}
