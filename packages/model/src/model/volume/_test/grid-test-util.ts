/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Tensor } from '@molstar/core/math/linear-algebra';
import type { Grid } from '../grid.js';

export const Dims = [13, 11, 9];
export const Orders = [
  [0, 1, 2],
  [0, 2, 1],
  [1, 0, 2],
  [1, 2, 0],
  [2, 0, 1],
  [2, 1, 0],
];

/** Smooth, not symmetric in any axis, with gradients that turn */
export function value(i: number, j: number, k: number) {
  return Math.sin(0.45 * i + 0.2 * j) * Math.cos(0.3 * j - 0.25 * k) + 0.02 * i * k + 0.1 * Math.sin(0.9 * k);
}

export function grid(
  order: number[],
  ctor: Tensor.ArrayCtor = Float32Array,
  dims = Dims,
  f: (i: number, j: number, k: number) => number = value,
): Grid {
  const space = Tensor.Space(dims, order, ctor);
  const data = space.create();
  for (let k = 0; k < dims[2]; k++) {
    for (let j = 0; j < dims[1]; j++) {
      for (let i = 0; i < dims[0]; i++) space.set(data, i, j, k, f(i, j, k));
    }
  }
  return {
    transform: { kind: 'matrix', matrix: Mat4.identity() },
    cells: Tensor.create(space, data),
    stats: { min: -1, max: 1, mean: 0, sigma: 1 },
  };
}

/** The same grid with a copied space: no strides, read through `get` and `dataOffset` */
export function viaAccessors(g: Grid): Grid {
  return { ...g, cells: Tensor.create({ ...g.cells.space }, g.cells.data) };
}

/** FNV-1a over the bit patterns of the values as float64 */
export function hash(values: ArrayLike<number>) {
  const bits = new Uint32Array(Float64Array.from(values).buffer);
  let h = 0x811c9dc5;
  for (let i = 0; i < bits.length; i++) h = Math.imul(h ^ bits[i], 0x01000193) >>> 0;
  return h.toString(16).padStart(8, '0');
}
