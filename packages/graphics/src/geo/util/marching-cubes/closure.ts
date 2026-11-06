/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Tensor } from '@molstar/core/math/linear-algebra';

/** Fields whose surface is closed against the ignored cells (id -2), and whether a vertex is on the closure */
export function getClosureFields(scalarField: Tensor, idField: Tensor, isoLevel: number, blurPasses = 2) {
  const { dimensions, axisOrderSlowToFast } = scalarField.space;
  const [nx, ny, nz] = dimensions;
  const space = Tensor.Space(dimensions as number[], axisOrderSlowToFast as number[], Float32Array);
  const idSpace = Tensor.Space(dimensions as number[], axisOrderSlowToFast as number[], Int32Array);
  const s0 = space.dataOffset(1, 0, 0) - space.dataOffset(0, 0, 0);
  const s1 = space.dataOffset(0, 1, 0) - space.dataOffset(0, 0, 0);
  const s2 = space.dataOffset(0, 0, 1) - space.dataOffset(0, 0, 0);
  const strides = Tensor.strides(space);
  const same = (t: Tensor) => {
    const ts = Tensor.strides(t.space);
    return !!ts && !!strides && ts.every((x, i) => x === strides[i]);
  };

  const field = space.create() as unknown as Float32Array;
  const id = idSpace.create() as unknown as Int32Array;
  if (same(scalarField) && same(idField)) {
    field.set(scalarField.data);
    id.set(idField.data);
  } else {
    const { get: getValue } = scalarField.space;
    const { get: getId } = idField.space;
    for (let k = 0; k < nz; ++k) {
      for (let j = 0; j < ny; ++j) {
        for (let i = 0; i < nx; ++i) {
          const o = i * s0 + j * s1 + k * s2;
          field[o] = getValue(scalarField.data, i, j, k);
          id[o] = getId(idField.data, i, j, k);
        }
      }
    }
  }

  const label = new Float32Array(field.length);
  let ignored = 0;
  let unknown = 0;
  let maxDepth = 0;
  let i0 = nx;
  let j0 = ny;
  let k0 = nz;
  let i1 = -1;
  let j1 = -1;
  let k1 = -1;
  for (let k = 0; k < nz; ++k) {
    for (let j = 0; j < ny; ++j) {
      for (let i = 0; i < nx; ++i) {
        const o = i * s0 + j * s1 + k * s2;
        const a = id[o];
        label[o] = a === -2 ? -1 : a === -1 ? 0 : 1;
        if (a === -2) id[o] = -1;
        const depth = field[o] - isoLevel;
        if (depth <= 0) continue;
        if (a === -2) ++ignored;
        else if (a === -1) ++unknown;
        if (depth > maxDepth) maxDepth = depth;
        if (i < i0) i0 = i;
        if (i > i1) i1 = i;
        if (j < j0) j0 = j;
        if (j > j1) j1 = j;
        if (k < k0) k0 = k;
        if (k > k1) k1 = k;
      }
    }
  }
  if (ignored === 0) return { scalarField, idField, isClosure: undefined };
  i0 = Math.max(1, i0 - 1);
  j0 = Math.max(1, j0 - 1);
  k0 = Math.max(1, k0 - 1);
  i1 = Math.min(nx - 2, i1 + 1);
  j1 = Math.min(ny - 2, j1 + 1);
  k1 = Math.min(nz - 2, k1 + 1);

  if (unknown > 0) {
    const distance = new Float32Array(label.length);
    for (let k = k0 - 1; k <= k1 + 1; ++k) {
      for (let j = j0 - 1; j <= j1 + 1; ++j) {
        for (let i = i0 - 1; i <= i1 + 1; ++i) {
          const o = i * s0 + j * s1 + k * s2;
          distance[o] = label[o] !== 0 ? 0 : field[o] > isoLevel ? Infinity : -1;
        }
      }
    }
    const offsets: number[] = [];
    const weights: number[] = [];
    for (let dk = -1; dk <= 0; ++dk) {
      for (let dj = -1; dj <= 1; ++dj) {
        for (let di = -1; di <= 1; ++di) {
          if (dk < 0 || dj < 0 || (dj === 0 && di < 0)) {
            offsets.push(di * s0 + dj * s1 + dk * s2);
            weights.push(Math.sqrt(di * di + dj * dj + dk * dk));
          }
        }
      }
    }
    const relax = (o: number, sign: number) => {
      let best = distance[o];
      if (!(best > 0)) return;
      let l = label[o];
      for (let q = 0; q < 13; ++q) {
        const n = o + sign * offsets[q];
        const d = distance[n] + weights[q];
        if (d < best && distance[n] >= 0) {
          best = d;
          l = label[n];
        }
      }
      distance[o] = best;
      label[o] = l;
    };
    for (let k = k0; k <= k1; ++k) {
      for (let j = j0; j <= j1; ++j) for (let i = i0; i <= i1; ++i) relax(i * s0 + j * s1 + k * s2, 1);
    }
    for (let k = k1; k >= k0; --k) {
      for (let j = j1; j >= j0; --j) for (let i = i1; i >= i0; --i) relax(i * s0 + j * s1 + k * s2, -1);
    }
  }

  const line = new Float32Array(Math.max(nx, ny, nz));
  const blur = (
    stride: number,
    b0: number,
    b1: number,
    sa: number,
    a0: number,
    a1: number,
    sc: number,
    c0: number,
    c1: number,
  ) => {
    for (let c = c0; c <= c1; ++c) {
      for (let a = a0; a <= a1; ++a) {
        const base = a * sa + c * sc;
        for (let b = b0; b <= b1; ++b) line[b] = label[base + b * stride];
        for (let b = b0; b <= b1; ++b) {
          const prev = line[b > b0 ? b - 1 : b];
          const next = line[b < b1 ? b + 1 : b];
          label[base + b * stride] = (prev + 2 * line[b] + next) / 4;
        }
      }
    }
  };
  for (let pass = 0; pass < blurPasses; ++pass) {
    blur(s0, i0, i1, s1, j0, j1, s2, k0, k1);
    blur(s1, j0, j1, s0, i0, i1, s2, k0, k1);
    blur(s2, k0, k1, s0, i0, i1, s1, j0, j1);
  }
  for (let k = k0; k <= k1; ++k) {
    for (let j = j0; j <= j1; ++j) {
      for (let i = i0; i <= i1; ++i) {
        const o = i * s0 + j * s1 + k * s2;
        if (label[o] < 0 && field[o] > isoLevel) field[o] = 2 * isoLevel - field[o];
      }
    }
  }

  const epsilon = maxDepth * 1e-4;
  const get = scalarField.space.get;
  const data = scalarField.data;
  const isClosure = (x: number, y: number, z: number) => {
    const i = Math.min(Math.floor(x), nx - 2);
    const j = Math.min(Math.floor(y), ny - 2);
    const k = Math.min(Math.floor(z), nz - 2);
    const fx = x - i;
    const fy = y - j;
    const fz = z - k;
    const c00 = get(data, i, j, k) * (1 - fx) + get(data, i + 1, j, k) * fx;
    const c10 = get(data, i, j + 1, k) * (1 - fx) + get(data, i + 1, j + 1, k) * fx;
    const c01 = get(data, i, j, k + 1) * (1 - fx) + get(data, i + 1, j, k + 1) * fx;
    const c11 = get(data, i, j + 1, k + 1) * (1 - fx) + get(data, i + 1, j + 1, k + 1) * fx;
    const v = (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz;
    return v - isoLevel > epsilon;
  };

  return {
    scalarField: Tensor.create(space, field as unknown as Tensor.Data),
    idField: Tensor.create(idSpace, id as unknown as Tensor.Data),
    isClosure,
  };
}
