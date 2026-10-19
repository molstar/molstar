/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { Grid } from '../grid.js';
import { Dims, Orders, grid, hash, viaAccessors } from './grid-test-util.js';

/** Gradients in (i, j, k) order, whatever the layout */
function gradients(g: Grid) {
  const { values, magnitude } = Grid.getGradients(g);
  const { dataOffset, dimensions: dims } = g.cells.space;
  const out: number[] = [];
  for (let k = 0; k < dims[2]; k++) {
    for (let j = 0; j < dims[1]; j++) {
      for (let i = 0; i < dims[0]; i++) {
        const o = dataOffset(i, j, k) * 3;
        out.push(values[o], values[o + 1], values[o + 2]);
      }
    }
  }
  return { values: hash(out), magnitude };
}

/** Interpolated gradients at seeded points, inside, on the faces and outside */
function interpolated(g: Grid) {
  const get = Grid.makeGetInterpolatedGradient(g);
  let seed = 3;
  const random = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  const p = Vec3(),
    out = Vec3();
  const values: number[] = [];
  for (let n = 0; n < 500; n++) {
    Vec3.set(p, random() * (Dims[0] + 2) - 1, random() * (Dims[1] + 2) - 1, random() * (Dims[2] + 2) - 1);
    if (n % 10 === 0) p[n % 3] = Dims[n % 3] - 1;
    Vec3.set(out, 0, 0, 0);
    const ok = get(p, out);
    values.push(ok ? 1 : 0, out[0], out[1], out[2]);
  }
  return hash(values);
}

describe('grid gradients', () => {
  it('are the same for every axis order', () => {
    const reference = gradients(grid([0, 1, 2]));
    expect(reference).toMatchInlineSnapshot(`
      {
        "magnitude": {
          "max": 0.6998230678393628,
          "mean": 0.12324878582119096,
          "min": 0.024600393648774803,
          "sigma": 0.15333258600528948,
        },
        "values": "f1b04078",
      }
    `);
    for (const order of Orders) {
      expect(Tensor.strides(grid(order).cells.space), `order ${order}`).toBeDefined();
      expect(gradients(grid(order)), `order ${order}`).toEqual(reference);
      expect(gradients(viaAccessors(grid(order))), `order ${order}, accessors`).toEqual(reference);
    }
  });

  it('read Float64 data', () => {
    const reference = gradients(grid([0, 1, 2], Float64Array));
    expect(reference).toMatchInlineSnapshot(`
      {
        "magnitude": {
          "max": 0.6998231203688433,
          "mean": 0.12324878642170492,
          "min": 0.024600404087209767,
          "sigma": 0.15333258762946228,
        },
        "values": "f26164c7",
      }
    `);
    for (const order of Orders) expect(gradients(grid(order, Float64Array)), `order ${order}`).toEqual(reference);
  });

  it('use one-sided differences on axes of size 1 and 2', () => {
    // size 1: both neighbours are the voxel itself, the divisor falls back to 1; size 2: one-sided differences
    const references = [];
    for (const dims of [
      [1, 2, 3],
      [3, 1, 2],
      [2, 3, 1],
    ]) {
      const reference = gradients(grid([0, 1, 2], Float32Array, dims));
      for (const order of Orders) {
        expect(gradients(grid(order, Float32Array, dims)), `${dims} order ${order}`).toEqual(reference);
        const g = viaAccessors(grid(order, Float32Array, dims));
        expect(gradients(g), `${dims} order ${order}, accessors`).toEqual(reference);
      }
      references.push(reference.values);
    }
    expect(references).toMatchInlineSnapshot(`
      [
        "6c48c935",
        "611e940c",
        "5f0377c7",
      ]
    `);
  });

  it('report a magnitude range of 0..1 when no gradient is a number', () => {
    for (const g of [
      grid([0, 1, 2], Float32Array, [1, 1, 1], () => Number.NaN),
      grid([2, 1, 0], Float32Array, [2, 2, 1], () => Number.NaN),
    ]) {
      const { magnitude } = Grid.getGradients(g);
      expect(magnitude.min).toBe(0);
      expect(magnitude.max).toBe(1);
      expect(Grid.getGradients(viaAccessors(g)).magnitude).toEqual(magnitude);
    }
  });

  it('read an overridden get', () => {
    const g = grid([1, 2, 0]);
    const { space, data } = g.cells;
    const doubled: Grid = {
      ...g,
      cells: Tensor.create({ ...space, get: (d, i, j, k) => 2 * space.get(d, i, j, k) }, data),
    };
    const a = Grid.getGradients(g).values;
    const b = Grid.getGradients(doubled).values;
    for (let i = 0; i < a.length; i++) expect(b[i]).toBe(2 * a[i]);
  });
});

describe('grid interpolated gradient', () => {
  it('is the same for every axis order', () => {
    const reference = interpolated(grid([0, 1, 2]));
    expect(reference).toMatchInlineSnapshot(`"e8736173"`);
    for (const order of Orders) {
      expect(interpolated(grid(order)), `order ${order}`).toEqual(reference);
      expect(interpolated(viaAccessors(grid(order))), `order ${order}, accessors`).toEqual(reference);
    }
  });
});
