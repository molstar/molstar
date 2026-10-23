/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { Grid } from '../grid.js';
import { Dims, Orders, grid, hash, value, viaAccessors } from './grid-test-util.js';

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

describe('grid isosurface bounding sphere', () => {
  /** The spheres for positive and negative iso values, ones that no voxel reaches included */
  function spheres(g: Grid) {
    return [0.5, 0, -0.5, 0.9, -0.9, 5, -5].map((iso) => {
      const s = Grid.getIsosurfaceBoundingSphere(g, iso);
      return [s.center[0], s.center[1], s.center[2], s.radius];
    });
  }

  it('is the same for every axis order', () => {
    const reference = spheres(grid([0, 1, 2]));
    expect(hash(reference.flat())).toMatchInlineSnapshot(`"2ff5e0b3"`);
    for (const order of Orders) {
      expect(spheres(grid(order)), `order ${order}`).toEqual(reference);
      expect(spheres(viaAccessors(grid(order))), `order ${order}, accessors`).toEqual(reference);
    }
  });

  it('is cached for the last iso value, returned as a copy', () => {
    const g = grid([2, 0, 1]);
    const first = Grid.getIsosurfaceBoundingSphere(g, 0.5);
    const expected = Sphere3D.clone(first);
    first.radius = -1;
    expect(Grid.getIsosurfaceBoundingSphere(g, 0.5)).toEqual(expected);
    expect(Grid.getIsosurfaceBoundingSphere(g, 0.5, Sphere3D())).toEqual(expected);
    expect(Grid.getIsosurfaceBoundingSphere(g, -0.5)).toEqual(Grid.getIsosurfaceBoundingSphere(grid([2, 0, 1]), -0.5));
  });

  it('is not carried over to a copy with another transform', () => {
    const g = grid([0, 1, 2]);
    const s = Grid.getIsosurfaceBoundingSphere(g, 0.5);
    const moved: Grid = {
      ...g,
      transform: { kind: 'matrix', matrix: Mat4.fromTranslation(Mat4(), Vec3.create(100, 0, 0)) },
    };
    const t = Grid.getIsosurfaceBoundingSphere(moved, 0.5);
    expect(t.center[0]).toBeCloseTo(s.center[0] + 100, 5);
    expect(t.radius).toBeCloseTo(s.radius, 5);
  });
});

describe('grid trilinear interpolation', () => {
  const [mi, mj, mk] = Dims;

  it('reads the periodic image outside a periodic grid and wraps at its last voxel (perf-audit-2 04/12/B3)', () => {
    const g: Grid = { ...grid([0, 1, 2]), periodicity: 'xyz' };
    const get = Grid.makeGetTrilinearlyInterpolated(g, 'none');
    const p = Vec3.create(2.25, 3.5, 4.75);
    const inside = get(p);
    for (const [di, dj, dk] of [
      [1, 0, 0],
      [-1, 2, 0],
      [0, -1, -3],
    ]) {
      expect(get(Vec3.create(p[0] + di * mi, p[1] + dj * mj, p[2] + dk * mk))).toBeCloseTo(inside, 5);
    }
    expect(get(Vec3.create(mi - 0.5, 0, 0))).toBeCloseTo(0.5 * (value(mi - 1, 0, 0) + value(0, 0, 0)), 5);
    expect(get(Vec3.create(-0.5, 0, 0))).toBeCloseTo(0.5 * (value(mi - 1, 0, 0) + value(0, 0, 0)), 5);
  });

  it('gives NaN outside a grid that is not periodic', () => {
    const get = Grid.makeGetTrilinearlyInterpolated(grid([0, 1, 2]), 'none');
    expect(get(Vec3.create(1, 1, 1))).toBeCloseTo(value(1, 1, 1), 5);
    expect(get(Vec3.create(mi + 1, 1, 1))).toBeNaN();
    expect(get(Vec3.create(1, -2, 1))).toBeNaN();
  });
});
