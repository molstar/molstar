/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Cylinder, type CylinderProps, DefaultCylinderProps } from '../../primitive/cylinder.js';
import type { Primitive } from '../../primitive/primitive.js';

/** The cylinder as it was built with a seam column and a center vertex per cap face (three.js, for UVs) */
function seamCylinder(props?: CylinderProps): Primitive {
  const p = { ...DefaultCylinderProps, ...props };
  const { radiusTop, radiusBottom, height, radialSegments, heightSegments, thetaStart, thetaLength } = p;
  const indices: number[] = [];
  const vertices: number[] = [];
  const normals: number[] = [];
  let index = 0;
  const halfHeight = height / 2;

  const rows: number[][] = [];
  const normal = Vec3.zero();
  const slope = (radiusBottom - radiusTop) / height;
  for (let y = 0; y <= heightSegments; ++y) {
    const row: number[] = [];
    const v = y / heightSegments;
    const radius = v * (radiusBottom - radiusTop) + radiusTop;
    for (let x = 0; x <= radialSegments; ++x) {
      const theta = (x / radialSegments) * thetaLength + thetaStart;
      vertices.push(radius * Math.sin(theta), -v * height + halfHeight, radius * Math.cos(theta));
      Vec3.normalize(normal, Vec3.set(normal, Math.sin(theta), slope, Math.cos(theta)));
      normals.push(...normal);
      row.push(index++);
    }
    rows.push(row);
  }
  for (let x = 0; x < radialSegments; ++x) {
    for (let y = 0; y < heightSegments; ++y) {
      const [a, b, c, d] = [rows[y][x], rows[y + 1][x], rows[y + 1][x + 1], rows[y][x + 1]];
      indices.push(a, b, d, b, c, d);
    }
  }

  const cap = (top: boolean) => {
    const radius = top ? radiusTop : radiusBottom;
    const sign = top ? 1 : -1;
    const centerStart = index;
    for (let x = 1; x <= radialSegments; ++x) {
      vertices.push(0, halfHeight * sign, 0);
      normals.push(0, sign, 0);
      ++index;
    }
    const centerEnd = index;
    for (let x = 0; x <= radialSegments; ++x) {
      const theta = (x / radialSegments) * thetaLength + thetaStart;
      vertices.push(radius * Math.sin(theta), halfHeight * sign, radius * Math.cos(theta));
      normals.push(0, sign, 0);
      ++index;
    }
    for (let x = 0; x < radialSegments; ++x) {
      const c = centerStart + x;
      const i = centerEnd + x;
      if (top) indices.push(i, i + 1, c);
      else indices.push(i + 1, i, c);
    }
  };
  if (p.topCap && radiusTop > 0) cap(true);
  if (p.bottomCap && radiusBottom > 0) cap(false);

  return {
    vertices: new Float32Array(vertices),
    normals: new Float32Array(normals),
    indices: new Uint32Array(indices),
  };
}

/** Per triangle corner, in order: what the rasterizer sees */
function soup(p: Primitive) {
  const positions: number[] = [];
  const normals: number[] = [];
  for (let i = 0; i < p.indices.length; ++i) {
    const k = p.indices[i] * 3;
    positions.push(p.vertices[k], p.vertices[k + 1], p.vertices[k + 2]);
    normals.push(p.normals[k], p.normals[k + 1], p.normals[k + 2]);
  }
  return { positions, normals };
}

function maxAbsDiff(a: number[], b: number[]) {
  let max = 0;
  for (let i = 0; i < a.length; ++i) max = Math.max(max, Math.abs(a[i] - b[i]));
  return max;
}

const cases: [string, CylinderProps][] = [
  ['default', {}],
  ['capped, 16 segments', { radiusTop: 0.1, radiusBottom: 0.1, radialSegments: 16, topCap: true, bottomCap: true }],
  [
    'cone, 7 segments, 3 rows',
    { radiusTop: 0, radiusBottom: 0.8, radialSegments: 7, heightSegments: 3, bottomCap: true },
  ],
  ['tapered, rotated start', { radiusTop: 0.5, radiusBottom: 1.2, radialSegments: 12, thetaStart: 0.4, topCap: true }],
  [
    '3 rows, both caps',
    { radiusTop: 0.7, radiusBottom: 0.4, radialSegments: 10, heightSegments: 3, topCap: true, bottomCap: true },
  ],
  ...[1, 2, 3].map((r): [string, CylinderProps] => [
    `${r} segments`,
    { radialSegments: r, topCap: true, bottomCap: true },
  ]),
];

const open: [string, CylinderProps][] = [
  ['half, capped', { radialSegments: 9, thetaLength: Math.PI, topCap: true, bottomCap: true }],
  ['more than a full turn', { radialSegments: 9, thetaLength: Math.PI * 2.5, topCap: true }],
  ['just off a full turn', { radialSegments: 9, thetaLength: Math.PI * 2 - 1e-12, topCap: true, bottomCap: true }],
];

describe('Cylinder primitive', () => {
  for (const [name, props] of cases) {
    it(`${name}: the same triangles as with the seam column and per-face cap centers`, () => {
      const actual = soup(Cylinder(props));
      const expected = soup(seamCylinder(props));
      expect(actual.positions.length).toBe(expected.positions.length);
      // only the seam differs: sin(2π) is -2.4e-16, not 0, and the closed cylinder reuses the θ = 0 column
      expect(maxAbsDiff(actual.positions, expected.positions)).toBeLessThanOrEqual(1e-15);
      expect(maxAbsDiff(actual.normals, expected.normals)).toBeLessThanOrEqual(1e-15);
    });
  }

  for (const [name, props] of open) {
    it(`${name}: keeps the seam column, the same triangles bit for bit`, () => {
      const actual = soup(Cylinder(props));
      const expected = soup(seamCylinder(props));
      expect(actual.positions).toEqual(expected.positions);
      expect(actual.normals).toEqual(expected.normals);
    });
  }

  it('every vertex is used', () => {
    for (const [, props] of [...cases, ...open]) {
      const p = Cylinder(props);
      const used = new Uint8Array(p.vertices.length / 3);
      for (let i = 0; i < p.indices.length; ++i) used[p.indices[i]] = 1;
      expect(used.every((u) => u === 1)).toBe(true);
    }
  });

  it('vertex counts: 16 capped segments 66 instead of 100, a half turn keeps its seam', () => {
    const capped = { radialSegments: 16, topCap: true, bottomCap: true };
    expect(seamCylinder(capped).vertices.length / 3).toBe(100);
    expect(Cylinder(capped).vertices.length / 3).toBe(66);
    // 2 rows of 10, 2 caps of 1 + 10
    expect(Cylinder({ ...capped, radialSegments: 9, thetaLength: Math.PI }).vertices.length / 3).toBe(42);
  });
});
