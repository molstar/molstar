/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// adapted from three.js, MIT License Copyright 2010-2018 three.js authors
// without UVs: a full turn has no seam column (it wraps to the first) and a cap has one center vertex, not one per face;
// the triangles are the same, in the same order (16 capped segments: 66 vertices instead of 100)

import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Primitive } from './primitive.js';

export const DefaultCylinderProps = {
  radiusTop: 1,
  radiusBottom: 1,
  height: 1,
  radialSegments: 8,
  heightSegments: 1,
  topCap: false,
  bottomCap: false,
  thetaStart: 0.0,
  thetaLength: Math.PI * 2,
};
export type CylinderProps = Partial<typeof DefaultCylinderProps>;

export function Cylinder(props?: CylinderProps): Primitive {
  const {
    radiusTop,
    radiusBottom,
    height,
    radialSegments,
    heightSegments,
    topCap,
    bottomCap,
    thetaStart,
    thetaLength,
  } = { ...DefaultCylinderProps, ...props };

  // buffers
  const indices: number[] = [];
  const vertices: number[] = [];
  const normals: number[] = [];

  // helper variables
  let index = 0;
  const indexArray: number[][] = [];
  const halfHeight = height / 2;
  // a full turn ends where it starts (as in `Circle`)
  const isFull = thetaLength === Math.PI * 2;
  const columnCount = isFull ? radialSegments : radialSegments + 1;

  // generate geometry
  generateTorso();

  if (topCap && radiusTop > 0) generateCap(true);
  if (bottomCap && radiusBottom > 0) generateCap(false);

  return {
    vertices: new Float32Array(vertices),
    normals: new Float32Array(normals),
    indices: new Uint32Array(indices),
  };

  function generateTorso() {
    const normal = Vec3.zero();

    // this will be used to calculate the normal
    const slope = (radiusBottom - radiusTop) / height;

    // generate vertices and normals

    for (let y = 0; y <= heightSegments; ++y) {
      const indexRow: number[] = [];
      const v = y / heightSegments;

      // calculate the radius of the current row
      const radius = v * (radiusBottom - radiusTop) + radiusTop;

      for (let x = 0; x < columnCount; ++x) {
        const u = x / radialSegments;
        const theta = u * thetaLength + thetaStart;
        const sinTheta = Math.sin(theta);
        const cosTheta = Math.cos(theta);

        // vertex
        vertices.push(radius * sinTheta, -v * height + halfHeight, radius * cosTheta);

        // normal
        Vec3.normalize(normal, Vec3.set(normal, sinTheta, slope, cosTheta));
        normals.push(...normal);

        // save index of vertex in respective row
        indexRow.push(index++);
      }
      // the last column of a full turn is the first
      if (isFull) indexRow.push(indexRow[0]);

      // now save vertices of the row in our index array
      indexArray.push(indexRow);
    }

    // generate indices

    for (let x = 0; x < radialSegments; ++x) {
      for (let y = 0; y < heightSegments; ++y) {
        // we use the index array to access the correct indices
        const a = indexArray[y][x];
        const b = indexArray[y + 1][x];
        const c = indexArray[y + 1][x + 1];
        const d = indexArray[y][x + 1];

        // faces
        indices.push(a, b, d);
        indices.push(b, c, d);
      }
    }
  }

  function generateCap(top: boolean) {
    const radius = top === true ? radiusTop : radiusBottom;
    const sign = top === true ? 1 : -1;

    // the center vertex, shared by all faces (no uvs)
    const c = index;
    vertices.push(0, halfHeight * sign, 0);
    normals.push(0, sign, 0);
    ++index;

    // save the index of the first ring vertex
    const ringStart = index;

    // now we generate the surrounding vertices and normals
    for (let x = 0; x < columnCount; ++x) {
      const u = x / radialSegments;
      const theta = u * thetaLength + thetaStart;

      const cosTheta = Math.cos(theta);
      const sinTheta = Math.sin(theta);

      // vertex
      vertices.push(radius * sinTheta, halfHeight * sign, radius * cosTheta);

      // normal
      normals.push(0, sign, 0);

      // increase index
      ++index;
    }

    // generate indices
    for (let x = 0; x < radialSegments; ++x) {
      const i = ringStart + x;
      // the last face of a full turn closes on the first ring vertex
      const j = isFull && x === radialSegments - 1 ? ringStart : i + 1;

      if (top === true) {
        indices.push(i, j, c); // face top
      } else {
        indices.push(j, i, c); // face bottom
      }
    }
  }
}
