/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { LocationIterator } from '../../util/location-iterator.js';
import { LinesBuilder } from '../lines/lines-builder.js';
import { Lines } from '../lines/lines.js';
import { Mesh } from '../mesh/mesh.js';
import {
  createPosesTransform,
  createTransform,
  getInstanceMatricesById,
  getInstanceRows,
  type TransformData,
} from '../transform-data.js';

/**
 * Per vertex and instance data (`createVertexInstanceColor`, `createVertexInstanceSize`) is written a slot an
 * instance by the position iterator's instance index and read in the shaders by `aInstance`, as the geo exports read
 * `tColor` and `tSize`: instance `i` must be placed by instance `i`'s transform, whatever order the instance grid
 * stores the rows in (`getInstanceRows`)
 */
describe('instances reordered by their grid, read by instance', () => {
  const count = 40;
  /** instances on a line in a scrambled order: the grid's cells reorder them */
  const original = new Float32Array(count * 16);
  for (let i = 0; i < count; ++i) {
    Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(((i * 17) % count) * 10, 0, 0)), original, i * 16);
  }
  const invariant = Sphere3D.create(Vec3.create(1, 0, 0), 1);
  const reordered = (data: TransformData) => {
    const order = data.aInstance.ref.value;
    expect(Array.from(order)).not.toEqual(Array.from({ length: count }, (_, i) => i));
    expect(getInstanceRows(data)).toBeDefined();
  };
  /** each instance's location its group's vertex (`local`) placed by `matrix(instance)` */
  const expectPlaced = (it: LocationIterator, local: Vec3 | ((group: number) => Vec3), matrix: (i: number) => Mat4) => {
    const expected = Vec3();
    while (it.hasNext) {
      const { location, instanceIndex, groupIndex } = it.move();
      const position = (location as { position: Vec3 }).position;
      Vec3.transformMat4(expected, typeof local === 'function' ? local(groupIndex) : local, matrix(instanceIndex));
      for (let k = 0; k < 3; ++k) expect(position[k], `instance ${instanceIndex}`).toBeCloseTo(expected[k], 4);
    }
  };
  const matrixOf = (i: number) => Mat4.fromArray(Mat4(), original, i * 16);

  it("places a mesh's vertices by its instance's matrix", () => {
    const mesh = Mesh.create(
      Float32Array.of(1, 0, 0),
      Uint32Array.of(0, 0, 0),
      Float32Array.of(0, 0, 1),
      Uint32Array.of(0),
      1,
      1,
    );
    const data = createTransform(original.slice(), count, invariant, 5, 4);
    reordered(data);
    expectPlaced(Mesh.Utils.createPositionIterator(mesh, data), Vec3.create(1, 0, 0), matrixOf);
    // the exports' matrices in the instances' order
    const byId = getInstanceMatricesById(data);
    for (let i = 0; i < count * 16; ++i) expect(byId[i]).toBe(original[i]);
  });

  it("places lines' ends by their instance's matrix (the matrices' readers)", () => {
    const builder = LinesBuilder.create();
    builder.add(1, 0, 0, 1, 1, 0, 0);
    const lines: Lines = builder.getLines();
    const data = createTransform(original.slice(), count, invariant, 5, 4);
    reordered(data);
    const ends = [Vec3.create(1, 0, 0), Vec3.create(1, 1, 0)];
    expectPlaced(Lines.Utils.createPositionIterator(lines, data), (g) => ends[g], matrixOf);
  });

  it("places a mesh's vertices by its instance's pose", () => {
    const positions = new Float32Array(count * 3);
    const rotations = new Float32Array(count * 4);
    for (let i = 0; i < count; ++i) {
      positions.set([((i * 17) % count) * 10, i % 3, 0], i * 3);
      Quat.toArray(Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), i * 0.3), rotations, i * 4);
    }
    const poses = InstanceTransforms.ofPoses({ positions, rotations, count, first: 0, center: Vec3() });
    const mesh = Mesh.create(
      Float32Array.of(1, 0, 0),
      Uint32Array.of(0, 0, 0),
      Float32Array.of(0, 0, 1),
      Uint32Array.of(0),
      1,
      1,
    );
    const data = createPosesTransform(poses, invariant, 5, 4);
    reordered(data);
    expectPlaced(Mesh.Utils.createPositionIterator(mesh, data), Vec3.create(1, 0, 0), (i) =>
      InstanceTransforms.getMatrix(poses, i, Mat4()),
    );
  });
});
