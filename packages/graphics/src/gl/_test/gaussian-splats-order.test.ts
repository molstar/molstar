/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { ValueCell } from '@molstar/core/util/value-cell';
import { GaussianSplats } from '@molstar/graphics/geo/geometry/gaussian-splats/gaussian-splats';
import { setSplat } from '@molstar/graphics/geo/geometry/gaussian-splats/splat-data';
import { int32ToFloatBits } from '@molstar/graphics/geo/geometry/center-group';
import {
  createSplatOrder,
  disposeSplatOrder,
  getSplatOrderWords,
  releaseForeignSplatOrder,
  sortSplats,
  updateSplatOrder,
} from '../renderable/gaussian-splats-order.js';
import type { ComputeContext, StorageBuffer } from '../compute.js';
import type { SplatSorter } from '../compute/splat-sort.js';
import type { GaussianSplatsValues } from '../renderable/gaussian-splats.js';

/** Splats at `centers` (xyz each), their values as an object of `transforms` (one instance each) has them */
function splatValues(centers: number[][], transforms?: Mat4[]): GaussianSplatsValues {
  const { positionGroup, splatData } = GaussianSplats.createArrays(centers.length);
  centers.forEach((c, i) => {
    positionGroup.set([c[0], c[1], c[2], int32ToFloatBits(i)], i * 4);
    setSplat(splatData, i, 0, 0, 0, 0, 0, 0, 1, 1);
  });
  const splats = GaussianSplats.create(positionGroup, splatData, centers.length);
  const values = GaussianSplats.Utils.createValuesSimple(
    splats,
    PD.getDefaultValues(GaussianSplats.Params),
    ColorNames.white,
    1,
  );
  if (transforms) {
    const array = new Float32Array(transforms.length * 16);
    transforms.forEach((t, i) => Mat4.toArray(t, array, i * 16));
    ValueCell.update(values.aTransform, array);
    ValueCell.update(values.instanceCount, transforms.length);
    ValueCell.update(values.aInstance, new Uint32Array(transforms.map((_, i) => i)));
  }
  return values;
}

/** The view of a camera at `eye` looking at the origin */
function lookAt(eye: Vec3) {
  return Mat4.lookAt(Mat4(), eye, Vec3(), Vec3.unitY);
}

function order(values: GaussianSplatsValues) {
  return Array.from(getSplatOrderWords(values.tSplatOrder.ref.value));
}

describe('gaussian splats order', () => {
  it('starts as the splats are stored, per instance a block', () => {
    expect(Array.from(getSplatOrderWords(createSplatOrder(3, 2)))).toEqual([0, 1, 2, 0, 1, 2]);
  });

  it('sorts back to front: the most negative view z first, equal depths kept in order', () => {
    // centers along x; a view row of +x: z = x
    const pg = new Float32Array([2, 0, 0, 0, -1, 0, 0, 0, 5, 0, 0, 0, 2, 0, 0, 0]);
    const out = new Uint32Array(4);
    sortSplats(out, 0, pg, 4, 1, 0, 0, new Uint32Array(4));
    expect(Array.from(out)).toEqual([1, 0, 3, 2]);
    sortSplats(out, 0, pg, 4, -1, 0, 0, new Uint32Array(4));
    expect(Array.from(out)).toEqual([2, 0, 3, 1]);
  });

  it('orders the splats for the camera, again where it turns, not where it moves', () => {
    const values = splatValues([
      [0, 0, -1],
      [0, 0, 2],
      [0, 0, 0],
    ]);
    // from +z: z = -2 the farthest
    expect(updateSplatOrder(values, lookAt(Vec3.create(0, 0, 10)))).toBe(true);
    expect(order(values)).toEqual([0, 2, 1]);
    const version = values.tSplatOrder.ref.version;
    // moved along its view, or a little sideways (turned below the threshold): the order kept, not written
    expect(updateSplatOrder(values, lookAt(Vec3.create(0, 0, 30)))).toBe(false);
    expect(updateSplatOrder(values, lookAt(Vec3.create(0.1, 0, 10)))).toBe(false);
    expect(values.tSplatOrder.ref.version).toBe(version);
    // from -z: reversed
    expect(updateSplatOrder(values, lookAt(Vec3.create(0, 0, -10)))).toBe(true);
    expect(order(values)).toEqual([1, 2, 0]);
  });

  it('sorts again where the splats or their transforms change', () => {
    const values = splatValues([
      [0, 0, -1],
      [0, 0, 1],
    ]);
    const view = lookAt(Vec3.create(0, 0, 10));
    updateSplatOrder(values, view);
    expect(order(values)).toEqual([0, 1]);
    // the centers swapped
    const pg = values.splatCenters.ref.value;
    pg.array.set([0, 0, 1], 0);
    pg.array.set([0, 0, -1], 4);
    ValueCell.update(values.splatCenters, pg);
    expect(updateSplatOrder(values, view)).toBe(true);
    expect(order(values)).toEqual([1, 0]);
    // the instance turned half around
    const turned = Mat4.fromQuat(Mat4(), Quat.setAxisAngle(Quat(), Vec3.unitY, Math.PI));
    ValueCell.update(values.aTransform, Mat4.toArray(turned, new Float32Array(16), 0) as Float32Array);
    expect(updateSplatOrder(values, view)).toBe(true);
    expect(order(values)).toEqual([0, 1]);
  });

  it('sorts each instance in its block, by its own transform and its instance id', () => {
    const values = splatValues(
      [
        [0, 0, -1],
        [0, 0, 1],
      ],
      [Mat4.identity(), Mat4.fromQuat(Mat4(), Quat.setAxisAngle(Quat(), Vec3.unitY, Math.PI))],
    );
    // rows in the reverse of the instances' order (as an instance grid stores them)
    ValueCell.update(values.aInstance, new Uint32Array([1, 0]));
    expect(updateSplatOrder(values, lookAt(Vec3.create(0, 0, 10)))).toBe(true);
    // instance 1 (row 0, not turned) first-to-last as stored; instance 0 (row 1, turned) reversed
    expect(order(values)).toEqual([1, 0, 0, 1]);
  });

  it("forgets an order on the GPU of another context (a backend switch): the CPU's again, sorted whole", () => {
    const values = splatValues([
      [0, 0, -1],
      [0, 0, 1],
    ]);
    const view = lookAt(Vec3.create(0, 0, 10));
    updateSplatOrder(values, view);
    const buffer = {} as StorageBuffer;
    ValueCell.update(values.sSplatOrder, buffer);
    ValueCell.update(values.dSplatGpuOrder, true);
    // the context's own: kept
    const own = { ownsBuffer: (b: StorageBuffer) => b === buffer } as ComputeContext;
    expect(releaseForeignSplatOrder(values, own)).toBe(false);
    expect(values.dSplatGpuOrder.ref.value).toBe(true);
    // another's (or no compute: WebGL): unset
    expect(releaseForeignSplatOrder(values, undefined)).toBe(true);
    expect(values.dSplatGpuOrder.ref.value).toBe(false);
    expect(values.sSplatOrder.ref.value).toBeUndefined();
    expect(releaseForeignSplatOrder(values, undefined)).toBe(false);
    // the view unchanged, still sorted again: the state forgotten
    expect(updateSplatOrder(values, view)).toBe(true);
    expect(order(values)).toEqual([0, 1]);
  });

  it("destroys the GPU sort's buffers on dispose: the centers and the order, once", () => {
    const destroyed: string[] = [];
    const compute = {
      createBuffer: (byteLength: number, label: string) =>
        ({ byteLength, write: () => {}, destroy: () => destroyed.push(label) }) as unknown as StorageBuffer,
    } as unknown as ComputeContext;
    const gpu = { compute, sorter: { sort: () => {}, destroy: () => {} } as SplatSorter };
    const values = splatValues([
      [0, 0, -1],
      [0, 0, 1],
    ]);
    expect(updateSplatOrder(values, lookAt(Vec3.create(0, 0, 10)), gpu)).toBe(true);
    expect(values.sSplatOrder.ref.value).toBeDefined();
    expect(destroyed).toEqual([]);
    disposeSplatOrder(values);
    expect(destroyed.sort()).toEqual(['splat centers', 'splat order']);
    // forgotten: nothing destroyed twice
    disposeSplatOrder(values);
    expect(destroyed.length).toBe(2);
  });
});
