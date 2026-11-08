/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { MeshBuilder } from '@molstar/graphics/geo/geometry/mesh/mesh-builder';
import { Box } from '@molstar/graphics/geo/primitive/box';
import * as Shape from '@molstar/graphics/geo/shape/shape';
import type { ShapeProvider } from '@molstar/graphics/geo/shape/provider';
import {
  getShapeProviderCenter,
  getShapeProviderTransform,
  isSameShapeProviderTransform,
  transformShapeProvider,
} from '../shape/transform.js';

/** A unit box at the origin, as one instance at each of `at` */
function instancedBoxProvider(at: Vec3[]): ShapeProvider<undefined, Mesh, Mesh.Params> {
  return {
    label: 'Boxes',
    data: undefined,
    params: Mesh.Params,
    geometryUtils: Mesh.Utils,
    getShape: () => {
      const state = MeshBuilder.createState(64, 64);
      MeshBuilder.addPrimitive(state, Mat4.identity(), Box());
      return Shape.create(
        'boxes',
        undefined,
        MeshBuilder.getMesh(state),
        () => Color(0x999999),
        () => 1,
        () => 'box',
        at.map((p) => Mat4.fromTranslation(Mat4(), p)),
      );
    },
  };
}

/** A unit box at `at`, as one instance */
function boxProvider(at: Vec3): ShapeProvider<undefined, Mesh, Mesh.Params> {
  return {
    label: 'Box',
    data: undefined,
    params: Mesh.Params,
    geometryUtils: Mesh.Utils,
    getShape: () => {
      const state = MeshBuilder.createState(64, 64);
      MeshBuilder.addPrimitive(state, Mat4.fromTranslation(Mat4(), at), Box());
      return Shape.create(
        'box',
        undefined,
        MeshBuilder.getMesh(state),
        () => Color(0x999999),
        () => 1,
        () => 'box',
      );
    },
  };
}

describe('shape transform', () => {
  it('moves every shape the provider makes by its instances, the geometry kept', async () => {
    const provider = boxProvider(Vec3.create(1, 2, 3));
    const move = Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0));
    const moved = transformShapeProvider(provider, move);
    const props = PD.getDefaultValues(provider.params);
    const before = await provider.getShape(SyncRuntimeContext, undefined, props);
    const after = await moved.getShape(SyncRuntimeContext, undefined, props);
    expect(after.transforms).toHaveLength(before.transforms.length);
    expect(Mat4.areEqual(after.transforms[0], Mat4.mul(Mat4(), move, before.transforms[0]), 1e-9)).toBe(true);
    expect(after.geometry.boundingSphere.center).toEqual(before.geometry.boundingSphere.center);
  });

  it("gives the shape's center, moved by a transform before it", async () => {
    const provider = boxProvider(Vec3.create(1, 2, 3));
    // the bounding sphere's center: about the box's
    const center = await getShapeProviderCenter(SyncRuntimeContext, provider);
    for (let k = 0; k < 3; ++k) expect(center[k]).toBeCloseTo([1, 2, 3][k], 0);
    const moved = transformShapeProvider(provider, Mat4.fromTranslation(Mat4(), Vec3.create(0, 5, 0)));
    const c = await getShapeProviderCenter(SyncRuntimeContext, moved);
    for (let k = 0; k < 3; ++k) expect(c[k]).toBeCloseTo(center[k] + [0, 5, 0][k], 5);
  });

  it("gives the center of all of a shape's instances, as a representation bounds them, not its first's", async () => {
    const provider = instancedBoxProvider([Vec3.create(0, 0, 0), Vec3.create(20, 0, 0), Vec3.create(10, 10, 0)]);
    const c = await getShapeProviderCenter(SyncRuntimeContext, provider);
    // within the three boxes' bounds, about their middle: not the first's at the origin
    expect(Vec3.magnitude(c)).toBeGreaterThan(5);
    expect(c[0]).toBeGreaterThan(5);
    expect(c[0]).toBeLessThan(15);
    expect(c[1]).toBeGreaterThan(1);
    expect(c[1]).toBeLessThan(10);
    expect(Math.abs(c[2])).toBeLessThan(1);
  });

  it('records the transform a provider moves by, composed over a nested one', () => {
    const provider = boxProvider(Vec3.create(1, 2, 3));
    expect(getShapeProviderTransform(provider)).toBeUndefined();
    const a = Mat4.fromTranslation(Mat4(), Vec3.create(1, 0, 0));
    const b = Mat4.fromScaling(Mat4(), Vec3.create(2, 2, 2));
    const moved = transformShapeProvider(provider, a);
    const twice = transformShapeProvider(moved, b);
    expect(Mat4.areEqual(getShapeProviderTransform(moved)!, a, 0)).toBe(true);
    expect(Mat4.areEqual(getShapeProviderTransform(twice)!, Mat4.mul(Mat4(), b, a), 0)).toBe(true);
  });

  it('tells providers moved alike from those moved otherwise, whatever their getters', () => {
    const provider = boxProvider(Vec3.create(1, 2, 3));
    // a provider's getter is built anew on every apply: not what tells them apart
    expect(isSameShapeProviderTransform(provider, boxProvider(Vec3.create(1, 2, 3)))).toBe(true);
    const move = Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0));
    const moved = transformShapeProvider(provider, move);
    expect(isSameShapeProviderTransform(provider, moved)).toBe(false);
    expect(isSameShapeProviderTransform(moved, transformShapeProvider(boxProvider(Vec3.unit), move))).toBe(true);
    const other = transformShapeProvider(provider, Mat4.fromTranslation(Mat4(), Vec3.create(0, 10, 0)));
    expect(isSameShapeProviderTransform(moved, other)).toBe(false);
  });
});
