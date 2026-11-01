/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util/value-cell';
import { createTransform } from '@molstar/graphics/geo/geometry/transform-data';
import type { GraphicsRenderObject } from '@molstar/graphics/gl/render-object';
import { Visual } from '../visual.js';

/** The values `Visual.setTransform` reads: the transform data of `count` identity instances and the spheres */
function renderObject(count: number) {
  const transforms = new Float32Array(count * 16);
  for (let i = 0; i < count; ++i) Mat4.toArray(Mat4.identity(), transforms, i * 16);
  const values = {
    ...createTransform(transforms, count, undefined, 0, 0),
    invariantBoundingSphere: ValueCell.create(Sphere3D.create(Vec3.zero(), 1)),
    boundingSphere: ValueCell.create(Sphere3D()),
  };
  return { values } as unknown as GraphicsRenderObject;
}

describe('Visual.setTransform', () => {
  it('allocates the extra transforms when they are first set', () => {
    const ro = renderObject(2);
    expect(ro.values.extraTransform.ref.value.length).toBe(0);

    const extra = new Float32Array(32);
    Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(1, 0, 0)), extra, 0);
    Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(0, 2, 0)), extra, 16);
    Visual.setTransform(ro, undefined, extra);

    expect(ro.values.hasExtraTransform.ref.value).toBe(true);
    expect(Array.from(ro.values.extraTransform.ref.value)).toEqual(Array.from(extra));
    // in the instance grid's order: `aInstance` names the instance of each
    const { aTransform, aInstance } = ro.values;
    const translation = (i: number) => {
      const p = Array.from(aInstance.ref.value).indexOf(i);
      return Array.from(aTransform.ref.value.subarray(p * 16 + 12, p * 16 + 15));
    };
    expect(translation(0)).toEqual([1, 0, 0]);
    expect(translation(1)).toEqual([0, 2, 0]);
  });
});
