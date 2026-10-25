/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import {
  AutoGpuCullingMinInstances,
  AutoGpuCullingThreshold,
  resolveGpuCulling,
} from '@molstar/graphics/geo/geometry/base';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { SpheresBuilder } from '@molstar/graphics/geo/geometry/spheres/spheres-builder';
import { createTransform } from '@molstar/graphics/geo/geometry/transform-data';
import {
  createRenderObject,
  isMergedRenderObject,
  type GraphicsRenderObject,
} from '@molstar/graphics/gl/render-object';
import { mergeRenderObjectsByType } from '../merged.js';

/** A spheres render object of `instances` instances (translated copies of a few spheres) */
function spheres(instances: number): GraphicsRenderObject {
  const builder = SpheresBuilder.create(4, 2);
  for (let i = 0; i < 4; ++i) builder.add(i, 0, 0, i);
  const geometry = builder.getSpheres();
  const transforms = new Float32Array(instances * 16);
  for (let i = 0; i < instances; ++i)
    Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(0, i * 10, 0)), transforms, i * 16);
  const transform = createTransform(transforms, instances, Sphere3D.create(Vec3(), 5), 20, 64);
  const props = PD.getDefaultValues(Spheres.Params);
  const values = Spheres.Utils.createValuesSimple(geometry, props, ColorNames.orange, 1, transform);
  return createRenderObject('spheres', values, Spheres.Utils.createRenderableState(props), -1);
}

describe('merging render objects', () => {
  it('culls on the GPU from many instances, and many per draw (a segment and level)', () => {
    const min = AutoGpuCullingMinInstances;
    expect(resolveGpuCulling('auto', min, 4)).toBe(true);
    expect(resolveGpuCulling('auto', min - 1, 4)).toBe(false);
    const draws = Math.ceil(min / AutoGpuCullingThreshold) * 2;
    expect(resolveGpuCulling('auto', AutoGpuCullingThreshold * draws, draws)).toBe(true);
    expect(resolveGpuCulling('auto', AutoGpuCullingThreshold * draws - 1, draws)).toBe(false);
    expect(resolveGpuCulling(true, 1, 4)).toBe(true);
    expect(resolveGpuCulling(false, 1e6, 1)).toBe(false);
  });

  it('merges an object alone to cull it on the GPU, where the device does and it pays', () => {
    const entries = [{ visual: 0, renderObject: spheres(AutoGpuCullingMinInstances) }];
    const out: GraphicsRenderObject[] = [];
    const merged = mergeRenderObjectsByType(entries, 'auto', new Map(), out, 'auto', true);
    expect(out.length).toBe(1);
    expect(isMergedRenderObject(out[0])).toBe(true);
    expect((out[0].values as { gpuCulling?: { ref: { value: boolean } } }).gpuCulling?.ref.value).toBe(true);
    expect(merged.size).toBe(1);

    // without GPU culling on the device (WebGL) it stays as it is
    const plain: GraphicsRenderObject[] = [];
    mergeRenderObjectsByType(entries, 'auto', new Map(), plain, 'auto', false);
    expect(plain).toEqual([entries[0].renderObject]);
    // nor with fewer instances
    const few: GraphicsRenderObject[] = [];
    const fewer = spheres(AutoGpuCullingMinInstances - 1);
    mergeRenderObjectsByType([{ visual: 0, renderObject: fewer }], 'auto', new Map(), few, 'auto', true);
    expect(isMergedRenderObject(few[0])).toBe(false);
  });

  it('merges by the number of units as before where GPU culling does not pay', () => {
    const entries = Array.from({ length: 31 }, (_, i) => ({ visual: i, renderObject: spheres(1) }));
    const out: GraphicsRenderObject[] = [];
    mergeRenderObjectsByType(entries, 'auto', new Map(), out, 'auto', true);
    expect(out.length).toBe(1);
    expect((out[0].values as { gpuCulling?: { ref: { value: boolean } } }).gpuCulling?.ref.value).toBe(false);
    // fewer units than `auto` merges: unmerged
    const small: GraphicsRenderObject[] = [];
    mergeRenderObjectsByType(entries.slice(0, 5), 'auto', new Map(), small, 'auto', true);
    expect(small.length).toBe(5);
  });
});
