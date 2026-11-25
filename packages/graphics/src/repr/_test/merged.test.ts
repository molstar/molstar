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
  AutoGpuCullingMinRasterSpheres,
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
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { mergeRenderObjectsByType, resolveMergedCulling } from '../merged.js';

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

  it('culls on the GPU from many spheres the raster could draw', () => {
    const min = AutoGpuCullingMinRasterSpheres;
    expect(resolveGpuCulling('auto', 100, 100, min)).toBe(true);
    expect(resolveGpuCulling('auto', 100, 100, min - 1)).toBe(false);
    expect(resolveGpuCulling(false, 100, 100, min)).toBe(false);
  });

  it('culls spheres merged anyway for the raster at any count, where the device has it', () => {
    const graphics = { capabilities: { atomicMinMax64: true } } as GraphicsContext;
    const records = (n: number) => Array.from({ length: n }, () => spheres(2).values as never);
    // a units visual's one record (merged whatever), and two objects merged: for the raster only
    expect(resolveMergedCulling('spheres', records(1), 'auto', true, graphics, true)).toEqual({
      merge: true,
      gpu: true,
      rasterOnly: true,
    });
    expect(resolveMergedCulling('spheres', records(2), 'auto', true, graphics)).toMatchObject({
      gpu: true,
      rasterOnly: true,
    });
    // one object not merged otherwise: not merged for it
    expect(resolveMergedCulling('spheres', records(1), 'auto', true, graphics)).toMatchObject({
      merge: false,
      gpu: false,
    });
    // nor without the raster, without GPU culling, with it off, or of another type
    const plain = { capabilities: { atomicMinMax64: false } } as GraphicsContext;
    expect(resolveMergedCulling('spheres', records(2), 'auto', true, plain).gpu).toBe(false);
    expect(resolveMergedCulling('spheres', records(2), 'auto', false, graphics).gpu).toBe(false);
    expect(resolveMergedCulling('spheres', records(2), false, true, graphics).gpu).toBe(false);
    expect(resolveMergedCulling('mesh', records(2), 'auto', true, graphics).gpu).toBe(false);
  });

  it('merges an object alone to cull it on the GPU, where the device does and it pays', () => {
    const entries = [{ visual: 0, renderObject: spheres(AutoGpuCullingMinInstances) }];
    const out: GraphicsRenderObject[] = [];
    const merged = mergeRenderObjectsByType(entries, new Map(), out, 'auto', true);
    expect(out.length).toBe(1);
    expect(isMergedRenderObject(out[0])).toBe(true);
    expect((out[0].values as { gpuCulling?: { ref: { value: boolean } } }).gpuCulling?.ref.value).toBe(true);
    expect(merged.size).toBe(1);

    // without GPU culling on the device (WebGL) it stays as it is
    const plain: GraphicsRenderObject[] = [];
    mergeRenderObjectsByType(entries, new Map(), plain, 'auto', false);
    expect(plain).toEqual([entries[0].renderObject]);
    // nor with fewer instances
    const few: GraphicsRenderObject[] = [];
    const fewer = spheres(AutoGpuCullingMinInstances - 1);
    mergeRenderObjectsByType([{ visual: 0, renderObject: fewer }], new Map(), few, 'auto', true);
    expect(isMergedRenderObject(few[0])).toBe(false);
  });

  it('marks an object culled on the GPU for the raster only, by its spheres and not its instances', () => {
    const graphics = { capabilities: { atomicMinMax64: true } } as GraphicsContext;
    const rasterOnly = (o: GraphicsRenderObject) =>
      (o.values as { gpuCullingRasterOnly?: { ref: { value: boolean } } }).gpuCullingRasterOnly?.ref.value;
    const instances = 30_000;
    const ro = spheres(instances);
    expect(instances * ro.values.uVertexCount.ref.value).toBeGreaterThanOrEqual(AutoGpuCullingMinRasterSpheres);
    const out: GraphicsRenderObject[] = [];
    mergeRenderObjectsByType([{ visual: 0, renderObject: ro }], new Map(), out, 'auto', true, graphics);
    expect((out[0].values as { gpuCulling?: { ref: { value: boolean } } }).gpuCulling?.ref.value).toBe(true);
    expect(rasterOnly(out[0])).toBe(true);
    // by its instances: culled on the GPU whatever the raster takes
    const many: GraphicsRenderObject[] = [];
    const entries = [{ visual: 0, renderObject: spheres(AutoGpuCullingMinInstances) }];
    mergeRenderObjectsByType(entries, new Map(), many, 'auto', true, graphics);
    expect(rasterOnly(many[0])).toBe(false);
    // asked for: not for the raster only
    const asked: GraphicsRenderObject[] = [];
    mergeRenderObjectsByType([{ visual: 0, renderObject: ro }], new Map(), asked, true, true, graphics);
    expect(rasterOnly(asked[0])).toBe(false);
  });

  it('merges two or more where GPU culling does not pay, one alone not', () => {
    const entries = Array.from({ length: 5 }, (_, i) => ({ visual: i, renderObject: spheres(1) }));
    const out: GraphicsRenderObject[] = [];
    mergeRenderObjectsByType(entries, new Map(), out, 'auto', true);
    expect(out.length).toBe(1);
    expect((out[0].values as { gpuCulling?: { ref: { value: boolean } } }).gpuCulling?.ref.value).toBe(false);
    const one: GraphicsRenderObject[] = [];
    mergeRenderObjectsByType(entries.slice(0, 1), new Map(), one, 'auto', true);
    expect(one).toEqual([entries[0].renderObject]);
  });
});
