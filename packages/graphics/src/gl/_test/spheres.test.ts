/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createRenderObject } from '../render-object.js';
import { Scene } from '../scene.js';
import { getGLContext, tryGetGLContext } from './gl.js';
import { setDebugMode } from '@molstar/core/util/debug';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';

export function createSpheres() {
  const spheres = Spheres.createEmpty();
  const props = PD.getDefaultValues(Spheres.Params);
  const values = Spheres.Utils.createValuesSimple(spheres, props, ColorNames.orange, 1);
  const state = Spheres.Utils.createRenderableState(props);
  return createRenderObject('spheres', values, state, -1);
}

describe('spheres', () => {
  // the radius of a single sphere's centers is 0 or a rounding error, by where it lies: the bounds are of its size
  it.each([
    [0, 0, 0],
    [10, 0, 0],
    [1, 2, 3],
  ])('bounds a single sphere by its size at %d, %d, %d', (x, y, z) => {
    const spheres = Spheres.create(new Float32Array([x, y, z]), new Float32Array(1), 1);
    const props = PD.getDefaultValues(Spheres.Params);
    const values = Spheres.Utils.createValuesSimple(spheres, props, ColorNames.orange, 2);
    expect(values.padding.ref.value).toBe(2);
    expect(values.boundingSphere.ref.value.radius).toBeCloseTo(2);
    expect(values.invariantBoundingSphere.ref.value.radius).toBeCloseTo(2);
  });

  it('bounds no spheres by nothing', () => {
    const props = PD.getDefaultValues(Spheres.Params);
    const values = Spheres.Utils.createValuesSimple(Spheres.createEmpty(), props, ColorNames.orange, 2);
    expect(values.padding.ref.value).toBe(0);
  });

  const ctx = tryGetGLContext(32, 32, { fragDepth: true });

  (ctx ? it : it.skip)('basic', { tags: ['medium'] }, async () => {
    const ctx = getGLContext(32, 32);
    const scene = Scene.create(ctx);
    const spheres = createSpheres();
    scene.add(spheres);
    setDebugMode(true);
    expect(() => scene.commit()).not.toThrow();
    setDebugMode(false);
    ctx.destroy();
  });
});
