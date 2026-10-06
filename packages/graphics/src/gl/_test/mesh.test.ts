/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderObject } from '../render-object.js';
import { Scene } from '../scene.js';
import { getGLContext, tryGetGLContext } from './gl.js';
import { setDebugMode } from '@molstar/core/util/debug';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';

export function createMesh() {
  const mesh = Mesh.createEmpty();
  const props = PD.getDefaultValues(Mesh.Params);
  const values = Mesh.Utils.createValuesSimple(mesh, props, ColorNames.orange, 1);
  const state = Mesh.Utils.createRenderableState(props);
  return createRenderObject('mesh', values, state, -1);
}

describe('mesh', () => {
  const ctx = tryGetGLContext(32, 32);

  (ctx ? it : it.skip)('basic', async () => {
    const ctx = getGLContext(32, 32);
    const scene = Scene.create(ctx);
    const mesh = createMesh();
    scene.add(mesh);
    setDebugMode(true);
    expect(() => scene.commit()).not.toThrow();
    setDebugMode(false);
    ctx.destroy();
  });
});
