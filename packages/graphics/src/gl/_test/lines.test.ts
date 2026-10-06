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
import { Lines } from '@molstar/graphics/geo/geometry/lines/lines';

export function createLines() {
  const lines = Lines.createEmpty();
  const props = PD.getDefaultValues(Lines.Params);
  const values = Lines.Utils.createValuesSimple(lines, props, ColorNames.orange, 1);
  const state = Lines.Utils.createRenderableState(props);
  return createRenderObject('lines', values, state, -1);
}

describe('lines', () => {
  const ctx = tryGetGLContext(32, 32);

  (ctx ? it : it.skip)('basic', async () => {
    const ctx = getGLContext(32, 32);
    const scene = Scene.create(ctx);
    const lines = createLines();
    scene.add(lines);
    setDebugMode(true);
    expect(() => scene.commit()).not.toThrow();
    setDebugMode(false);
    ctx.destroy();
  });
});
