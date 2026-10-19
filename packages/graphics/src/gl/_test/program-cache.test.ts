/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { createRenderObject } from '../render-object.js';
import { Scene } from '../scene.js';
import { getGLContext, tryGetGLContext } from './gl.js';

function createMesh() {
  const mesh = Mesh.createEmpty();
  const props = PD.getDefaultValues(Mesh.Params);
  const values = Mesh.Utils.createValuesSimple(mesh, props, ColorNames.orange, 1);
  const state = Mesh.Utils.createRenderableState(props);
  return createRenderObject('mesh', values, state, -1);
}

// the program cache keys the generated declarations by schema and defines (memoized): equal ones share programs,
// defines that change the declarations don't
describe('WebGL program cache', () => {
  const ctx = tryGetGLContext(32, 32);

  (ctx ? it : it.skip)('shares programs of equal defines, not of others', () => {
    const ctx = getGLContext(32, 32);
    const { resourceCounts } = ctx.stats;
    const scene = Scene.create(ctx);
    const programs = () => {
      scene.commit();
      ctx.resources.linkPrograms();
      return resourceCounts.program;
    };

    scene.add(createMesh());
    const first = programs();
    expect(first).toBeGreaterThan(0);

    scene.add(createMesh());
    expect(programs()).toBe(first);

    const other = createMesh();
    ValueCell.update(other.values.dFlatShaded, true);
    scene.add(other);
    expect(programs()).toBeGreaterThan(first);
    ctx.destroy();
  });
});
