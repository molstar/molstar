/**
 * Copyright (c) 2021-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderObject } from '../render-object';
import { Scene } from '../scene';
import { getGLContext, tryGetGLContext } from './gl';
import { setDebugMode } from '../../mol-util/debug';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Spheres } from '../../mol-geo/geometry/spheres/spheres';

export function createSpheres() {
    const spheres = Spheres.createEmpty();
    const props = PD.getDefaultValues(Spheres.Params);
    const values = Spheres.Utils.createValuesSimple(spheres, props, ColorNames.orange, 1);
    const state = Spheres.Utils.createRenderableState(props);
    return createRenderObject('spheres', values, state, -1);
}

describe('spheres', () => {
    const ctx = tryGetGLContext(32, 32, { fragDepth: true, textureFloat: true });

    (ctx ? it : it.skip)('basic', async () => {
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