/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderObject } from '../render-object';
import { Scene } from '../scene';
import { getGLContext, tryGetGLContext } from './gl';
import { setDebugMode } from '../../mol-util/debug';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Cylinders } from '../../mol-geo/geometry/cylinders/cylinders';

export function createCylinders() {
    const cylinders = Cylinders.createEmpty();
    const props = PD.getDefaultValues(Cylinders.Params);
    const values = Cylinders.Utils.createValuesSimple(cylinders, props, ColorNames.orange, 1);
    const state = Cylinders.Utils.createRenderableState(props);
    return createRenderObject('cylinders', values, state, -1);
}

describe('cylinders', () => {
    const ctx = tryGetGLContext(32, 32, { fragDepth: true });

    (ctx ? it : it.skip)('basic', async () => {
        const ctx = getGLContext(32, 32);
        const scene = Scene.create(ctx);
        const cylinders = createCylinders();
        scene.add(cylinders);
        setDebugMode(true);
        expect(() => scene.commit()).not.toThrow();
        setDebugMode(false);
        ctx.destroy();
    });
});