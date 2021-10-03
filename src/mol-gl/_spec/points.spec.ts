/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createRenderObject } from '../render-object';
import { Scene } from '../scene';
import getGLContext from 'gl';
import { setDebugMode } from '../../mol-util/debug';
import { createRenderer } from './renderer.spec';
import { ColorNames } from '../../mol-util/color/names';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Points } from '../../mol-geo/geometry/points/points';

export function createPoints() {
    const points = Points.createEmpty();
    const props = PD.getDefaultValues(Points.Params);
    const values = Points.Utils.createValuesSimple(points, props, ColorNames.orange, 1);
    const state = Points.Utils.createRenderableState(props);
    return createRenderObject('points', values, state, -1);
}

describe('points', () => {
    it('basic', async () => {
        const gl = getGLContext(32, 32);
        const { ctx } = createRenderer(gl);
        const scene = Scene.create(ctx);
        const points = createPoints();
        scene.add(points);
        setDebugMode(true);
        expect(() => scene.commit()).not.toThrow();
        setDebugMode(false);
        gl.getExtension('STACKGL_destroy_context')?.destroy();
    });
});