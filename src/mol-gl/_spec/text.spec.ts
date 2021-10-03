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
import { Text } from '../../mol-geo/geometry/text/text';

export function createText() {
    const text = Text.createEmpty();
    const props = PD.getDefaultValues(Text.Params);
    const values = Text.Utils.createValuesSimple(text, props, ColorNames.orange, 1);
    const state = Text.Utils.createRenderableState(props);
    return createRenderObject('text', values, state, -1);
}

describe('text', () => {
    it('basic', async () => {
        const gl = getGLContext(32, 32);
        const { ctx } = createRenderer(gl);
        const scene = Scene.create(ctx);
        const text = createText();
        scene.add(text);
        setDebugMode(true);
        expect(() => scene.commit()).not.toThrow();
        setDebugMode(false);
        gl.getExtension('STACKGL_destroy_context')?.destroy();
    });
});