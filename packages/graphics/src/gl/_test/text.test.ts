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
import { Text } from '@molstar/graphics/geo/geometry/text/text';

export function createText() {
    const text = Text.createEmpty();
    const props = PD.getDefaultValues(Text.Params);
    const values = Text.Utils.createValuesSimple(text, props, ColorNames.orange, 1);
    const state = Text.Utils.createRenderableState(props);
    return createRenderObject('text', values, state, -1);
}

describe('text', () => {
    const ctx = tryGetGLContext(32, 32);

    (ctx ? it : it.skip)('basic', async () => {
        const ctx = getGLContext(32, 32);
        const scene = Scene.create(ctx);
        const text = createText();
        scene.add(text);
        setDebugMode(true);
        expect(() => scene.commit()).not.toThrow();
        setDebugMode(false);
        ctx.destroy();
    });
});