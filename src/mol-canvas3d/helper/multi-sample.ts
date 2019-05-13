/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from 'mol-gl/compute/util';
import { TextureSpec, UniformSpec, Values } from 'mol-gl/renderable/schema';
import { Texture } from 'mol-gl/webgl/texture';
import { WebGLContext } from 'mol-gl/webgl/context';
import { ValueCell } from 'mol-util';
import { Vec2 } from 'mol-math/linear-algebra';
import { ShaderCode } from 'mol-gl/shader-code';
import { createComputeRenderItem } from 'mol-gl/webgl/render-item';
import { createComputeRenderable } from 'mol-gl/renderable';

const ComposeSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),
    uWeight: UniformSpec('f'),
}

export function getComposeRenderable(ctx: WebGLContext, colorTexture: Texture) {
    const values: Values<typeof ComposeSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.width, colorTexture.height)),
        uWeight: ValueCell.create(1.0),
    }

    const schema = { ...ComposeSchema }
    const shaderCode = ShaderCode(
        require('mol-gl/shader/quad.vert').default,
        require('mol-gl/shader/compose.frag').default
    )
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values)

    return createComputeRenderable(renderItem, values)
}

export const JitterVectors = [
    [
        [ 0, 0 ]
    ],
    [
        [ 4, 4 ], [ -4, -4 ]
    ],
    [
        [ -2, -6 ], [ 6, -2 ], [ -6, 2 ], [ 2, 6 ]
    ],
    [
        [ 1, -3 ], [ -1, 3 ], [ 5, 1 ], [ -3, -5 ],
        [ -5, 5 ], [ -7, -1 ], [ 3, 7 ], [ 7, -7 ]
    ],
    [
        [ 1, 1 ], [ -1, -3 ], [ -3, 2 ], [ 4, -1 ],
        [ -5, -2 ], [ 2, 5 ], [ 5, 3 ], [ 3, -5 ],
        [ -2, 6 ], [ 0, -7 ], [ -4, -6 ], [ -6, 4 ],
        [ -8, 0 ], [ 7, -4 ], [ 6, 7 ], [ -7, -8 ]
    ],
    [
        [ -4, -7 ], [ -7, -5 ], [ -3, -5 ], [ -5, -4 ],
        [ -1, -4 ], [ -2, -2 ], [ -6, -1 ], [ -4, 0 ],
        [ -7, 1 ], [ -1, 2 ], [ -6, 3 ], [ -3, 3 ],
        [ -7, 6 ], [ -3, 6 ], [ -5, 7 ], [ -1, 7 ],
        [ 5, -7 ], [ 1, -6 ], [ 6, -5 ], [ 4, -4 ],
        [ 2, -3 ], [ 7, -2 ], [ 1, -1 ], [ 4, -1 ],
        [ 2, 1 ], [ 6, 2 ], [ 0, 4 ], [ 4, 4 ],
        [ 2, 5 ], [ 7, 5 ], [ 5, 6 ], [ 3, 7 ]
    ]
]
  
JitterVectors.forEach(offsetList => {
    offsetList.forEach(offset => {
        // 0.0625 = 1 / 16
        offset[0] *= 0.0625
        offset[1] *= 0.0625
    })
})