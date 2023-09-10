/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { Vec2 } from '../../mol-math/linear-algebra/3d/vec2';
import { ValueCell } from '../../mol-util/value-cell';
import { ComputeRenderable, createComputeRenderable } from '../renderable';
import { TextureSpec, UniformSpec, Values } from '../renderable/schema';
import { ShaderCode } from '../shader-code';
import { hiZ_frag } from '../shader/hi-z.frag';
import { quad_vert } from '../shader/quad.vert';
import { createComputeRenderItem } from '../webgl/render-item';
import { Texture } from '../webgl/texture';
import { QuadSchema, QuadValues } from './util';


const HiZSchema = {
    ...QuadSchema,
    tPreviousLevel: TextureSpec('texture', 'alpha', 'float', 'nearest'),
    uInvSize: UniformSpec('v2'),
    uOffset: UniformSpec('v2'),
};
const HiZShaderCode = ShaderCode('hi-z', quad_vert, hiZ_frag);
export type HiZRenderable = ComputeRenderable<Values<typeof HiZSchema>>

export function createHiZRenderable(ctx: WebGLContext, previousLevel: Texture): HiZRenderable {
    const values: Values<typeof HiZSchema> = {
        ...QuadValues,
        tPreviousLevel: ValueCell.create(previousLevel),
        uInvSize: ValueCell.create(Vec2()),
        uOffset: ValueCell.create(Vec2()),
    };

    const schema = { ...HiZSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', HiZShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}