/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from 'mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec } from 'mol-gl/renderable/schema';
import { ShaderCode } from 'mol-gl/shader-code';
import { WebGLContext } from 'mol-gl/webgl/context';
import { Texture } from 'mol-gl/webgl/texture';
import { ValueCell } from 'mol-util';
import { createComputeRenderItem } from 'mol-gl/webgl/render-item';
import { createComputeRenderable } from 'mol-gl/renderable';
import { Vec2 } from 'mol-math/linear-algebra';
import { ParamDefinition as PD } from 'mol-util/param-definition';

const SSAOPassSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    uEnable: UniformSpec('i'),
    uKernelSize: UniformSpec('i'),
    uBias: UniformSpec('f'),
    uRadius: UniformSpec('f'),
}

export const SSAOPassParams = {
    enable: PD.Boolean(true),
    kernelSize: PD.Numeric(4, { min: 1, max: 100, step: 1 }),
    bias: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
    radius: PD.Numeric(128, { min: 0, max: 256, step: 1 }),
}
export type SSAOPassProps = PD.Values<typeof SSAOPassParams>

export function getSSAOPassRenderable(ctx: WebGLContext, colorTexture: Texture, depthTexture: Texture, props: Partial<SSAOPassProps>) {
    const p = { ...PD.getDefaultValues(SSAOPassParams), props }
    const values: Values<typeof SSAOPassSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        tDepth: ValueCell.create(depthTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.width, colorTexture.height)),

        uEnable: ValueCell.create(p.enable ? 1 : 0),
        uKernelSize: ValueCell.create(p.kernelSize),
        uBias: ValueCell.create(p.bias),
        uRadius: ValueCell.create(p.radius),
    }

    const schema = { ...SSAOPassSchema }
    const shaderCode = ShaderCode(
        require('mol-gl/shader/quad.vert').default,
        require('mol-gl/shader/passes/ssao.frag').default
    )
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values)

    return createComputeRenderable(renderItem, values)
}