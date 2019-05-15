/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from 'mol-gl/compute/util';
import { TextureSpec, Values, UniformSpec, DefineSpec } from 'mol-gl/renderable/schema';
import { ShaderCode } from 'mol-gl/shader-code';
import { WebGLContext } from 'mol-gl/webgl/context';
import { Texture } from 'mol-gl/webgl/texture';
import { ValueCell } from 'mol-util';
import { createComputeRenderItem } from 'mol-gl/webgl/render-item';
import { createComputeRenderable, ComputeRenderable } from 'mol-gl/renderable';
import { Vec2 } from 'mol-math/linear-algebra';
import { ParamDefinition as PD } from 'mol-util/param-definition';

const PostprocessingSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    dOcclusionEnable: DefineSpec('boolean'),
    dOcclusionKernelSize: DefineSpec('number'),
    uOcclusionBias: UniformSpec('f'),
    uOcclusionRadius: UniformSpec('f'),

    dOutlineEnable: DefineSpec('boolean'),
    uOutlineScale: UniformSpec('f'),
    uOutlineThreshold: UniformSpec('f'),
}

export const PostprocessingParams = {
    occlusionEnable: PD.Boolean(false),
    occlusionKernelSize: PD.Numeric(4, { min: 1, max: 32, step: 1 }),
    occlusionBias: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
    occlusionRadius: PD.Numeric(64, { min: 0, max: 256, step: 1 }),

    outlineEnable: PD.Boolean(false),
    outlineScale: PD.Numeric(1, { min: 0, max: 10, step: 1 }),
    outlineThreshold: PD.Numeric(0.8, { min: 0, max: 1, step: 0.01 }),
}
export type PostprocessingProps = PD.Values<typeof PostprocessingParams>

type PostprocessingRenderable = ComputeRenderable<Values<typeof PostprocessingSchema>>

export function getPostprocessingRenderable(ctx: WebGLContext, colorTexture: Texture, depthTexture: Texture, props: Partial<PostprocessingProps>): PostprocessingRenderable {
    const p = { ...PD.getDefaultValues(PostprocessingParams), props }
    const values: Values<typeof PostprocessingSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        tDepth: ValueCell.create(depthTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.width, colorTexture.height)),

        dOcclusionEnable: ValueCell.create(p.occlusionEnable),
        dOcclusionKernelSize: ValueCell.create(p.occlusionKernelSize),
        uOcclusionBias: ValueCell.create(p.occlusionBias),
        uOcclusionRadius: ValueCell.create(p.occlusionRadius),

        dOutlineEnable: ValueCell.create(p.outlineEnable),
        uOutlineScale: ValueCell.create(p.outlineScale * ctx.pixelRatio),
        uOutlineThreshold: ValueCell.create(p.outlineThreshold),
    }

    const schema = { ...PostprocessingSchema }
    const shaderCode = ShaderCode(
        require('mol-gl/shader/quad.vert').default,
        require('mol-gl/shader/postprocessing.frag').default
    )
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values)

    return createComputeRenderable(renderItem, values)
}

export function setPostprocessingProps(props: Partial<PostprocessingProps>, postprocessing: PostprocessingRenderable, currentProps: PostprocessingProps, webgl: WebGLContext) {
    if (props.occlusionEnable !== undefined) {
        currentProps.occlusionEnable = props.occlusionEnable
        ValueCell.update(postprocessing.values.dOcclusionEnable, props.occlusionEnable)
    }
    if (props.occlusionKernelSize !== undefined) {
        currentProps.occlusionKernelSize = props.occlusionKernelSize
        ValueCell.update(postprocessing.values.dOcclusionKernelSize, props.occlusionKernelSize)
    }
    if (props.occlusionBias !== undefined) {
        currentProps.occlusionBias = props.occlusionBias
        ValueCell.update(postprocessing.values.uOcclusionBias, props.occlusionBias)
    }
    if (props.occlusionRadius !== undefined) {
        currentProps.occlusionRadius = props.occlusionRadius
        ValueCell.update(postprocessing.values.uOcclusionRadius, props.occlusionRadius)
    }

    if (props.outlineEnable !== undefined) {
        currentProps.outlineEnable = props.outlineEnable
        ValueCell.update(postprocessing.values.dOutlineEnable, props.outlineEnable)
    }
    if (props.outlineScale !== undefined) {
        currentProps.outlineScale = props.outlineScale
        ValueCell.update(postprocessing.values.uOutlineScale, props.outlineScale * webgl.pixelRatio)
    }
    if (props.outlineThreshold !== undefined) {
        currentProps.outlineThreshold = props.outlineThreshold
        ValueCell.update(postprocessing.values.uOutlineThreshold, props.outlineThreshold)
    }

    postprocessing.update()
}