/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createComputeRenderable } from '../../renderable'
import { WebGLContext } from '../../webgl/context';
import { createComputeRenderItem } from '../../webgl/render-item';
import { AttributeSpec, Values, TextureSpec, ValueSpec } from '../../renderable/schema';
import { Texture, createTexture } from 'mol-gl/webgl/texture';
import { ShaderCode } from 'mol-gl/shader-code';
import { ValueCell } from 'mol-util';
import { decodeFloatRGB } from 'mol-util/float-packing';
import { QuadPositions, readTexture } from '../util';

const HistopyramidSumSchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),
    aPosition: AttributeSpec('float32', 2, 0),

    tTexture: TextureSpec('texture', 'rgba', 'float', 'nearest'),
}

function getHistopyramidSumRenderable(ctx: WebGLContext, texture: Texture) {
    const values: Values<typeof HistopyramidSumSchema> = {
        drawCount: ValueCell.create(6),
        instanceCount: ValueCell.create(1),
        aPosition: ValueCell.create(QuadPositions),

        tTexture: ValueCell.create(texture),
    }

    const schema = { ...HistopyramidSumSchema }
    const shaderCode = ShaderCode(
        require('mol-gl/shader/quad.vert').default,
        require('mol-gl/shader/histogram-pyramid/sum.frag').default,
        { standardDerivatives: false, fragDepth: false }
    )
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values)

    return createComputeRenderable(renderItem, values);
}

/** name for shared framebuffer used for histogram-pyramid operations */
const FramebufferName = 'histogram-pyramid'

export function getHistopyramidSum(ctx: WebGLContext, pyramidTopTexture: Texture) {
    const { gl, framebufferCache } = ctx

    const framebuffer = framebufferCache.get(FramebufferName).value

    const encodeFloatRenderable = getHistopyramidSumRenderable(ctx, pyramidTopTexture)
    encodeFloatRenderable.update()
    encodeFloatRenderable.use()

    const encodedFloatTexture = createTexture(ctx, 'image-uint8', 'rgba', 'ubyte', 'nearest')
    encodedFloatTexture.define(1, 1)
    encodedFloatTexture.attachFramebuffer(framebuffer, 0)

    gl.viewport(0, 0, 1, 1)
    encodeFloatRenderable.render()

    const sumImage = readTexture(ctx, encodedFloatTexture)
    return decodeFloatRGB(sumImage.array[0], sumImage.array[1], sumImage.array[2])
}