/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createComputeRenderable, ComputeRenderable } from '../../renderable'
import { WebGLContext } from '../../webgl/context';
import { createComputeRenderItem } from '../../webgl/render-item';
import { Values, TextureSpec } from '../../renderable/schema';
import { Texture, createTexture } from 'mol-gl/webgl/texture';
import { ShaderCode } from 'mol-gl/shader-code';
import { ValueCell } from 'mol-util';
import { decodeFloatRGB } from 'mol-util/float-packing';
import { QuadSchema, QuadValues } from '../util';

const HistopyramidSumSchema = {
    ...QuadSchema,
    tTexture: TextureSpec('texture', 'rgba', 'float', 'nearest'),
}

let HistopyramidSumRenderable: ComputeRenderable<Values<typeof HistopyramidSumSchema>>
function getHistopyramidSumRenderable(ctx: WebGLContext, texture: Texture) {
    if (HistopyramidSumRenderable) {
        ValueCell.update(HistopyramidSumRenderable.values.tTexture, texture)
        HistopyramidSumRenderable.update()
        return HistopyramidSumRenderable
    } else {
        const values: Values<typeof HistopyramidSumSchema> = {
            ...QuadValues,
            tTexture: ValueCell.create(texture),
        }

        const schema = { ...HistopyramidSumSchema }
        const shaderCode = ShaderCode(
            require('mol-gl/shader/quad.vert').default,
            require('mol-gl/shader/histogram-pyramid/sum.frag').default
        )
        const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values)

        HistopyramidSumRenderable = createComputeRenderable(renderItem, values)
        return HistopyramidSumRenderable
    }
}

let SumTexture: Texture
function getSumTexture(ctx: WebGLContext) {
    if (SumTexture) return SumTexture
    SumTexture = createTexture(ctx, 'image-uint8', 'rgba', 'ubyte', 'nearest')
    SumTexture.define(1, 1)
    return SumTexture
}

/** name for shared framebuffer used for histogram-pyramid operations */
const FramebufferName = 'histogram-pyramid-sum'

const sumArray = new Uint8Array(4)
export function getHistopyramidSum(ctx: WebGLContext, pyramidTopTexture: Texture) {
    const { gl, framebufferCache } = ctx

    const renderable = getHistopyramidSumRenderable(ctx, pyramidTopTexture)
    const framebuffer = framebufferCache.get(FramebufferName).value
    const sumTexture = getSumTexture(ctx)
    sumTexture.attachFramebuffer(framebuffer, 0)

    gl.viewport(0, 0, 1, 1)
    renderable.render()
    ctx.readPixels(0, 0, 1, 1, sumArray)
    ctx.unbindFramebuffer()

    return decodeFloatRGB(sumArray[0], sumArray[1], sumArray[2])
}