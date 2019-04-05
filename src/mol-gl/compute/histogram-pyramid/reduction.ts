/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createComputeRenderable } from '../../renderable'
import { WebGLContext } from '../../webgl/context';
import { createComputeRenderItem } from '../../webgl/render-item';
import { Values, TextureSpec, UniformSpec } from '../../renderable/schema';
import { Texture, createTexture } from 'mol-gl/webgl/texture';
import { ShaderCode } from 'mol-gl/shader-code';
import { ValueCell } from 'mol-util';
import { GLRenderingContext } from 'mol-gl/webgl/compat';
import { QuadSchema, QuadValues } from '../util';
import { Vec2 } from 'mol-math/linear-algebra';
import { getHistopyramidSum } from './sum';

const HistopyramidReductionSchema = {
    ...QuadSchema,
    tPreviousLevel: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    uSize: UniformSpec('f'),
}

function getHistopyramidReductionRenderable(ctx: WebGLContext, initialTexture: Texture) {
    const values: Values<typeof HistopyramidReductionSchema> = {
        ...QuadValues,
        tPreviousLevel: ValueCell.create(initialTexture),
        uSize: ValueCell.create(0),
    }

    const schema = { ...HistopyramidReductionSchema }
    const shaderCode = ShaderCode(
        require('mol-gl/shader/quad.vert').default,
        require('mol-gl/shader/histogram-pyramid/reduction.frag').default
    )
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values)

    return createComputeRenderable(renderItem, values);
}

/** name for shared framebuffer used for histogram-pyramid operations */
const FramebufferName = 'histogram-pyramid-reduction'

function setRenderingDefaults(gl: GLRenderingContext) {
    gl.disable(gl.CULL_FACE)
    gl.disable(gl.BLEND)
    gl.disable(gl.DEPTH_TEST)
    gl.depthMask(false)
}

export interface HistogramPyramid {
    pyramidTex: Texture
    totalTex: Texture
    initialTex: Texture
    count: number
    height: number
    levels: number
    scale: Vec2
}

export function createHistogramPyramid(ctx: WebGLContext, inputTexture: Texture): HistogramPyramid {
    const { gl, framebufferCache } = ctx

    const inputTextureMaxDim = Math.max(inputTexture.width, inputTexture.height)

    // This part set the levels
    const levels = Math.ceil(Math.log(inputTextureMaxDim) / Math.log(2))
    // console.log('levels', levels)

    const initialTexture = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    initialTexture.load({ array: new Float32Array(4), width: 1, height: 1 })
    initialTexture.define(Math.pow(2, levels), Math.pow(2, levels))

    const framebuffer = framebufferCache.get(FramebufferName).value
    inputTexture.attachFramebuffer(framebuffer, 0)
    initialTexture.define(Math.pow(2, levels), Math.pow(2, levels))
    // TODO need to initialize texSubImage2D to make Firefox happy
    gl.copyTexSubImage2D(gl.TEXTURE_2D, 0, 0, 0, 0, 0, inputTexture.width, inputTexture.height);

    const initialTextureMaxDim = Math.max(initialTexture.width, initialTexture.height)

    const pyramidTexture = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    pyramidTexture.define(Math.pow(2, levels), Math.pow(2, levels))

    // TODO cache globally for reuse
    const levelTextures: Texture[] = []
    for (let i = 0; i < levels; ++i) {
        const tex = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
        tex.define(Math.pow(2, i), Math.pow(2, i))
        levelTextures.push(tex)
    }

    const renderable = getHistopyramidReductionRenderable(ctx, initialTexture)
    renderable.update()
    renderable.use()

    let offset = 0;
    for (let i = 0; i < levels; i++) {
        const currLevel = levels - 1 - i
        levelTextures[currLevel].attachFramebuffer(framebuffer, 0)

        const size = Math.pow(2, currLevel)
        // console.log('size', size, 'draw-level', currLevel, 'read-level', levels - i)
        gl.clear(gl.COLOR_BUFFER_BIT)

        ValueCell.update(renderable.values.uSize, Math.pow(2, i + 1) / initialTextureMaxDim)
        const readTex = i === 0 ? initialTexture : levelTextures[levels - i]
        // console.log(readTex.width, readTex.height)
        ValueCell.update(renderable.values.tPreviousLevel, readTex)

        renderable.update()
        renderable.use()
        setRenderingDefaults(gl)
        gl.viewport(0, 0, size, size)
        renderable.render()

        pyramidTexture.bind(0)
        // TODO need to initialize texSubImage2D to make Firefox happy
        gl.copyTexSubImage2D(gl.TEXTURE_2D, 0, offset, 0, 0, 0, size, size);
        pyramidTexture.unbind(0)

        // if (i >= levels - 4) {
        //     console.log('==============', i)
        //     const rt = readTexture(ctx, levelTextures[currLevel])
        //     console.log('array', rt.array)
        //     for (let i = 0, il = rt.width * rt.height; i < il; ++i) {
        //         // const v = decodeFloatRGB(rt.array[i * 4], rt.array[i * 4 + 1], rt.array[i * 4 + 2])
        //         // console.log(i, 'v', v, 'rgb', rt.array[i * 4], rt.array[i * 4 + 1], rt.array[i * 4 + 2])
        //         console.log(i, 'f', rt.array[i * 4])
        //     }
        // }

        offset += size;
    }

    // printTexture(ctx, pyramidTexture, 3)

    //

    const finalCount = getHistopyramidSum(ctx, levelTextures[0])
    const height = Math.ceil(finalCount / Math.pow(2, levels))
    // console.log('height', height, 'finalCount', finalCount)

    //

    const scale = Vec2.create(
        initialTexture.width / inputTexture.width,
        initialTexture.height / inputTexture.height
    )
    return {
        pyramidTex: pyramidTexture,
        totalTex: levelTextures[0],
        initialTex: initialTexture,
        count: finalCount,
        height,
        levels,
        scale
    }
}