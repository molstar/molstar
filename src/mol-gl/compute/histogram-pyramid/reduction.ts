/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createComputeRenderable, ComputeRenderable } from '../../renderable'
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

let HistopyramidReductionRenderable: ComputeRenderable<Values<typeof HistopyramidReductionSchema>>
function getHistopyramidReductionRenderable(ctx: WebGLContext, initialTexture: Texture) {
    if (HistopyramidReductionRenderable) {
        ValueCell.update(HistopyramidReductionRenderable.values.tPreviousLevel, initialTexture)
        HistopyramidReductionRenderable.update()
        return HistopyramidReductionRenderable
    } else {
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

        HistopyramidReductionRenderable = createComputeRenderable(renderItem, values);
        return HistopyramidReductionRenderable
    }
}

/** name for shared framebuffer used for histogram-pyramid operations */
const FramebufferName = 'histogram-pyramid-reduction'

const LevelTextures: Texture[] = []
function getLevelTexture(ctx: WebGLContext, level: number) {
    let tex = LevelTextures[level]
    const size = Math.pow(2, level)
    if (tex === undefined) {
        tex = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
        LevelTextures[level] = tex
    }
    tex.define(size, size) // always call to set size AND clear
    return tex
}

function setRenderingDefaults(gl: GLRenderingContext) {
    gl.disable(gl.CULL_FACE)
    gl.disable(gl.BLEND)
    gl.disable(gl.DEPTH_TEST)
    gl.depthMask(false)
}

export interface HistogramPyramid {
    pyramidTex: Texture
    count: number
    height: number
    levels: number
    scale: Vec2
}

export function createHistogramPyramid(ctx: WebGLContext, inputTexture: Texture): HistogramPyramid {
    const { gl, framebufferCache } = ctx

    // printTexture(ctx, inputTexture, 2)
    const inputTextureMaxDim = Math.max(inputTexture.width, inputTexture.height)

    // This part set the levels
    const levels = Math.ceil(Math.log(inputTextureMaxDim) / Math.log(2))
    const maxSize = Math.pow(2, levels)
    // console.log('levels', levels, 'maxSize', maxSize)

    const initialTexture = getLevelTexture(ctx, levels)

    const framebuffer = framebufferCache.get(FramebufferName).value
    inputTexture.attachFramebuffer(framebuffer, 0)
    // TODO need to initialize texSubImage2D to make Firefox happy
    gl.copyTexSubImage2D(gl.TEXTURE_2D, 0, 0, 0, 0, 0, inputTexture.width, inputTexture.height);

    const pyramidTexture = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    pyramidTexture.define(maxSize, maxSize)

    const levelTextures: Texture[] = []
    for (let i = 0; i < levels; ++i) levelTextures.push(getLevelTexture(ctx, i))

    const renderable = getHistopyramidReductionRenderable(ctx, initialTexture)
    setRenderingDefaults(gl)

    let offset = 0;
    for (let i = 0; i < levels; i++) {
        const currLevel = levels - 1 - i
        levelTextures[currLevel].attachFramebuffer(framebuffer, 0)

        const size = Math.pow(2, currLevel)
        // console.log('size', size, 'draw-level', currLevel, 'read-level', levels - i)
        gl.clear(gl.COLOR_BUFFER_BIT)
        gl.viewport(0, 0, size, size)

        ValueCell.update(renderable.values.uSize, Math.pow(2, i + 1) / maxSize)
        if (i > 0) {
            ValueCell.update(renderable.values.tPreviousLevel, levelTextures[levels - i])
            renderable.update()
        }
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

    // printTexture(ctx, pyramidTexture, 2)

    //

    const finalCount = getHistopyramidSum(ctx, levelTextures[0])
    const height = Math.ceil(finalCount / Math.pow(2, levels))
    const scale = Vec2.create(maxSize / inputTexture.width, maxSize / inputTexture.height)
    // console.log('height', height, 'finalCount', finalCount, 'scale', scale)

    return {
        pyramidTex: pyramidTexture,
        count: finalCount,
        height,
        levels,
        scale
    }
}