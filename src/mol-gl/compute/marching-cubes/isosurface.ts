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
import { Vec3, Vec2, Mat4 } from 'mol-math/linear-algebra';
import { QuadSchema, QuadValues } from '../util';
import { HistogramPyramid } from '../histogram-pyramid/reduction';
import { getTriIndices } from './tables';

/** name for shared framebuffer used for gpu marching cubes operations */
const FramebufferName = 'marching-cubes-isosurface'

const IsosurfaceSchema = {
    ...QuadSchema,

    tTriIndices: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    tActiveVoxelsPyramid: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    tActiveVoxelsBase: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    tVolumeData: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tActiveVoxelsTotal: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    uIsoValue: UniformSpec('f'),

    uSize: UniformSpec('f'),
    uLevels: UniformSpec('f'),

    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v3'),
    uGridTransform: UniformSpec('m4'),

    uScale: UniformSpec('v2'),
}

function getIsosurfaceRenderable(ctx: WebGLContext, activeVoxelsPyramid: Texture, activeVoxelsBase: Texture, volumeData: Texture, activeVoxelsTotal: Texture, gridDimensions: Vec3, transform: Mat4, isoValue: number, levels: number, scale: Vec2) {
    // console.log('uSize', Math.pow(2, levels))
    const values: Values<typeof IsosurfaceSchema> = {
        ...QuadValues,

        tTriIndices: ValueCell.create(getTriIndices()),
        tActiveVoxelsPyramid: ValueCell.create(activeVoxelsPyramid),
        tActiveVoxelsBase: ValueCell.create(activeVoxelsBase),
        tVolumeData: ValueCell.create(volumeData),
        tActiveVoxelsTotal: ValueCell.create(activeVoxelsTotal),
        uIsoValue: ValueCell.create(isoValue),

        uSize: ValueCell.create(Math.pow(2, levels)),
        uLevels: ValueCell.create(levels),

        uGridDim: ValueCell.create(gridDimensions),
        uGridTexDim: ValueCell.create(Vec3.create(volumeData.width, volumeData.height, 0)),
        uGridTransform: ValueCell.create(transform),

        uScale: ValueCell.create(scale),
    }

    const schema = { ...IsosurfaceSchema }
    const shaderCode = ShaderCode(
        require('mol-gl/shader/quad.vert').default,
        require('mol-gl/shader/marching-cubes/isosurface.frag').default,
        { drawBuffers: true }
    )
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values)

    return createComputeRenderable(renderItem, values);
}

function setRenderingDefaults(gl: GLRenderingContext) {
    gl.disable(gl.CULL_FACE)
    gl.disable(gl.BLEND)
    gl.disable(gl.DEPTH_TEST)
    gl.depthMask(false)
}

export function createIsosurfaceBuffers(ctx: WebGLContext, activeVoxelsBase: Texture, volumeData: Texture, histogramPyramid: HistogramPyramid, gridDimensions: Vec3, transform: Mat4, isoValue: number) {
    const { gl, framebufferCache } = ctx
    const { pyramidTex, totalTex, height, levels, scale, count } = histogramPyramid

    const framebuffer = framebufferCache.get(FramebufferName).value
    framebuffer.bind()

    const vertexGroupTexture = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    vertexGroupTexture.define(pyramidTex.width, pyramidTex.height)

    const normalTexture = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    normalTexture.define(pyramidTex.width, pyramidTex.height)

    // const infoTex = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    // infoTex.define(pyramidTex.width, pyramidTex.height)

    // const pointTexA = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    // pointTexA.define(pyramidTex.width, pyramidTex.height)

    // const pointTexB = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    // pointTexB.define(pyramidTex.width, pyramidTex.height)

    // const coordTex = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    // coordTex.define(pyramidTex.width, pyramidTex.height)

    // const indexTex = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    // indexTex.define(pyramidTex.width, pyramidTex.height)

    const pr = getIsosurfaceRenderable(ctx, pyramidTex, activeVoxelsBase, volumeData, totalTex, gridDimensions, transform, isoValue, levels, scale)
    pr.update()
    pr.use()

    vertexGroupTexture.attachFramebuffer(framebuffer, 0)
    normalTexture.attachFramebuffer(framebuffer, 1)
    // infoTex.attachFramebuffer(framebuffer, 1)
    // pointTexA.attachFramebuffer(framebuffer, 2)
    // pointTexB.attachFramebuffer(framebuffer, 3)
    // coordTex.attachFramebuffer(framebuffer, 4)
    // indexTex.attachFramebuffer(framebuffer, 5)

    const { drawBuffers } = ctx.extensions
    if (!drawBuffers) throw new Error('need draw buffers')

    drawBuffers.drawBuffers([
        drawBuffers.COLOR_ATTACHMENT0,
        drawBuffers.COLOR_ATTACHMENT1,
        // drawBuffers.COLOR_ATTACHMENT2,
        // drawBuffers.COLOR_ATTACHMENT3,
        // drawBuffers.COLOR_ATTACHMENT4,
        // drawBuffers.COLOR_ATTACHMENT5
    ])

    setRenderingDefaults(gl)
    gl.viewport(0, 0, pyramidTex.width, pyramidTex.height)
    gl.scissor(0, 0, pyramidTex.width, height)
    pr.render()
    gl.disable(gl.SCISSOR_TEST)

    // const vgt = readTexture(ctx, vertexGroupTexture, pyramidTex.width, height)
    // console.log('vertexGroupTexture', vgt.array.subarray(0, 4 * count))

    // const vt = readTexture(ctx, verticesTex, pyramidTex.width, height)
    // console.log('vt', vt)
    // const vertices = new Float32Array(3 * compacted.count)
    // for (let i = 0; i < compacted.count; ++i) {
    //     vertices[i * 3] = vt.array[i * 4]
    //     vertices[i * 3 + 1] = vt.array[i * 4 + 1]
    //     vertices[i * 3 + 2] = vt.array[i * 4 + 2]
    // }
    // console.log('vertices', vertices)

    // const it = readTexture(ctx, infoTex, pyramidTex.width, height)
    // console.log('info', it.array.subarray(0, 4 * compacted.count))

    // const pat = readTexture(ctx, pointTexA, pyramidTex.width, height)
    // console.log('point a', pat.array.subarray(0, 4 * compacted.count))

    // const pbt = readTexture(ctx, pointTexB, pyramidTex.width, height)
    // console.log('point b', pbt.array.subarray(0, 4 * compacted.count))

    // const ct = readTexture(ctx, coordTex, pyramidTex.width, height)
    // console.log('coord', ct.array.subarray(0, 4 * compacted.count))

    // const idxt = readTexture(ctx, indexTex, pyramidTex.width, height)
    // console.log('index', idxt.array.subarray(0, 4 * compacted.count))

    // const { field, idField } = await fieldFromTexture2d(ctx, volumeData, gridDimensions)
    // console.log({ field, idField })

    // const valuesA = new Float32Array(compacted.count)
    // const valuesB = new Float32Array(compacted.count)
    // for (let i = 0; i < compacted.count; ++i) {
    //     valuesA[i] = field.space.get(field.data, pat.array[i * 4], pat.array[i * 4 + 1], pat.array[i * 4 + 2])
    //     valuesB[i] = field.space.get(field.data, pbt.array[i * 4], pbt.array[i * 4 + 1], pbt.array[i * 4 + 2])
    // }
    // console.log('valuesA', valuesA)
    // console.log('valuesB', valuesB)

    return { vertexGroupTexture, normalTexture, vertexCount: count }
}