/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createComputeRenderable } from '../../renderable';
import { WebGLContext } from '../../webgl/context';
import { createComputeRenderItem } from '../../webgl/render-item';
import { Values, TextureSpec, UniformSpec } from '../../renderable/schema';
import { Texture } from '../../../mol-gl/webgl/texture';
import { ShaderCode } from '../../../mol-gl/shader-code';
import { ValueCell } from '../../../mol-util';
import { Vec3, Vec2, Mat4 } from '../../../mol-math/linear-algebra';
import { QuadSchema, QuadValues } from '../util';
import { HistogramPyramid } from '../histogram-pyramid/reduction';
import { getTriIndices } from './tables';
import quad_vert from '../../../mol-gl/shader/quad.vert';
import isosurface_frag from '../../../mol-gl/shader/marching-cubes/isosurface.frag';

const IsosurfaceSchema = {
    ...QuadSchema,

    tTriIndices: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    tActiveVoxelsPyramid: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    tActiveVoxelsBase: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    tVolumeData: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uIsoValue: UniformSpec('f'),

    uSize: UniformSpec('f'),
    uLevels: UniformSpec('f'),
    uCount: UniformSpec('f'),

    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v3'),
    uGridTransform: UniformSpec('m4'),

    uScale: UniformSpec('v2'),
};

function getIsosurfaceRenderable(ctx: WebGLContext, activeVoxelsPyramid: Texture, activeVoxelsBase: Texture, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, transform: Mat4, isoValue: number, levels: number, scale: Vec2, count: number, height: number) {
    // console.log('uSize', Math.pow(2, levels))
    const values: Values<typeof IsosurfaceSchema> = {
        ...QuadValues,
        uQuadScale: ValueCell.create(Vec2.create(1, height / Math.pow(2, levels))),

        tTriIndices: ValueCell.create(getTriIndices()),
        tActiveVoxelsPyramid: ValueCell.create(activeVoxelsPyramid),
        tActiveVoxelsBase: ValueCell.create(activeVoxelsBase),
        tVolumeData: ValueCell.create(volumeData),
        uIsoValue: ValueCell.create(isoValue),

        uSize: ValueCell.create(Math.pow(2, levels)),
        uLevels: ValueCell.create(levels),
        uCount: ValueCell.create(count),

        uGridDim: ValueCell.create(gridDim),
        uGridTexDim: ValueCell.create(gridTexDim),
        uGridTransform: ValueCell.create(transform),

        uScale: ValueCell.create(scale),
    };

    const schema = { ...IsosurfaceSchema };
    const shaderCode = ShaderCode('isosurface', quad_vert, isosurface_frag, { drawBuffers: true });
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function setRenderingDefaults(ctx: WebGLContext) {
    const { gl, state } = ctx;
    state.disable(gl.CULL_FACE);
    state.disable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.disable(gl.SCISSOR_TEST);
    state.depthMask(false);
    state.colorMask(true, true, true, true);
    state.clearColor(0, 0, 0, 0);
}

export function createIsosurfaceBuffers(ctx: WebGLContext, activeVoxelsBase: Texture, volumeData: Texture, histogramPyramid: HistogramPyramid, gridDim: Vec3, gridTexDim: Vec3, transform: Mat4, isoValue: number, vertexGroupTexture?: Texture, normalTexture?: Texture) {
    const { gl, resources } = ctx;
    const { pyramidTex, height, levels, scale, count } = histogramPyramid;

    // console.log('iso', 'gridDim', gridDim, 'scale', scale, 'gridTexDim', gridTexDim)
    // console.log('iso volumeData', volumeData)

    const framebuffer = resources.framebuffer();

    let needsClear = false;

    if (!vertexGroupTexture) {
        vertexGroupTexture = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        vertexGroupTexture.define(pyramidTex.getWidth(), pyramidTex.getHeight());
    } else if (vertexGroupTexture.getWidth() !== pyramidTex.getWidth() || vertexGroupTexture.getHeight() !== pyramidTex.getHeight()) {
        vertexGroupTexture.define(pyramidTex.getWidth(), pyramidTex.getHeight());
    } else {
        needsClear = true;
    }

    if (!normalTexture) {
        normalTexture = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        normalTexture.define(pyramidTex.getWidth(), pyramidTex.getHeight());
    } else if (normalTexture.getWidth() !== pyramidTex.getWidth() || normalTexture.getHeight() !== pyramidTex.getHeight()) {
        normalTexture.define(pyramidTex.getWidth(), pyramidTex.getHeight());
    } else {
        needsClear = true;
    }

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

    const renderable = getIsosurfaceRenderable(ctx, pyramidTex, activeVoxelsBase, volumeData, gridDim, gridTexDim, transform, isoValue, levels, scale, count, height);
    ctx.state.currentRenderItemId = -1;

    vertexGroupTexture.attachFramebuffer(framebuffer, 0);
    normalTexture.attachFramebuffer(framebuffer, 1);
    // infoTex.attachFramebuffer(framebuffer, 1)
    // pointTexA.attachFramebuffer(framebuffer, 2)
    // pointTexB.attachFramebuffer(framebuffer, 3)
    // coordTex.attachFramebuffer(framebuffer, 4)
    // indexTex.attachFramebuffer(framebuffer, 5)

    const { drawBuffers } = ctx.extensions;
    if (!drawBuffers) throw new Error('need WebGL draw buffers');

    drawBuffers.drawBuffers([
        drawBuffers.COLOR_ATTACHMENT0,
        drawBuffers.COLOR_ATTACHMENT1,
        // drawBuffers.COLOR_ATTACHMENT2,
        // drawBuffers.COLOR_ATTACHMENT3,
        // drawBuffers.COLOR_ATTACHMENT4,
        // drawBuffers.COLOR_ATTACHMENT5
    ]);

    setRenderingDefaults(ctx);
    gl.viewport(0, 0, pyramidTex.getWidth(), pyramidTex.getHeight());
    if (needsClear) gl.clear(gl.COLOR_BUFFER_BIT);
    renderable.render();

    gl.finish();

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

    return { vertexGroupTexture, normalTexture, vertexCount: count };
}