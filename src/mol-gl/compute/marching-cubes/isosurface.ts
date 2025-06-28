/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ComputeRenderable, createComputeRenderable } from '../../renderable';
import { WebGLContext } from '../../webgl/context';
import { createComputeRenderItem } from '../../webgl/render-item';
import { Values, TextureSpec, UniformSpec, DefineSpec } from '../../renderable/schema';
import { Texture } from '../../../mol-gl/webgl/texture';
import { ShaderCode } from '../../../mol-gl/shader-code';
import { ValueCell } from '../../../mol-util';
import { Vec3, Vec2, Mat4, Mat3 } from '../../../mol-math/linear-algebra';
import { QuadSchema, QuadValues } from '../util';
import { createHistogramPyramid, HistogramPyramid } from '../histogram-pyramid/reduction';
import { getTriIndices } from './tables';
import { quad_vert } from '../../../mol-gl/shader/quad.vert';
import { isosurface_frag } from '../../../mol-gl/shader/marching-cubes/isosurface.frag';
import { calcActiveVoxels } from './active-voxels';
import { isWebGL2 } from '../../webgl/compat';
import { isTimingMode } from '../../../mol-util/debug';

const IsosurfaceSchema = {
    ...QuadSchema,

    tTriIndices: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    tActiveVoxelsPyramid: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    tActiveVoxelsBase: TextureSpec('texture', 'rgba', 'float', 'nearest'),
    tVolumeData: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    dValueChannel: DefineSpec('string', ['red', 'alpha']),
    uIsoValue: UniformSpec('f'),

    uSize: UniformSpec('f'),
    uLevels: UniformSpec('f'),
    uCount: UniformSpec('f'),
    uInvert: UniformSpec('b'),

    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v3'),
    uGridTransform: UniformSpec('m4'),
    uGridTransformAdjoint: UniformSpec('m3'),
    uScale: UniformSpec('v2'),

    dPackedGroup: DefineSpec('boolean'),
    dAxisOrder: DefineSpec('string', ['012', '021', '102', '120', '201', '210']),
    dConstantGroup: DefineSpec('boolean'),
};
type IsosurfaceValues = Values<typeof IsosurfaceSchema>

const IsosurfaceName = 'isosurface';

function valueChannel(ctx: WebGLContext, volumeData: Texture) {
    return isWebGL2(ctx.gl) && volumeData.format === ctx.gl.RED ? 'red' : 'alpha';
}

function getIsosurfaceRenderable(ctx: WebGLContext, activeVoxelsPyramid: Texture, activeVoxelsBase: Texture, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, transform: Mat4, isoValue: number, levels: number, scale: Vec2, count: number, invert: boolean, packedGroup: boolean, axisOrder: Vec3, constantGroup: boolean): ComputeRenderable<IsosurfaceValues> {
    if (ctx.namedComputeRenderables[IsosurfaceName]) {
        const v = ctx.namedComputeRenderables[IsosurfaceName].values as IsosurfaceValues;

        ValueCell.update(v.tActiveVoxelsPyramid, activeVoxelsPyramid);
        ValueCell.update(v.tActiveVoxelsBase, activeVoxelsBase);
        ValueCell.update(v.tVolumeData, volumeData);
        ValueCell.update(v.dValueChannel, valueChannel(ctx, volumeData));

        ValueCell.updateIfChanged(v.uIsoValue, isoValue);
        ValueCell.updateIfChanged(v.uSize, Math.pow(2, levels));
        ValueCell.updateIfChanged(v.uLevels, levels);
        ValueCell.updateIfChanged(v.uCount, count);
        ValueCell.updateIfChanged(v.uInvert, invert);

        ValueCell.update(v.uGridDim, gridDim);
        ValueCell.update(v.uGridTexDim, gridTexDim);
        ValueCell.update(v.uGridTransform, transform);
        ValueCell.update(v.uGridTransformAdjoint, Mat3.adjointFromMat4(Mat3(), transform));
        ValueCell.update(v.uScale, scale);

        ValueCell.updateIfChanged(v.dPackedGroup, packedGroup);
        ValueCell.updateIfChanged(v.dAxisOrder, axisOrder.join(''));
        ValueCell.updateIfChanged(v.dConstantGroup, constantGroup);

        ctx.namedComputeRenderables[IsosurfaceName].update();
    } else {
        ctx.namedComputeRenderables[IsosurfaceName] = createIsosurfaceRenderable(ctx, activeVoxelsPyramid, activeVoxelsBase, volumeData, gridDim, gridTexDim, transform, isoValue, levels, scale, count, invert, packedGroup, axisOrder, constantGroup);
    }
    return ctx.namedComputeRenderables[IsosurfaceName];
}

function createIsosurfaceRenderable(ctx: WebGLContext, activeVoxelsPyramid: Texture, activeVoxelsBase: Texture, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, transform: Mat4, isoValue: number, levels: number, scale: Vec2, count: number, invert: boolean, packedGroup: boolean, axisOrder: Vec3, constantGroup: boolean) {
    // console.log('uSize', Math.pow(2, levels))
    const values: IsosurfaceValues = {
        ...QuadValues,
        tTriIndices: ValueCell.create(getTriIndices()),

        tActiveVoxelsPyramid: ValueCell.create(activeVoxelsPyramid),
        tActiveVoxelsBase: ValueCell.create(activeVoxelsBase),
        tVolumeData: ValueCell.create(volumeData),
        dValueChannel: ValueCell.create(valueChannel(ctx, volumeData)),

        uIsoValue: ValueCell.create(isoValue),
        uSize: ValueCell.create(Math.pow(2, levels)),
        uLevels: ValueCell.create(levels),
        uCount: ValueCell.create(count),
        uInvert: ValueCell.create(invert),

        uGridDim: ValueCell.create(gridDim),
        uGridTexDim: ValueCell.create(gridTexDim),
        uGridTransform: ValueCell.create(transform),
        uGridTransformAdjoint: ValueCell.create(Mat3.adjointFromMat4(Mat3(), transform)),
        uScale: ValueCell.create(scale),

        dPackedGroup: ValueCell.create(packedGroup),
        dAxisOrder: ValueCell.create(axisOrder.join('')),
        dConstantGroup: ValueCell.create(constantGroup),
    };

    const schema = { ...IsosurfaceSchema };
    const shaderCode = ShaderCode('isosurface', quad_vert, isosurface_frag, { drawBuffers: 'required' });
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

export function createIsosurfaceBuffers(ctx: WebGLContext, activeVoxelsBase: Texture, volumeData: Texture, histogramPyramid: HistogramPyramid, gridDim: Vec3, gridTexDim: Vec3, transform: Mat4, isoValue: number, invert: boolean, packedGroup: boolean, axisOrder: Vec3, constantGroup: boolean, vertexTexture?: Texture, groupTexture?: Texture, normalTexture?: Texture) {
    const { drawBuffers } = ctx.extensions;
    if (!drawBuffers) throw new Error('need WebGL draw buffers');

    if (isTimingMode) ctx.timer.mark('createIsosurfaceBuffers');
    const { gl, state, resources, extensions } = ctx;
    const { pyramidTex, height, levels, scale, count } = histogramPyramid;
    const width = pyramidTex.getWidth();

    // console.log('width', width, 'height', height);
    // console.log('iso', 'gridDim', gridDim, 'scale', scale, 'gridTexDim', gridTexDim);
    // console.log('iso volumeData', volumeData);

    if (!ctx.namedFramebuffers[IsosurfaceName]) {
        ctx.namedFramebuffers[IsosurfaceName] = resources.framebuffer();
    }
    const framebuffer = ctx.namedFramebuffers[IsosurfaceName];

    if (isWebGL2(gl)) {
        if (!vertexTexture) {
            vertexTexture = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        }

        if (!groupTexture) {
            groupTexture = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        }

        if (!normalTexture) {
            normalTexture = extensions.colorBufferHalfFloat && extensions.textureHalfFloat
                ? resources.texture('image-float16', 'rgba', 'fp16', 'nearest')
                : resources.texture('image-float32', 'rgba', 'float', 'nearest');
        }
    } else {
        // webgl1 requires consistent bit plane counts
        // this is quite wasteful but good enough for medium size meshes

        if (!vertexTexture) {
            vertexTexture = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        }

        if (!groupTexture) {
            groupTexture = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        }

        if (!normalTexture) {
            normalTexture = resources.texture('image-float32', 'rgba', 'float', 'nearest');
        }
    }

    vertexTexture.define(width, height);
    groupTexture.define(width, height);
    normalTexture.define(width, height);

    vertexTexture.attachFramebuffer(framebuffer, 0);
    groupTexture.attachFramebuffer(framebuffer, 1);
    normalTexture.attachFramebuffer(framebuffer, 2);

    const renderable = getIsosurfaceRenderable(ctx, pyramidTex, activeVoxelsBase, volumeData, gridDim, gridTexDim, transform, isoValue, levels, scale, count, invert, packedGroup, axisOrder, constantGroup);
    ctx.state.currentRenderItemId = -1;

    framebuffer.bind();
    drawBuffers.drawBuffers([
        drawBuffers.COLOR_ATTACHMENT0,
        drawBuffers.COLOR_ATTACHMENT1,
        drawBuffers.COLOR_ATTACHMENT2,
    ]);

    setRenderingDefaults(ctx);
    state.viewport(0, 0, width, height);
    gl.clear(gl.COLOR_BUFFER_BIT);
    renderable.render();

    gl.finish();
    if (isTimingMode) ctx.timer.markEnd('createIsosurfaceBuffers');

    // printTextureImage(readTexture(ctx, vertexTexture, new Float32Array(width * height * 4)), { scale: 0.75, normalize: true });
    // printTextureImage(readTexture(ctx, groupTexture, new Uint8Array(width * height * 4)), { scale: 0.75, normalize: true });
    // printTextureImage(readTexture(ctx, normalTexture, new Float32Array(width * height * 4)), { scale: 0.75, normalize: true });

    return { vertexTexture, groupTexture, normalTexture, vertexCount: count };
}

//

/**
 * GPU isosurface extraction
 *
 * Algorithm from "High‐speed Marching Cubes using HistoPyramids"
 * by C Dyken, G Ziegler, C Theobalt, HP Seidel
 * https://doi.org/10.1111/j.1467-8659.2008.01182.x
 *
 * Implementation based on http://www.miaumiau.cat/2016/10/stream-compaction-in-webgl/
 */
export function extractIsosurface(ctx: WebGLContext, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, gridTexScale: Vec2, transform: Mat4, isoValue: number, invert: boolean, packedGroup: boolean, axisOrder: Vec3, constantGroup: boolean, vertexTexture?: Texture, groupTexture?: Texture, normalTexture?: Texture) {
    if (isTimingMode) ctx.timer.mark('extractIsosurface');
    const activeVoxelsTex = calcActiveVoxels(ctx, volumeData, gridDim, gridTexDim, isoValue, gridTexScale);
    const compacted = createHistogramPyramid(ctx, activeVoxelsTex, gridTexScale, gridTexDim);
    const gv = createIsosurfaceBuffers(ctx, activeVoxelsTex, volumeData, compacted, gridDim, gridTexDim, transform, isoValue, invert, packedGroup, axisOrder, constantGroup, vertexTexture, groupTexture, normalTexture);
    if (isTimingMode) ctx.timer.markEnd('extractIsosurface');

    return gv;
}