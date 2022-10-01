/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ComputeRenderable, createComputeRenderable } from '../../renderable';
import { WebGLContext } from '../../webgl/context';
import { createComputeRenderItem } from '../../webgl/render-item';
import { Values, TextureSpec, UniformSpec } from '../../renderable/schema';
import { Texture } from '../../../mol-gl/webgl/texture';
import { ShaderCode } from '../../../mol-gl/shader-code';
import { ValueCell } from '../../../mol-util';
import { Vec3, Vec2 } from '../../../mol-math/linear-algebra';
import { QuadSchema, QuadValues } from '../util';
import { getTriCount } from './tables';
import { quad_vert } from '../../../mol-gl/shader/quad.vert';
import { activeVoxels_frag } from '../../../mol-gl/shader/marching-cubes/active-voxels.frag';
import { isTimingMode } from '../../../mol-util/debug';

const ActiveVoxelsSchema = {
    ...QuadSchema,

    tTriCount: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    tVolumeData: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uIsoValue: UniformSpec('f'),

    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v3'),

    uScale: UniformSpec('v2'),
};
type ActiveVoxelsValues = Values<typeof ActiveVoxelsSchema>

const ActiveVoxelsName = 'active-voxels';

function getActiveVoxelsRenderable(ctx: WebGLContext, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, isoValue: number, scale: Vec2): ComputeRenderable<ActiveVoxelsValues> {
    if (ctx.namedComputeRenderables[ActiveVoxelsName]) {
        const v = ctx.namedComputeRenderables[ActiveVoxelsName].values as ActiveVoxelsValues;

        ValueCell.update(v.uQuadScale, scale);
        ValueCell.update(v.tVolumeData, volumeData);
        ValueCell.updateIfChanged(v.uIsoValue, isoValue);
        ValueCell.update(v.uGridDim, gridDim);
        ValueCell.update(v.uGridTexDim, gridTexDim);
        ValueCell.update(v.uScale, scale);

        ctx.namedComputeRenderables[ActiveVoxelsName].update();
    } else {
        ctx.namedComputeRenderables[ActiveVoxelsName] = createActiveVoxelsRenderable(ctx, volumeData, gridDim, gridTexDim, isoValue, scale);
    }
    return ctx.namedComputeRenderables[ActiveVoxelsName];
}

function createActiveVoxelsRenderable(ctx: WebGLContext, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, isoValue: number, scale: Vec2) {
    const values: ActiveVoxelsValues = {
        ...QuadValues,
        tTriCount: ValueCell.create(getTriCount()),

        uQuadScale: ValueCell.create(scale),
        tVolumeData: ValueCell.create(volumeData),
        uIsoValue: ValueCell.create(isoValue),
        uGridDim: ValueCell.create(gridDim),
        uGridTexDim: ValueCell.create(gridTexDim),
        uScale: ValueCell.create(scale),
    };

    const schema = { ...ActiveVoxelsSchema };
    const shaderCode = ShaderCode('active-voxels', quad_vert, activeVoxels_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

function setRenderingDefaults(ctx: WebGLContext) {
    const { gl, state } = ctx;
    state.disable(gl.CULL_FACE);
    state.disable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.enable(gl.SCISSOR_TEST);
    state.depthMask(false);
    state.colorMask(true, true, true, true);
    state.clearColor(0, 0, 0, 0);
}

export function calcActiveVoxels(ctx: WebGLContext, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, isoValue: number, gridScale: Vec2) {
    if (isTimingMode) ctx.timer.mark('calcActiveVoxels');
    const { gl, state, resources } = ctx;
    const width = volumeData.getWidth();
    const height = volumeData.getHeight();

    if (!ctx.namedFramebuffers[ActiveVoxelsName]) {
        ctx.namedFramebuffers[ActiveVoxelsName] = resources.framebuffer();
    }
    const framebuffer = ctx.namedFramebuffers[ActiveVoxelsName];
    framebuffer.bind();

    if (!ctx.namedTextures[ActiveVoxelsName]) {
        ctx.namedTextures[ActiveVoxelsName] = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    }
    const activeVoxelsTex = ctx.namedTextures[ActiveVoxelsName];
    activeVoxelsTex.define(width, height);

    const renderable = getActiveVoxelsRenderable(ctx, volumeData, gridDim, gridTexDim, isoValue, gridScale);
    ctx.state.currentRenderItemId = -1;

    activeVoxelsTex.attachFramebuffer(framebuffer, 0);
    setRenderingDefaults(ctx);
    state.viewport(0, 0, width, height);
    state.scissor(0, 0, width, height);
    gl.clear(gl.COLOR_BUFFER_BIT);
    state.scissor(0, 0, gridTexDim[0], gridTexDim[1]);
    renderable.render();

    // console.log('gridScale', gridScale, 'gridTexDim', gridTexDim, 'gridDim', gridDim);
    // console.log('volumeData', volumeData);
    // console.log('at', readTexture(ctx, activeVoxelsTex));
    // printTextureImage(readTexture(ctx, activeVoxelsTex), { scale: 0.75 });

    gl.finish();
    if (isTimingMode) ctx.timer.markEnd('calcActiveVoxels');

    return activeVoxelsTex;
}