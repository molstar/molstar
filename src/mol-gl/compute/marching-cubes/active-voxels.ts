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
import { Vec3, Vec2 } from '../../../mol-math/linear-algebra';
import { QuadSchema, QuadValues } from '../util';
import { getTriCount } from './tables';
import quad_vert from '../../../mol-gl/shader/quad.vert';
import active_voxels_frag from '../../../mol-gl/shader/marching-cubes/active-voxels.frag';

const ActiveVoxelsSchema = {
    ...QuadSchema,

    tTriCount: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    tVolumeData: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uIsoValue: UniformSpec('f'),

    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v3'),

    uScale: UniformSpec('v2'),
};

function getActiveVoxelsRenderable(ctx: WebGLContext, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, isoValue: number, scale: Vec2) {
    const values: Values<typeof ActiveVoxelsSchema> = {
        ...QuadValues,
        uQuadScale: ValueCell.create(scale),

        tTriCount: ValueCell.create(getTriCount()),
        tVolumeData: ValueCell.create(volumeData),
        uIsoValue: ValueCell.create(isoValue),

        uGridDim: ValueCell.create(gridDim),
        uGridTexDim: ValueCell.create(gridTexDim),

        uScale: ValueCell.create(scale),
    };

    const schema = { ...ActiveVoxelsSchema };
    const shaderCode = ShaderCode('active-voxels', quad_vert, active_voxels_frag);
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

export function calcActiveVoxels(ctx: WebGLContext, volumeData: Texture, gridDim: Vec3, gridTexDim: Vec3, isoValue: number, gridScale: Vec2) {
    const { gl, resources } = ctx;
    const width = volumeData.getWidth();
    const height = volumeData.getHeight();

    const framebuffer = resources.framebuffer();
    framebuffer.bind();

    const activeVoxelsTex = resources.texture('image-float32', 'rgba', 'float', 'nearest');
    activeVoxelsTex.define(width, height);

    const renderable = getActiveVoxelsRenderable(ctx, volumeData, gridDim, gridTexDim, isoValue, gridScale);
    ctx.state.currentRenderItemId = -1;

    activeVoxelsTex.attachFramebuffer(framebuffer, 0);
    setRenderingDefaults(ctx);
    gl.viewport(0, 0, width, height);
    renderable.render();

    // console.log('gridScale', gridScale, 'gridTexDim', gridTexDim, 'gridDim', gridDim)
    // console.log('volumeData', volumeData)
    // console.log('at', readTexture(ctx, activeVoxelsTex))

    gl.finish();

    return activeVoxelsTex;
}