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
import { Vec3 } from 'mol-math/linear-algebra';
import { QuadSchema, QuadValues } from '../util';
import { getTriCount } from './tables';

/** name for shared framebuffer used for gpu marching cubes operations */
const FramebufferName = 'marching-cubes-active-voxels'

const ActiveVoxelsSchema = {
    ...QuadSchema,

    tTriCount: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    tVolumeData: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uIsoValue: UniformSpec('f'),

    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v3'),
}

function getActiveVoxelsRenderable(ctx: WebGLContext, volumeData: Texture, gridDimensions: Vec3, isoValue: number) {
    const values: Values<typeof ActiveVoxelsSchema> = {
        ...QuadValues,

        tTriCount: ValueCell.create(getTriCount()),
        tVolumeData: ValueCell.create(volumeData),
        uIsoValue: ValueCell.create(isoValue),

        uGridDim: ValueCell.create(gridDimensions),
        uGridTexDim: ValueCell.create(Vec3.create(volumeData.width, volumeData.height, 0)),
    }

    const schema = { ...ActiveVoxelsSchema }
    const shaderCode = ShaderCode(
        require('mol-gl/shader/quad.vert').default,
        require('mol-gl/shader/marching-cubes/active-voxels.frag').default
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

export function calcActiveVoxels(ctx: WebGLContext, cornerTex: Texture, gridDimensions: Vec3, isoValue: number) {
    const { gl, framebufferCache } = ctx
    const { width, height } = cornerTex

    const framebuffer = framebufferCache.get(FramebufferName).value
    framebuffer.bind()

    const activeVoxelsTex = createTexture(ctx, 'image-float32', 'rgba', 'float', 'nearest')
    activeVoxelsTex.define(width, height)

    const renderable = getActiveVoxelsRenderable(ctx, cornerTex, gridDimensions, isoValue)
    renderable.update()

    activeVoxelsTex.attachFramebuffer(framebuffer, 0)
    setRenderingDefaults(gl)
    gl.viewport(0, 0, width, height)
    renderable.render()

    // const at = readTexture(ctx, activeVoxelsTex)
    // console.log('at', at)

    return activeVoxelsTex
}