/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ComputeRenderable, createComputeRenderable } from '../../renderable';
import { WebGLContext } from '../../webgl/context';
import { createComputeRenderItem } from '../../webgl/render-item';
import { Values, TextureSpec } from '../../renderable/schema';
import { Texture } from '../../../mol-gl/webgl/texture';
import { ShaderCode } from '../../../mol-gl/shader-code';
import { ValueCell } from '../../../mol-util';
import { unpackRGBToInt } from '../../../mol-util/number-packing';
import { QuadSchema, QuadValues } from '../util';
import { quad_vert } from '../../../mol-gl/shader/quad.vert';
import { sum_frag } from '../../../mol-gl/shader/histogram-pyramid/sum.frag';
import { isWebGL2 } from '../../webgl/compat';
import { isTimingMode } from '../../../mol-util/debug';

const HistopyramidSumSchema = {
    ...QuadSchema,
    tTexture: TextureSpec('texture', 'rgba', 'float', 'nearest'),
};
type HistopyramidSumValues = Values<typeof HistopyramidSumSchema>

const HistopyramidSumName = 'histopyramid-sum';

function getHistopyramidSumRenderable(ctx: WebGLContext, texture: Texture): ComputeRenderable<HistopyramidSumValues> {
    if (ctx.namedComputeRenderables[HistopyramidSumName]) {
        const v = ctx.namedComputeRenderables[HistopyramidSumName].values as HistopyramidSumValues;

        ValueCell.update(v.tTexture, texture);

        ctx.namedComputeRenderables[HistopyramidSumName].update();
    } else {
        ctx.namedComputeRenderables[HistopyramidSumName] = createHistopyramidSumRenderable(ctx, texture);
    }
    return ctx.namedComputeRenderables[HistopyramidSumName];
}

function createHistopyramidSumRenderable(ctx: WebGLContext, texture: Texture) {
    const values: HistopyramidSumValues = {
        ...QuadValues,
        tTexture: ValueCell.create(texture),
    };

    const schema = { ...HistopyramidSumSchema };
    const shaderCode = ShaderCode('sum', quad_vert, sum_frag, {}, { 0: 'ivec4' });
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

const sumBytes = new Uint8Array(4);
const sumInts = new Int32Array(4);

export function getHistopyramidSum(ctx: WebGLContext, pyramidTopTexture: Texture) {
    if (isTimingMode) ctx.timer.mark('getHistopyramidSum');
    const { gl, state, resources } = ctx;

    const renderable = getHistopyramidSumRenderable(ctx, pyramidTopTexture);
    ctx.state.currentRenderItemId = -1;

    if (!ctx.namedFramebuffers[HistopyramidSumName]) {
        ctx.namedFramebuffers[HistopyramidSumName] = resources.framebuffer();
    }
    const framebuffer = ctx.namedFramebuffers[HistopyramidSumName];

    if (!ctx.namedTextures[HistopyramidSumName]) {
        ctx.namedTextures[HistopyramidSumName] = isWebGL2(gl)
            ? resources.texture('image-int32', 'rgba', 'int', 'nearest')
            : resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
        ctx.namedTextures[HistopyramidSumName].define(1, 1);
    }
    const sumTexture = ctx.namedTextures[HistopyramidSumName];
    sumTexture.attachFramebuffer(framebuffer, 0);

    setRenderingDefaults(ctx);

    state.viewport(0, 0, 1, 1);
    renderable.render();
    gl.finish();

    ctx.readPixels(0, 0, 1, 1, isWebGL2(gl) ? sumInts : sumBytes);
    ctx.unbindFramebuffer();
    if (isTimingMode) ctx.timer.markEnd('getHistopyramidSum');

    return isWebGL2(gl)
        ? sumInts[0]
        : unpackRGBToInt(sumBytes[0], sumBytes[1], sumBytes[2]);
}