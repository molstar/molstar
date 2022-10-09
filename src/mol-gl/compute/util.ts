/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { createNullTexture, Texture } from '../../mol-gl/webgl/texture';
import { ValueCell } from '../../mol-util';
import { ValueSpec, AttributeSpec, UniformSpec, Values, TextureSpec } from '../../mol-gl/renderable/schema';
import { Vec2 } from '../../mol-math/linear-algebra';
import { ShaderCode } from '../shader-code';
import { copy_frag } from '../shader/copy.frag';
import { quad_vert } from '../shader/quad.vert';
import { createComputeRenderItem } from '../webgl/render-item';
import { ComputeRenderable, createComputeRenderable } from '../renderable';

export const QuadPositions = new Float32Array([
    1.0, 1.0, -1.0, 1.0, -1.0, -1.0, // First triangle
    -1.0, -1.0, 1.0, -1.0, 1.0, 1.0 // Second triangle
]);

export const QuadSchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),
    aPosition: AttributeSpec('float32', 2, 0),
    uQuadScale: UniformSpec('v2'),
};

export const QuadValues: Values<typeof QuadSchema> = {
    drawCount: ValueCell.create(6),
    instanceCount: ValueCell.create(1),
    aPosition: ValueCell.create(QuadPositions),
    uQuadScale: ValueCell.create(Vec2.create(1, 1)),
};

//

const CopySchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),
};
const CopyShaderCode = ShaderCode('copy', quad_vert, copy_frag);
export type CopyRenderable = ComputeRenderable<Values<typeof CopySchema>>

export function createCopyRenderable(ctx: WebGLContext, texture: Texture): CopyRenderable {
    const values: Values<typeof CopySchema> = {
        ...QuadValues,
        tColor: ValueCell.create(texture),
        uTexSize: ValueCell.create(Vec2.create(texture.getWidth(), texture.getHeight())),
    };

    const schema = { ...CopySchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', CopyShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

const SharedCopyName = 'shared-copy';

export function getSharedCopyRenderable(ctx: WebGLContext, texture: Texture) {
    if (!ctx.namedComputeRenderables[SharedCopyName]) {
        ctx.namedComputeRenderables[SharedCopyName] = createCopyRenderable(ctx, createNullTexture());
    }
    const copy = ctx.namedComputeRenderables[SharedCopyName] as CopyRenderable;
    ValueCell.update(copy.values.tColor, texture);
    ValueCell.update(copy.values.uTexSize, Vec2.set(copy.values.uTexSize.ref.value, texture.getWidth(), texture.getHeight()));
    copy.update();
    return copy;
}

//

const ReadTextureName = 'read-texture';
const ReadAlphaTextureName = 'read-alpha-texture';

export function readTexture<T extends Uint8Array | Float32Array | Int32Array = Uint8Array>(ctx: WebGLContext, texture: Texture, array?: T) {
    const { gl, resources } = ctx;
    if (!array && texture.type !== gl.UNSIGNED_BYTE) throw new Error('unsupported texture type');

    if (!ctx.namedFramebuffers[ReadTextureName]) {
        ctx.namedFramebuffers[ReadTextureName] = resources.framebuffer();
    }
    const framebuffer = ctx.namedFramebuffers[ReadTextureName];

    const width = texture.getWidth();
    const height = texture.getHeight();
    if (!array) array = new Uint8Array(width * height * 4) as T;
    framebuffer.bind();
    texture.attachFramebuffer(framebuffer, 0);
    ctx.readPixels(0, 0, width, height, array);

    return { array, width, height };
}

export function readAlphaTexture(ctx: WebGLContext, texture: Texture) {
    const { gl, state, resources } = ctx;
    if (texture.type !== gl.UNSIGNED_BYTE) throw new Error('unsupported texture type');

    const width = texture.getWidth();
    const height = texture.getHeight();

    const copy = getSharedCopyRenderable(ctx, texture);
    state.currentRenderItemId = -1;

    if (!ctx.namedFramebuffers[ReadAlphaTextureName]) {
        ctx.namedFramebuffers[ReadAlphaTextureName] = resources.framebuffer();
    }
    const framebuffer = ctx.namedFramebuffers[ReadAlphaTextureName];
    framebuffer.bind();

    if (!ctx.namedTextures[ReadAlphaTextureName]) {
        ctx.namedTextures[ReadAlphaTextureName] = resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    }
    const copyTex = ctx.namedTextures[ReadAlphaTextureName];
    copyTex.define(width, height);
    copyTex.attachFramebuffer(framebuffer, 0);

    state.disable(gl.CULL_FACE);
    state.enable(gl.BLEND);
    state.disable(gl.DEPTH_TEST);
    state.enable(gl.SCISSOR_TEST);
    state.depthMask(false);
    state.clearColor(0, 0, 0, 0);
    state.blendFunc(gl.ONE, gl.ONE);
    state.blendEquation(gl.FUNC_ADD);
    state.viewport(0, 0, width, height);
    state.scissor(0, 0, width, height);
    gl.clear(gl.COLOR_BUFFER_BIT);
    copy.render();

    const array = new Uint8Array(width * height * 4);
    ctx.readPixels(0, 0, width, height, array);

    return { array, width, height };
}