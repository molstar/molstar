/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { printTextureImage } from '../../mol-gl/renderable/util';
import { defaults, ValueCell } from '../../mol-util';
import { ValueSpec, AttributeSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { Vec2 } from '../../mol-math/linear-algebra';
import { GLRenderingContext } from '../../mol-gl/webgl/compat';

export const QuadPositions = new Float32Array([
    1.0,  1.0,  -1.0,  1.0,  -1.0, -1.0, // First triangle
    -1.0, -1.0,   1.0, -1.0,   1.0,  1.0  // Second triangle
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

function getArrayForTexture(gl: GLRenderingContext, texture: Texture, size: number) {
    switch (texture.type) {
        case gl.UNSIGNED_BYTE: return new Uint8Array(size);
        case gl.FLOAT: return new Float32Array(size);
    }
    throw new Error('unknown/unsupported texture type');
}

export function readTexture(ctx: WebGLContext, texture: Texture, width?: number, height?: number) {
    const { gl, resources } = ctx;
    width = defaults(width, texture.getWidth());
    height = defaults(height, texture.getHeight());
    const size = width * height * 4;
    const framebuffer = resources.framebuffer();
    const array = getArrayForTexture(gl, texture, size);
    framebuffer.bind();
    texture.attachFramebuffer(framebuffer, 0);
    ctx.readPixels(0, 0, width, height, array);

    return { array, width, height };
}

export function printTexture(ctx: WebGLContext, texture: Texture, scale: number) {
    printTextureImage(readTexture(ctx, texture), scale);
}