/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { WebGLContext } from 'mol-gl/webgl/context';
import { Texture } from 'mol-gl/webgl/texture';
import { printTextureImage } from 'mol-gl/renderable/util';
import { defaults, ValueCell } from 'mol-util';
import { ValueSpec, AttributeSpec, UniformSpec } from 'mol-gl/renderable/schema';
import { Vec2 } from 'mol-math/linear-algebra';

export const QuadPositions = new Float32Array([
     1.0,  1.0,  -1.0,  1.0,  -1.0, -1.0, // First triangle
    -1.0, -1.0,   1.0, -1.0,   1.0,  1.0  // Second triangle
])

export const QuadSchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),
    aPosition: AttributeSpec('float32', 2, 0),
    uScale: UniformSpec('v2'),
}

export const QuadValues = {
    drawCount: ValueCell.create(6),
    instanceCount: ValueCell.create(1),
    aPosition: ValueCell.create(QuadPositions),
    uScale: ValueCell.create(Vec2.create(1, 1)),
}

//

export function readTexture(ctx: WebGLContext, texture: Texture, width?: number, height?: number) {
    const { gl, framebufferCache } = ctx
    width = defaults(width, texture.width)
    height = defaults(height, texture.height)
    const framebuffer = framebufferCache.get('read-texture').value
    const array = texture.type === gl.UNSIGNED_BYTE ? new Uint8Array(width * height * 4) : new Float32Array(width * height * 4)
    framebuffer.bind()
    texture.attachFramebuffer(framebuffer, 0)
    ctx.readPixels(0, 0, width, height, array)

    return { array, width, height }
}

export function printTexture(ctx: WebGLContext, texture: Texture, scale: number) {
    printTextureImage(readTexture(ctx, texture), scale)
}