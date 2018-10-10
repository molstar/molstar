/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context'
import { TextureImage, TextureVolume } from '../renderable/util';
import { ValueCell } from 'mol-util';
import { RenderableSchema } from '../renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { Framebuffer } from './framebuffer';

const getNextTextureId = idFactory()

export type TextureKindValue = {
    'image-uint8': TextureImage<Uint8Array>
    'image-float32': TextureImage<Float32Array>
    'volume-uint8': TextureVolume<Uint8Array>
    'volume-float32': TextureVolume<Float32Array>
}
export type TextureKind = keyof TextureKindValue
export type TextureType = 'ubyte' | 'float'
export type TextureFormat = 'alpha' | 'rgb' | 'rgba'
export type TextureAttachment = 'depth' | 'stencil' | 'color0'
export type TextureFilter = 'nearest' | 'linear'

export function getTarget(ctx: Context, kind: TextureKind): number {
    const { gl } = ctx
    switch (kind) {
        case 'image-uint8': return gl.TEXTURE_2D
        case 'image-float32': return gl.TEXTURE_2D
        case 'volume-uint8': return (gl as WebGL2RenderingContext).TEXTURE_3D
        case 'volume-float32': return (gl as WebGL2RenderingContext).TEXTURE_3D
    }
}

export function getFormat(ctx: Context, format: TextureFormat): number {
    const { gl } = ctx
    switch (format) {
        case 'alpha': return gl.ALPHA
        case 'rgb': return gl.RGB
        case 'rgba': return gl.RGBA
    }
}

export function getInternalFormat(ctx: Context, format: TextureFormat, type: TextureType): number {
    const { gl, isWebGL2 } = ctx
    if (isWebGL2) {
        switch (format) {
            case 'alpha':
                switch (type) {
                    case 'ubyte': return gl.ALPHA
                    case 'float': throw new Error('invalid format/type combination alpha/float')
                }
            case 'rgb':
                switch (type) {
                    case 'ubyte': return gl.RGB
                    case 'float': return (gl as WebGL2RenderingContext).RGB32F
                }
            case 'rgba':
                switch (type) {
                    case 'ubyte': return gl.RGBA
                    case 'float': return (gl as WebGL2RenderingContext).RGBA32F
                }
        }
    }
    return getFormat(ctx, format)
}

export function getType(ctx: Context, type: TextureType): number {
    const { gl } = ctx
    switch (type) {
        case 'ubyte': return gl.UNSIGNED_BYTE
        case 'float': return gl.FLOAT
    }
}

export function getFilter(ctx: Context, type: TextureFilter): number {
    const { gl } = ctx
    switch (type) {
        case 'nearest': return gl.NEAREST
        case 'linear': return gl.LINEAR
    }
}

export function getAttachment(ctx: Context, attachment: TextureAttachment): number {
    const { gl } = ctx
    switch (attachment) {
        case 'depth': return gl.DEPTH_ATTACHMENT
        case 'stencil': return gl.STENCIL_ATTACHMENT
        case 'color0': return gl.COLOR_ATTACHMENT0
    }
}

export interface Texture {
    readonly id: number
    readonly target: number
    readonly format: number
    readonly internalFormat: number
    readonly type: number

    load: (image: TextureImage<any>) => void
    bind: (id: TextureId) => void
    unbind: (id: TextureId) => void
    attachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment) => void
    destroy: () => void
}

export type TextureId = 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15

export type TextureValues = { [k: string]: ValueCell<TextureImage<any>> }
export type Textures = { [k: string]: Texture }

export function createTexture(ctx: Context, kind: TextureKind, _format: TextureFormat, _type: TextureType, _filter: TextureFilter): Texture {
    const id = getNextTextureId()
    const { gl } = ctx
    const texture = gl.createTexture()
    if (texture === null) {
        throw new Error('Could not create WebGL texture')
    }

    const target = getTarget(ctx, kind)
    const filter = getFilter(ctx, _filter)
    const format = getFormat(ctx, _format)
    const internalFormat = getInternalFormat(ctx, _format, _type)
    const type = getType(ctx, _type)

    let destroyed = false
    ctx.textureCount += 1

    return {
        id,
        target,
        format,
        internalFormat,
        type,

        load: (data: TextureImage<any> | TextureVolume<any>) => {
            gl.bindTexture(target, texture)
            // unpack alignment of 1 since we use textures only for data
            gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
            gl.pixelStorei(gl.UNPACK_COLORSPACE_CONVERSION_WEBGL, gl.NONE);
            if (target === gl.TEXTURE_2D) {
                const { array, width, height } = data as TextureImage<any>;
                // TODO remove cast when webgl2 types are fixed
                (gl as WebGLRenderingContext).texImage2D(target, 0, internalFormat, width, height, 0, format, type, array)
            } else if (target === (gl as WebGL2RenderingContext).TEXTURE_3D) {
                const { array, width, height, depth } = data as TextureVolume<any>;
                (gl as WebGL2RenderingContext).texImage3D(target, 0, internalFormat, width, height, depth, 0, format, type, array)
            } else {
                throw new Error('unknown texture target')
            }
            gl.texParameteri(target, gl.TEXTURE_MAG_FILTER, filter)
            gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, filter)
            // clamp-to-edge needed for non-power-of-two textures
            gl.texParameteri(target, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(target, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
            gl.bindTexture(target, null)
        },
        bind: (id: TextureId) => {
            gl.activeTexture(gl.TEXTURE0 + id)
            gl.bindTexture(target, texture)
        },
        unbind: (id: TextureId) => {
            gl.activeTexture(gl.TEXTURE0 + id)
            gl.bindTexture(target, null)
        },
        attachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment) => {
            if (target !== gl.TEXTURE_2D) throw new Error('framebuffer texture must be 2d')
            framebuffer.bind()
            gl.framebufferTexture2D(gl.FRAMEBUFFER, getAttachment(ctx, attachment), gl.TEXTURE_2D, texture, 0)
        },
        destroy: () => {
            if (destroyed) return
            gl.deleteTexture(texture)
            destroyed = true
            ctx.textureCount -= 1
        }
    }
}

export function createTextures(ctx: Context, schema: RenderableSchema, values: TextureValues) {
    const textures: Textures = {}
    Object.keys(schema).forEach((k, i) => {
        const spec = schema[k]
        if (spec.type === 'texture') {
            const texture = createTexture(ctx, spec.kind, spec.format, spec.dataType, spec.filter)
            texture.load(values[k].ref.value)
            textures[k] = texture
        }
    })
    return textures
}