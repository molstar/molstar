/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context'
import { TextureImage } from '../renderable/util';
import { ValueCell } from 'mol-util';
import { RenderableSchema } from '../renderable/schema';
import { idFactory } from 'mol-util/id-factory';
import { Framebuffer } from './framebuffer';

const getNextTextureId = idFactory()

export type TextureKindValue = {
    'image-uint8': TextureImage<Uint8Array>
    'image-float32': TextureImage<Float32Array>
}
export type TextureKind = keyof TextureKindValue
export type TextureType = 'ubyte' | 'float'
export type TextureFormat = 'alpha' | 'rgb' | 'rgba'
export type TextureAttachment = 'depth' | 'stencil' | 'color0'
export type TextureFilter = 'nearest' | 'linear'

export function getFormat(ctx: Context, format: TextureFormat) {
    const { gl } = ctx
    switch (format) {
        case 'alpha': return gl.ALPHA
        case 'rgb': return gl.RGB
        case 'rgba': return gl.RGBA
    }
}

export function getInternalFormat(ctx: Context, format: TextureFormat, type: TextureType) {
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

export function getType(ctx: Context, type: TextureType) {
    const { gl } = ctx
    switch (type) {
        case 'ubyte': return gl.UNSIGNED_BYTE
        case 'float': return gl.FLOAT
    }
}

export function getFilter(ctx: Context, type: TextureFilter) {
    const { gl } = ctx
    switch (type) {
        case 'nearest': return gl.NEAREST
        case 'linear': return gl.LINEAR
    }
}

export function getAttachment(ctx: Context, attachment: TextureAttachment) {
    const { gl } = ctx
    switch (attachment) {
        case 'depth': return gl.DEPTH_ATTACHMENT
        case 'stencil': return gl.STENCIL_ATTACHMENT
        case 'color0': return gl.COLOR_ATTACHMENT0
    }
}

export interface Texture {
    readonly id: number
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

export function createTexture(ctx: Context, _format: TextureFormat, _type: TextureType, _filter: TextureFilter): Texture {
    const id = getNextTextureId()
    const { gl } = ctx
    const texture = gl.createTexture()
    if (texture === null) {
        throw new Error('Could not create WebGL texture')
    }



    const filter = getFilter(ctx, _filter)
    const format = getFormat(ctx, _format)
    const internalFormat = getInternalFormat(ctx, _format, _type)
    const type = getType(ctx, _type)

    let destroyed = false
    ctx.textureCount += 1

    return {
        id,
        format,
        internalFormat,
        type,

        load: (image: TextureImage<any>) => {
            const { array, width, height } = image
            gl.bindTexture(gl.TEXTURE_2D, texture)
            // unpack alignment of 1 since we use textures only for data
            gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
            gl.pixelStorei(gl.UNPACK_COLORSPACE_CONVERSION_WEBGL, gl.NONE);
            (gl as WebGLRenderingContext).texImage2D(gl.TEXTURE_2D, 0, internalFormat, width, height, 0, format, type, array) // TODO remove cast when webgl2 types are fixed
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter)
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filter)
            // clamp-to-edge needed for non-power-of-two textures
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
            gl.bindTexture(gl.TEXTURE_2D, null)
        },
        bind: (id: TextureId) => {
            gl.activeTexture(gl.TEXTURE0 + id)
            gl.bindTexture(gl.TEXTURE_2D, texture)
        },
        unbind: (id: TextureId) => {
            gl.activeTexture(gl.TEXTURE0 + id)
            gl.bindTexture(gl.TEXTURE_2D, null)
        },
        attachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment) => {
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
            const texture = createTexture(ctx, spec.format, spec.dataType, spec.filter)
            texture.load(values[k].ref.value)
            textures[k] = texture
        }
    })
    return textures
}