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

const getNextTextureId = idFactory()

export type TextureFormat = 'rgb' | 'rgba'
export type TextureType = 'ubyte' | 'uint'

export function getFormat(ctx: Context, format: TextureFormat) {
    const { gl } = ctx
    switch (format) {
        case 'rgb': return gl.RGB
        case 'rgba': return gl.RGBA
    }
}

export function getType(ctx: Context, type: TextureType) {
    const { gl } = ctx
    switch (type) {
        case 'ubyte': return gl.UNSIGNED_BYTE
        case 'uint': return gl.UNSIGNED_INT
    }
}

export interface Texture {
    readonly id: number
    readonly texture: WebGLTexture
    readonly format: number
    readonly type: number

    readonly width: number
    readonly height: number

    load: (image: TextureImage) => void
    bind: (id: TextureId) => void
    unbind: (id: TextureId) => void
    destroy: () => void
}

export type TextureId = 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15

export type TextureValues = { [k: string]: ValueCell<TextureImage> }
export type Textures = { [k: string]: Texture }

export function createTexture(ctx: Context, _format: TextureFormat, _type: TextureType): Texture {
    const id = getNextTextureId()
    const { gl } = ctx
    const texture = gl.createTexture()
    if (texture === null) {
        throw new Error('Could not create WebGL texture')
    }

    const magFilter = gl.NEAREST
    const minFilter = gl.NEAREST
    const format = getFormat(ctx, _format)
    const type = getType(ctx, _type)

    let _width = 0
    let _height = 0

    let destroyed = false
    ctx.textureCount += 1

    return {
        id,
        texture,
        format,
        type,

        get width () { return _width },
        get height () { return _height },

        load: (image: TextureImage) => {
            const { array, width, height } = image
            gl.bindTexture(gl.TEXTURE_2D, texture)
            // unpack alignment of 1 since we use textures only for data
            gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
            gl.texImage2D(gl.TEXTURE_2D, 0, format, width, height, 0, format, type, array)
            _width = width
            _height = height
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, magFilter)
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, minFilter)
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
            const texture = createTexture(ctx, spec.format, spec.dataType)
            texture.load(values[k].ref.value)
            textures[k] = texture
        }
    })
    return textures
}