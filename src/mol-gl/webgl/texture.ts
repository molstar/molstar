/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context'
import { TextureImage } from '../renderable/util';

export interface Texture {
    load: (image: TextureImage) => void
    bind: (id: TextureId) => void
    unbind: (id: TextureId) => void
    destroy: () => void
}

export type TextureId = 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15

export type TextureDefs = { [k: string]: true }
export type TextureUniformDefs = { [k: string]: 't2' }
export type TextureValues = { [k: string]: TextureImage }
export type Textures = { [k: string]: Texture }

export function createTexture(ctx: Context): Texture {
    const { gl } = ctx
    const texture = gl.createTexture()
    if (texture === null) {
        throw new Error('Could not create WebGL texture')
    }

    const _textureType = gl.TEXTURE_2D
    const _magFilter = gl.NEAREST
    const _minFilter = gl.NEAREST
    const _format = gl.RGB
    const _arrayType = gl.UNSIGNED_BYTE

    let destroyed = false
    ctx.textureCount += 1

    return {
        load: (image: TextureImage) => {
            const { array, width, height } = image
            gl.bindTexture(_textureType, texture)
            // unpack alignment of 1 since we use textures only for data
            gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
            gl.texImage2D(_textureType, 0, _format, width, height, 0, _format, _arrayType, array)
            gl.texParameteri(_textureType, gl.TEXTURE_MAG_FILTER, _magFilter)
            gl.texParameteri(_textureType, gl.TEXTURE_MIN_FILTER, _minFilter)
            // clamp-to-edge needed for non-power-of-two textures
            gl.texParameteri(_textureType, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
            gl.texParameteri(_textureType, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
            gl.bindTexture(_textureType, null)
        },
        bind: (id: TextureId) => {
            gl.activeTexture(gl.TEXTURE0 + id)
            gl.bindTexture(_textureType, texture)
        },
        unbind: (id: TextureId) => {
            gl.activeTexture(gl.TEXTURE0 + id)
            gl.bindTexture(_textureType, null)
        },

        destroy: () => {
            if (destroyed) return
            gl.deleteTexture(texture)
            destroyed = true
            ctx.textureCount -= 1
        }
    }
}

export function createTextures(ctx: Context, props: TextureDefs, state: TextureValues) {
    const textures: Textures = {}
    Object.keys(props).forEach((k, i) => {
        const texture = createTexture(ctx)
        texture.load(state[k])
        textures[k] = texture
    })
    return textures
}