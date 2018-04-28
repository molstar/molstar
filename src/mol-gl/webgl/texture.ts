/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// import { Context } from './context'

export interface Texture {
    load: (image: ImageData) => void
    bind: (id: TextureId) => void
}

export type TextureId = 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15
export type TextureTarget = 'TEXTURE0' | 'TEXTURE1' | 'TEXTURE2' | 'TEXTURE3' | 'TEXTURE4' | 'TEXTURE5' | 'TEXTURE6' | 'TEXTURE7' | 'TEXTURE8' | 'TEXTURE9' | 'TEXTURE10' | 'TEXTURE11' | 'TEXTURE12' | 'TEXTURE13' | 'TEXTURE14' | 'TEXTURE15'

export type TextureDefs = { [k: string]: '' }
export type TextureUniforms<T extends TextureDefs> = { [k in keyof T]: 't2' }
export type TextureValues<T extends TextureDefs> = { [k in keyof T]: ImageData }
export type Textures<T extends TextureDefs> = { [k in keyof T]: Texture }

export function createTexture(gl: WebGLRenderingContext): Texture {
    const texture = gl.createTexture()
    if (texture === null) {
        throw new Error('Could not create WebGL texture')
    }

    return {
        load: (image: ImageData) => {
            gl.bindTexture(gl.TEXTURE_2D, texture)
            gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true)
            gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, image)
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST)
            gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST)
            gl.bindTexture(gl.TEXTURE_2D, null)
        },
        bind: (id: TextureId) => {
            gl.activeTexture(gl[`TEXTURE${id}` as TextureTarget])
            gl.bindTexture(gl.TEXTURE_2D, texture)
        }
    }
}

export function createTextures<T extends TextureDefs>(gl: WebGLRenderingContext, props: T, state: TextureValues<T>) {
    const textures: Partial<Textures<T>> = {}
    Object.keys(props).forEach(k => {
        const texture = createTexture(gl)
        texture.load(state[k])
        textures[k] = texture
    })
    return textures as Textures<T>
}
