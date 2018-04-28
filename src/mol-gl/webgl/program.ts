/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Shaders } from '../shaders'
import { getShader } from './shader'
import { Context } from './context';
import { getUniformSetters, UniformDefs, UniformValues } from './uniform';
import {AttributeDefs, AttributeBuffers } from './buffer';
import { TextureId, TextureDefs, TextureUniforms, Textures } from './texture';

export interface Program<U extends UniformDefs, A extends AttributeDefs, T extends TextureDefs> {
    readonly id: number

    setUniforms: (uniformValues: Partial<UniformValues<U>>) => void
    bindAttributes: (attribueBuffers: AttributeBuffers<A>) => void
    bindTextures: (textures: Textures<T>) => void

    destroy: () => void
}

type AttributeLocations<T extends AttributeDefs> = { [K in keyof T]: number }

function getAttributeLocations<A extends AttributeDefs>(gl: WebGLRenderingContext, program: WebGLProgram, attributes: A) {
    gl.useProgram(program)
    const locations: Partial<AttributeLocations<A>> = {}
    Object.keys(attributes).forEach(k => {
        const loc = gl.getAttribLocation(program, k)
        gl.enableVertexAttribArray(loc)
        locations[k] = loc
    })
    return locations as AttributeLocations<A>
}

function getTextureUniforms<T extends TextureDefs>(textures: T) {
    const textureUniforms: Partial<TextureUniforms<T>> = {}
    Object.keys(textureUniforms).forEach(k => textureUniforms[k] = 't2')
    return textureUniforms as TextureUniforms<T>
}

export function createProgram<U extends UniformDefs, A extends AttributeDefs, T extends TextureDefs>(ctx: Context, shaders: Shaders, uniformDefs: U, attributeDefs: A, textureDefs: T): Program<U, A, T> {
    const { gl } = ctx

    const program = gl.createProgram()
    if (program === null) {
        throw new Error('Could not create WebGL program')
    }

    const glVertShader = getShader(ctx, 'vert', shaders.vert)
    const glFragShader = getShader(ctx, 'frag', shaders.frag)

    gl.attachShader(program, glVertShader.value)
    gl.attachShader(program, glFragShader.value)
    gl.linkProgram(program)

    const uniformSetters = getUniformSetters(gl, program, uniformDefs)
    const attributeLocations = getAttributeLocations(gl, program, attributeDefs)
    const textureUniforms = getTextureUniforms(textureDefs)
    const textureUniformSetters = getUniformSetters(gl, program, textureUniforms)

    let destroyed = false

    return {
        id: 0,

        setUniforms: (uniformValues: Partial<UniformValues<U>>) => {
            Object.keys(uniformValues).forEach(k => {
                const value = uniformValues[k]
                if (value !== undefined) uniformSetters[k](value)
            })
        },
        bindAttributes: (attribueBuffers: AttributeBuffers<A>) => {
            Object.keys(attribueBuffers).forEach(k => {
                attribueBuffers[k].bind(attributeLocations[k], 0, 0)
            })
        },
        bindTextures: (textures: Textures<T>) => {
            Object.keys(textures).forEach((k, i) => {
                textures[k].bind(i as TextureId)
                textureUniformSetters[k](i)
            })
        },

        destroy: () => {
            if (destroyed) return
            glVertShader.free()
            glFragShader.free()
            gl.deleteProgram(program)
            destroyed = true
        }
    }
}