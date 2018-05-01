/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ShaderCode } from '../shader-code'
import { Context } from './context';
import { getUniformSetters, UniformDefs, UniformValues } from './uniform';
import {AttributeDefs, AttributeBuffers } from './buffer';
import { TextureId, TextureDefs, TextureUniformDefs, Textures } from './texture';
import { createReferenceCache, ReferenceCache } from 'mol-util/reference-cache';
import { idFactory } from 'mol-util/id-factory';

const getNextProgramId = idFactory()

export interface Program {
    readonly id: number

    use: () => void
    setUniforms: (uniformValues: UniformValues) => void
    bindAttributes: (attribueBuffers: AttributeBuffers) => void
    bindTextures: (textures: Textures) => void

    destroy: () => void
}

type AttributeLocations = { [k: string]: number }

function getAttributeLocations(ctx: Context, program: WebGLProgram, attributeDefs: AttributeDefs) {
    const { gl } = ctx
    const locations: AttributeLocations = {}
    gl.useProgram(program)
    Object.keys(attributeDefs).forEach(k => {
        const loc = gl.getAttribLocation(program, k)
        if (loc === -1) {
            console.info(`Could not get attribute location for '${k}'`)
        }
        locations[k] = loc
    })
    return locations
}

function getTextureUniformDefs(textureDefs: TextureDefs) {
    const textureUniformDefs: TextureUniformDefs = {}
    Object.keys(textureDefs).forEach(k => textureUniformDefs[k] = 't2')
    return textureUniformDefs
}

export interface ProgramProps {
    shaderCode: ShaderCode,
    uniformDefs: UniformDefs,
    attributeDefs: AttributeDefs,
    textureDefs: TextureDefs
}

export function createProgram(ctx: Context, props: ProgramProps): Program {
    const { gl, shaderCache } = ctx
    const { shaderCode, uniformDefs, attributeDefs, textureDefs } = props

    const program = gl.createProgram()
    if (program === null) {
        throw new Error('Could not create WebGL program')
    }

    const vertShaderRef = shaderCache.get(ctx, { type: 'vert', source: shaderCode.vert })
    const fragShaderRef = shaderCache.get(ctx, { type: 'frag', source: shaderCode.frag })

    vertShaderRef.value.attach(program)
    fragShaderRef.value.attach(program)
    gl.linkProgram(program)

    const uniformSetters = getUniformSetters(ctx, program, uniformDefs)
    const attributeLocations = getAttributeLocations(ctx, program, attributeDefs)
    const textureUniformDefs = getTextureUniformDefs(textureDefs)
    const textureUniformSetters = getUniformSetters(ctx, program, textureUniformDefs)

    let destroyed = false

    return {
        id: getNextProgramId(),

        use: () => {
            gl.useProgram(program)
        },
        setUniforms: (uniformValues: UniformValues) => {
            Object.keys(uniformValues).forEach(k => {
                const value = uniformValues[k]
                if (value !== undefined) uniformSetters[k](value)
            })
        },
        bindAttributes: (attribueBuffers: AttributeBuffers) => {
            Object.keys(attribueBuffers).forEach(k => {
                const loc = attributeLocations[k]
                if (loc !== -1) attribueBuffers[k].bind(loc)
            })
        },
        bindTextures: (textures: Textures) => {
            Object.keys(textures).forEach((k, i) => {
                textures[k].bind(i as TextureId)
                textureUniformSetters[k](i)
            })
        },

        destroy: () => {
            if (destroyed) return
            vertShaderRef.free()
            fragShaderRef.free()
            gl.deleteProgram(program)
            destroyed = true
        }
    }
}

export type ProgramCache = ReferenceCache<Program, ProgramProps, Context>

export function createProgramCache(): ProgramCache {
    return createReferenceCache(
        (props: ProgramProps) => JSON.stringify(props),
        (ctx: Context, props: ProgramProps) => createProgram(ctx, props),
        (program: Program) => { program.destroy() }
    )
}