/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ShaderCode, DefineValues, addShaderDefines } from '../shader-code'
import { WebGLContext } from './context';
import { getUniformUpdaters, getTextureUniformUpdaters, UniformValues } from './uniform';
import { AttributeBuffers } from './buffer';
import { TextureId, Textures } from './texture';
import { createReferenceCache, ReferenceCache } from 'mol-util/reference-cache';
import { idFactory } from 'mol-util/id-factory';
import { RenderableSchema } from '../renderable/schema';
import { hashFnv32a, hashString } from 'mol-data/util';

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

function getAttributeLocations(ctx: WebGLContext, program: WebGLProgram, schema: RenderableSchema) {
    const { gl } = ctx
    const locations: AttributeLocations = {}
    gl.useProgram(program)
    Object.keys(schema).forEach(k => {
        const spec = schema[k]
        if (spec.type === 'attribute') {
            const loc = gl.getAttribLocation(program, k)
            // if (loc === -1) {
            //     console.info(`Could not get attribute location for '${k}'`)
            // }
            locations[k] = loc
        }
    })
    return locations
}

export interface ProgramProps {
    defineValues: DefineValues,
    shaderCode: ShaderCode,
    schema: RenderableSchema
}

export function createProgram(ctx: WebGLContext, props: ProgramProps): Program {
    const { gl, shaderCache } = ctx
    const { defineValues, shaderCode: _shaderCode, schema } = props

    const program = gl.createProgram()
    if (program === null) {
        throw new Error('Could not create WebGL program')
    }

    const shaderCode = addShaderDefines(ctx, defineValues, _shaderCode)
    const vertShaderRef = shaderCache.get(ctx, { type: 'vert', source: shaderCode.vert })
    const fragShaderRef = shaderCache.get(ctx, { type: 'frag', source: shaderCode.frag })

    vertShaderRef.value.attach(program)
    fragShaderRef.value.attach(program)
    gl.linkProgram(program)

    const uniformUpdaters = getUniformUpdaters(ctx, program, schema)
    const attributeLocations = getAttributeLocations(ctx, program, schema)
    const textureUniformUpdaters = getTextureUniformUpdaters(ctx, program, schema)

    let destroyed = false

    return {
        id: getNextProgramId(),

        use: () => {
            Object.keys(uniformUpdaters).forEach(k => uniformUpdaters[k].clear())
            Object.keys(textureUniformUpdaters).forEach(k => textureUniformUpdaters[k].clear())
            gl.useProgram(program)
        },
        setUniforms: (uniformValues: UniformValues) => {
            Object.keys(uniformValues).forEach(k => {
                const uv = uniformValues[k]
                if (uv !== undefined) uniformUpdaters[k].set(uv.ref.value)
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
                textureUniformUpdaters[k].set(i)
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

export type ProgramCache = ReferenceCache<Program, ProgramProps, WebGLContext>

function defineValueHash(v: boolean | number | string): number {
    return typeof v === 'boolean' ? (v ? 1 : 0) :
        typeof v === 'number' ? v : hashString(v)
}

export function createProgramCache(): ProgramCache {
    return createReferenceCache(
        (props: ProgramProps) => {
            const array = [ props.shaderCode.id ]
            Object.keys(props.defineValues).forEach(k => {
                const v = props.defineValues[k].ref.value
                array.push(hashString(k), defineValueHash(v))
            })
            return hashFnv32a(array).toString()
        },
        (ctx: WebGLContext, props: ProgramProps) => createProgram(ctx, props),
        (program: Program) => { program.destroy() }
    )
}