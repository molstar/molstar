/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ShaderCode, DefineValues, addShaderDefines } from '../shader-code'
import { WebGLContext } from './context';
import { UniformValues, getUniformSetters } from './uniform';
import { AttributeBuffers } from './buffer';
import { Textures, TextureId } from './texture';
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

type Locations = { [k: string]: number }

function getLocations(ctx: WebGLContext, program: WebGLProgram, schema: RenderableSchema) {
    const { gl } = ctx
    const locations: Locations = {}
    Object.keys(schema).forEach(k => {
        const spec = schema[k]
        if (spec.type === 'attribute') {
            const loc = gl.getAttribLocation(program, k)
            // if (loc === -1) console.info(`Could not get attribute location for '${k}'`)
            locations[k] = loc
        } else if (spec.type === 'uniform' || spec.type === 'texture') {
            const loc = gl.getUniformLocation(program, k)
            // if (loc === null) console.info(`Could not get uniform location for '${k}'`)
            locations[k] = loc as number
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
    const programId = getNextProgramId()

    const shaderCode = addShaderDefines(ctx, defineValues, _shaderCode)
    const vertShaderRef = shaderCache.get(ctx, { type: 'vert', source: shaderCode.vert })
    const fragShaderRef = shaderCache.get(ctx, { type: 'frag', source: shaderCode.frag })

    vertShaderRef.value.attach(program)
    fragShaderRef.value.attach(program)
    gl.linkProgram(program)
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        throw new Error(`Could not compile WebGL program. \n\n${gl.getProgramInfoLog(program)}`);
    }

    const locations = getLocations(ctx, program, schema)
    const uniformSetters = getUniformSetters(schema)

    let destroyed = false

    return {
        id: programId,

        use: () => {
            // console.log('use', programId)
            ctx.currentProgramId = programId
            gl.useProgram(program)
        },
        setUniforms: (uniformValues: UniformValues) => {
            const uniformKeys = Object.keys(uniformValues)
            for (let i = 0, il = uniformKeys.length; i < il; ++i) {
                const k = uniformKeys[i]
                const l = locations[k]
                const v = uniformValues[k]
                if (v) uniformSetters[k](gl, l, v.ref.value)
            }
        },
        bindAttributes: (attribueBuffers: AttributeBuffers) => {
            const attributeKeys = Object.keys(attribueBuffers)
            for (let i = 0, il = attributeKeys.length; i < il; ++i) {
                const k = attributeKeys[i]
                const l = locations[k]
                if (l !== -1) attribueBuffers[k].bind(l)
            }
        },
        bindTextures: (textures: Textures) => {
            const textureKeys = Object.keys(textures)
            for (let i = 0, il = textureKeys.length; i < il; ++i) {
                const k = textureKeys[i]
                const l = locations[k]
                textures[k].bind(i as TextureId)
                uniformSetters[k](gl, l, i as TextureId)
            }
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