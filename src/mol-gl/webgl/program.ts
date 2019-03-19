/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ShaderCode, DefineValues, addShaderDefines } from '../shader-code'
import { WebGLContext } from './context';
import { getUniformSetters, UniformsList } from './uniform';
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
    setUniforms: (uniformValues: UniformsList) => void
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
        setUniforms: (uniformValues: UniformsList) => {
            for (let i = 0, il = uniformValues.length; i < il; ++i) {
                const [k, v] = uniformValues[i]
                if (v) uniformSetters[k](gl, locations[k], v.ref.value)
            }
        },
        bindAttributes: (attribueBuffers: AttributeBuffers) => {
            for (let i = 0, il = attribueBuffers.length; i < il; ++i) {
                const [k, buffer] = attribueBuffers[i]
                const l = locations[k]
                if (l !== -1) buffer.bind(l)
            }
        },
        bindTextures: (textures: Textures) => {
            for (let i = 0, il = textures.length; i < il; ++i) {
                const [k, texture] = textures[i]
                texture.bind(i as TextureId)
                // TODO if the order and count of textures in a material can be made invariant
                //      this needs to be set only when the material changes
                uniformSetters[k](gl, locations[k], i as TextureId)
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