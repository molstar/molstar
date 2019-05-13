/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ShaderCode, DefineValues, addShaderDefines } from '../shader-code'
import { WebGLState } from './state';
import { WebGLExtensions } from './extensions';
import { getUniformSetters, UniformsList, getUniformType } from './uniform';
import { AttributeBuffers, getAttribType } from './buffer';
import { TextureId, Textures } from './texture';
import { createReferenceCache, ReferenceCache } from 'mol-util/reference-cache';
import { idFactory } from 'mol-util/id-factory';
import { RenderableSchema } from '../renderable/schema';
import { hashFnv32a, hashString } from 'mol-data/util';
import { isDebugMode } from 'mol-util/debug';
import { GLRenderingContext } from './compat';
import { ShaderCache } from './shader';

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

function getLocations(gl: GLRenderingContext, program: WebGLProgram, schema: RenderableSchema) {
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

function checkActiveAttributes(gl: GLRenderingContext, program: WebGLProgram, schema: RenderableSchema) {
    const attribCount = gl.getProgramParameter(program, gl.ACTIVE_ATTRIBUTES);
    for (let i = 0; i < attribCount; ++i) {
        const info = gl.getActiveAttrib(program, i);
        if (info) {
            const { name, type } = info
            if (name.startsWith('__activeAttribute')) {
                // name assigned by `gl.shim.ts`, ignore for checks
                continue
            }
            const spec = schema[name]
            if (spec === undefined) {
                throw new Error(`missing 'uniform' or 'texture' with name '${name}' in schema`)
            }
            if (spec.type !== 'attribute') {
                throw new Error(`'${name}' must be of type 'attribute' but is '${spec.type}'`)
            }
            const attribType = getAttribType(gl, spec.kind, spec.itemSize)
            if (attribType !== type) {
                throw new Error(`unexpected attribute type for ${name}`)
            }
        }
    }
}

function checkActiveUniforms(gl: GLRenderingContext, program: WebGLProgram, schema: RenderableSchema) {
    const attribCount = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);
    for (let i = 0; i < attribCount; ++i) {
        const info = gl.getActiveUniform(program, i);
        if (info) {
            const { name, type } = info
            if (name.startsWith('__activeUniform')) {
                // name assigned by `gl.shim.ts`, ignore for checks
                continue
            }
            const spec = schema[name]
            if (spec === undefined) {
                throw new Error(`missing 'uniform' or 'texture' with name '${name}' in schema`)
            }
            if (spec.type === 'uniform') {
                const uniformType = getUniformType(gl, spec.kind)
                if (uniformType !== type) {
                    throw new Error(`unexpected uniform type for ${name}`)
                }
            } else if (spec.type === 'texture') {
                if (spec.kind === 'image-float32' || spec.kind === 'image-uint8') {
                    if (type !== gl.SAMPLER_2D) {
                        throw new Error(`unexpected sampler type for '${name}'`)
                    }
                } else if (spec.kind === 'volume-float32' || spec.kind === 'volume-uint8') {
                    if (type !== (gl as WebGL2RenderingContext).SAMPLER_3D) {
                        throw new Error(`unexpected sampler type for '${name}'`)
                    }
                } else {
                    // TODO
                }
            } else {
                throw new Error(`'${name}' must be of type 'uniform' or 'texture' but is '${spec.type}'`)
            }
        }
    }
}

export interface ProgramProps {
    defineValues: DefineValues,
    shaderCode: ShaderCode,
    schema: RenderableSchema
}

export function createProgram(gl: GLRenderingContext, state: WebGLState, extensions: WebGLExtensions, shaderCache: ShaderCache, props: ProgramProps): Program {
    const { defineValues, shaderCode: _shaderCode, schema } = props

    const program = gl.createProgram()
    if (program === null) {
        throw new Error('Could not create WebGL program')
    }
    const programId = getNextProgramId()

    const shaderCode = addShaderDefines(gl, extensions, defineValues, _shaderCode)
    const vertShaderRef = shaderCache.get({ type: 'vert', source: shaderCode.vert })
    const fragShaderRef = shaderCache.get({ type: 'frag', source: shaderCode.frag })

    vertShaderRef.value.attach(program)
    fragShaderRef.value.attach(program)
    gl.linkProgram(program)
    if (isDebugMode) {
        // no-op in FF on Mac, see https://bugzilla.mozilla.org/show_bug.cgi?id=1284425
        // gl.validateProgram(program)
        if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
            throw new Error(`Could not compile WebGL program. \n\n${gl.getProgramInfoLog(program)}`);
        }
    }

    const locations = getLocations(gl, program, schema)
    const uniformSetters = getUniformSetters(schema)

    if (isDebugMode) {
        checkActiveAttributes(gl, program, schema)
        checkActiveUniforms(gl, program, schema)
    }

    let destroyed = false

    return {
        id: programId,

        use: () => {
            // console.log('use', programId)
            state.currentProgramId = programId
            gl.useProgram(program)
        },
        setUniforms: (uniformValues: UniformsList) => {
            for (let i = 0, il = uniformValues.length; i < il; ++i) {
                const [k, v] = uniformValues[i]
                if (v) {
                    const l = locations[k]
                    if (l !== null) uniformSetters[k](gl, l, v.ref.value)
                }
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
                const l = locations[k]
                if (l !== null) {
                    // TODO if the order and count of textures in a material can be made invariant
                    //      bind needs to be called only when the material changes
                    texture.bind(i as TextureId)
                    uniformSetters[k](gl, l, i as TextureId)
                }
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

export type ProgramCache = ReferenceCache<Program, ProgramProps>

function defineValueHash(v: boolean | number | string): number {
    return typeof v === 'boolean' ? (v ? 1 : 0) :
        typeof v === 'number' ? v : hashString(v)
}

export function createProgramCache(gl: GLRenderingContext, state: WebGLState, extensions: WebGLExtensions, shaderCache: ShaderCache): ProgramCache {
    return createReferenceCache(
        (props: ProgramProps) => {
            const array = [ props.shaderCode.id ]
            Object.keys(props.defineValues).forEach(k => {
                const v = props.defineValues[k].ref.value
                array.push(hashString(k), defineValueHash(v))
            })
            return hashFnv32a(array).toString()
        },
        (props: ProgramProps) => createProgram(gl, state, extensions, shaderCache, props),
        (program: Program) => { program.destroy() }
    )
}