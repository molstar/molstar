/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, Mat4, Vec2, Vec3, Vec4 } from 'mol-math/linear-algebra'
import { WebGLContext } from './context';
import { ValueCell, arrayEqual } from 'mol-util';
import { RenderableSchema } from '../renderable/schema';

export type UniformKindValue = {
    'f': number
    'i': number
    'v2': Vec2
    'v3': Vec3
    'v4': Vec4
    'm3': Mat3
    'm4': Mat4
    't': number
}
export type UniformKind = keyof UniformKindValue
export type UniformType = number | Vec2 | Vec3 | Vec4 | Mat3 | Mat4
export interface UniformUpdater {
    set: (value: UniformType) => void,
    clear: () => void
}

export type UniformValues = { [k: string]: ValueCell<UniformType> }
export type UniformUpdaters = { [k: string]: UniformUpdater }

function createUniformSetter(ctx: WebGLContext, program: WebGLProgram, name: string, kind: UniformKind): (value: any) => void {
    const { gl } = ctx
    const location = gl.getUniformLocation(program, name)
    if (location === null) {
        // console.info(`Could not get WebGL uniform location for '${name}'`)
    }
    switch (kind) {
        case 'f': return (value: number) => gl.uniform1f(location, value)
        case 'i': case 't': return (value: number) => gl.uniform1i(location, value)
        case 'v2': return (value: Vec2) => (gl as WebGLRenderingContext).uniform2fv(location, value) // TODO remove cast when webgl2 types are fixed
        case 'v3': return (value: Vec3) => (gl as WebGLRenderingContext).uniform3fv(location, value)
        case 'v4': return (value: Vec4) => (gl as WebGLRenderingContext).uniform4fv(location, value)
        case 'm3': return (value: Mat3) => (gl as WebGLRenderingContext).uniformMatrix3fv(location, false, value)
        case 'm4': return (value: Mat4) => (gl as WebGLRenderingContext).uniformMatrix4fv(location, false, value)
    }
}

function createUniformUpdater(ctx: WebGLContext, program: WebGLProgram, name: string, kind: UniformKind): UniformUpdater {
    const setter = createUniformSetter(ctx, program, name, kind)
    let _value: UniformType | undefined = undefined
    return {
        set: value => {
            if (_value !== value || (Array.isArray(_value) && Array.isArray(value) && arrayEqual(_value, value))) {
                setter(value)
                _value = value 
            }
        },
        clear: () => { _value = undefined }
    }
}

export function getUniformUpdaters(ctx: WebGLContext, program: WebGLProgram, schema: RenderableSchema) {
    const updaters: UniformUpdaters = {}
    Object.keys(schema).forEach(k => {
        const spec = schema[k]
        if (spec.type === 'uniform') {
            updaters[k] = createUniformUpdater(ctx, program, k, spec.kind)
        }
    })
    return updaters
}

export function getTextureUniformUpdaters(ctx: WebGLContext, program: WebGLProgram, schema: RenderableSchema) {
    const updaters: UniformUpdaters = {}
    Object.keys(schema).forEach(k => {
        const spec = schema[k]
        if (spec.type === 'texture') {
            updaters[k] = createUniformUpdater(ctx, program, k, 't')
        }
    })
    return updaters
}