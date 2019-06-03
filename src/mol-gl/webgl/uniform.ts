/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, Mat4, Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra'
import { ValueCell } from '../../mol-util';
import { GLRenderingContext } from './compat';
import { RenderableSchema } from '../../mol-gl/renderable/schema';

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

export type UniformValues = { [k: string]: ValueCell<UniformType> }
export type UniformsList = [string, ValueCell<UniformType>][]

export function getUniformType(gl: GLRenderingContext, kind: UniformKind) {
    switch (kind) {
        case 'f': return gl.FLOAT
        case 'i': return gl.INT
        case 'v2': return gl.FLOAT_VEC2
        case 'v3': return gl.FLOAT_VEC3
        case 'v4': return gl.FLOAT_VEC4
        case 'm3': return gl.FLOAT_MAT3
        case 'm4': return gl.FLOAT_MAT4
        default: console.error(`unknown uniform kind '${kind}'`)
    }
}

export function setUniform(gl: GLRenderingContext, location: WebGLUniformLocation | null, kind: UniformKind, value: any) {
    switch (kind) {
        case 'f': gl.uniform1f(location, value); break
        case 'i': case 't': gl.uniform1i(location, value); break
        case 'v2': gl.uniform2fv(location, value); break
        case 'v3': gl.uniform3fv(location, value); break
        case 'v4': gl.uniform4fv(location, value); break
        case 'm3': gl.uniformMatrix3fv(location, false, value); break
        case 'm4': gl.uniformMatrix4fv(location, false, value); break
        default: console.error(`unknown uniform kind '${kind}'`)
    }
}

export type UniformSetter = (gl: GLRenderingContext, location: number, value: any) => void
export type UniformSetters = { [k: string]: UniformSetter }

function uniform1f (gl: GLRenderingContext, location: number, value: any) { gl.uniform1f(location, value) }
function uniform1i (gl: GLRenderingContext, location: number, value: any) { gl.uniform1i(location, value) }
function uniform2fv (gl: GLRenderingContext, location: number, value: any) { gl.uniform2fv(location, value) }
function uniform3fv (gl: GLRenderingContext, location: number, value: any) { gl.uniform3fv(location, value) }
function uniform4fv (gl: GLRenderingContext, location: number, value: any) { gl.uniform4fv(location, value) }
function uniformMatrix3fv (gl: GLRenderingContext, location: number, value: any) { gl.uniformMatrix3fv(location, false, value) }
function uniformMatrix4fv (gl: GLRenderingContext, location: number, value: any) { gl.uniformMatrix4fv(location, false, value) }

function getUniformSetter(kind: UniformKind) {
    switch (kind) {
        case 'f': return uniform1f
        case 'i': case 't': return uniform1i
        case 'v2': return uniform2fv
        case 'v3': return uniform3fv
        case 'v4': return uniform4fv
        case 'm3': return uniformMatrix3fv
        case 'm4': return uniformMatrix4fv
    }
    throw new Error(`unknown uniform kind '${kind}'`)
}

export function getUniformSetters(schema: RenderableSchema) {
    const setters: UniformSetters = {}
    Object.keys(schema).forEach(k => {
        const spec = schema[k]
        if (spec.type === 'uniform') {
            setters[k] = getUniformSetter(spec.kind as UniformKind)
        } else if (spec.type === 'texture') {
            setters[k] = getUniformSetter('t')
        }
    })
    return setters
}