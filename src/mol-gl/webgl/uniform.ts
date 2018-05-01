/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, Mat4, Vec2, Vec3, Vec4 } from 'mol-math/linear-algebra'
import { Context } from './context';
import { TextureImage } from '../renderable/util';

export type UniformKindValue = {
    'f': number
    'i': number
    'v2': Vec2
    'v3': Vec3
    'v4': Vec4
    'm3': Mat3
    'm4': Mat4
    't2': number
}
export type UniformKind = keyof UniformKindValue
export type UniformType = number | Vec2 | Vec3 | Vec4 | Mat3 | Mat4 | TextureImage

export type UniformDefs = { [k: string]: UniformKind }
export type UniformValues = { [k: string]: UniformType }
export type UniformSetters = { [k: string]: (value: UniformType) => void }

export function createUniformSetter(ctx: Context, program: WebGLProgram, name: string, kind: UniformKind): (value: any) => void {
    const { gl } = ctx
    const location = gl.getUniformLocation(program, name)
    if (location === null) {
        console.info(`Could not get WebGL uniform location for '${name}'`)
    }
    switch (kind) {
        case 'f': return (value: number) => gl.uniform1f(location, value)
        case 'i': case 't2': return (value: number) => gl.uniform1i(location, value)
        case 'v2': return (value: Vec2) => gl.uniform2fv(location, value)
        case 'v3': return (value: Vec3) => gl.uniform3fv(location, value)
        case 'v4': return (value: Vec4) => gl.uniform4fv(location, value)
        case 'm3': return (value: Mat3) => gl.uniformMatrix3fv(location, false, value)
        case 'm4': return (value: Mat4) => gl.uniformMatrix4fv(location, false, value)
    }
}

export function getUniformSetters(ctx: Context, program: WebGLProgram, uniforms: UniformDefs) {
    const setters: UniformSetters = {}
    Object.keys(uniforms).forEach(k => {
        setters[k] = createUniformSetter(ctx, program, k, uniforms[k])
    })
    return setters
}