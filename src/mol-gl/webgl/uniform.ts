/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, Mat4, Vec2, Vec3, Vec4 } from 'mol-math/linear-algebra'

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

export type UniformDefs = { [k: string]: UniformKind }
export type UniformValues<T extends UniformDefs> = { [K in keyof T]: UniformKindValue[T[K]] }
export type UniformSetters<T extends UniformDefs> = { [K in keyof T]: (value: UniformKindValue[T[K]]) => void }

export function createUniformSetter<K extends UniformKind, V = UniformKindValue[K]>(gl: WebGLRenderingContext, program: WebGLProgram, name: string, kind: K): (value: V) => void {
    const location = gl.getUniformLocation(program, name)
    switch (kind) {
        case 'f' as K: return (value: V) => gl.uniform1f(location, value as any as number)
        case 'i': case 't2': return (value: V) => gl.uniform1i(location, value as any as number)
        case 'v2': return (value: V) => gl.uniform2fv(location, value as any as Vec2)
        case 'v3': return (value: V) => gl.uniform3fv(location, value as any as Vec3)
        case 'v4': return (value: V) => gl.uniform4fv(location, value as any as Vec4)
        case 'm3': return (value: V) => gl.uniformMatrix3fv(location, false, value as any as Mat3)
        case 'm4': return (value: V) => gl.uniformMatrix4fv(location, false, value as any as Mat4)
    }
    throw new Error('Should never happen')
}

export function getUniformSetters<T extends UniformDefs, K = keyof T>(gl: WebGLRenderingContext, program: WebGLProgram, uniforms: T) {
    const setters: Partial<UniformSetters<T>> = {}
    Object.keys(uniforms).forEach(k => {
        setters[k] = createUniformSetter(gl, program, k, uniforms[k])
    })
    return setters as UniformSetters<T>
}