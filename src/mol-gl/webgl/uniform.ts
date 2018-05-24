/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, Mat4, Vec2, Vec3, Vec4 } from 'mol-math/linear-algebra'
import { Context } from './context';
import { TextureImage } from '../renderable/util';
import { ValueCell } from 'mol-util';

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
export interface UniformUpdater {
    set: (value: UniformType, version: number) => void,
    clear: () => void
}

export type UniformDefs = { [k: string]: UniformKind }
export type UniformValues = { [k: string]: ValueCell<UniformType> }
export type UniformUpdaters = { [k: string]: UniformUpdater }

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

export function getUniformUpdaters(ctx: Context, program: WebGLProgram, uniforms: UniformDefs) {
    const updaters: UniformUpdaters = {}
    Object.keys(uniforms).forEach(k => {
        const setter = createUniformSetter(ctx, program, k, uniforms[k])
        let _version = -1
        updaters[k] = {
            set: (value, version) => {
                if (_version !== version) {
                    setter(value)
                    _version = version
                }
            },
            clear: () => {
                _version = -1
            }
        }
    })
    return updaters
}