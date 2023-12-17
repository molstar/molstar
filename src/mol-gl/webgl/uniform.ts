/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, Mat4, Vec2, Vec3, Vec4 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { GLRenderingContext } from './compat';
import { RenderableSchema } from '../../mol-gl/renderable/schema';
import { ValueOf } from '../../mol-util/type-helpers';
import { deepClone } from '../../mol-util/object';

export type UniformKindValue = {
    'b': boolean; 'b[]': boolean[]
    'f': number; 'f[]': number[]
    'i': number; 'i[]': number[]
    'v2': Vec2; 'v2[]': number[]
    'v3': Vec3; 'v3[]': number[]
    'v4': Vec4; 'v4[]': number[]
    'iv2': Vec2; 'iv2[]': number[]
    'iv3': Vec3; 'iv3[]': number[]
    'iv4': Vec4; 'iv4[]': number[]
    'm3': Mat3; 'm3[]': number[]
    'm4': Mat4; 'm4[]': number[]
    't': number; 't[]': number[]
}
export type UniformKind = keyof UniformKindValue
export type UniformType = ValueOf<UniformKindValue>

export type UniformValues = { [k: string]: ValueCell<UniformType> }
export type UniformsList = [string, ValueCell<UniformType>][]

export function getUniformType(gl: GLRenderingContext, kind: UniformKind) {
    switch (kind) {
        case 'b': case 'b[]': return gl.BOOL;
        case 'f': case 'f[]': return gl.FLOAT;
        case 'i': case 'i[]': return gl.INT;
        case 'v2': case 'v2[]': return gl.FLOAT_VEC2;
        case 'v3': case 'v3[]': return gl.FLOAT_VEC3;
        case 'v4': case 'v4[]': return gl.FLOAT_VEC4;
        case 'iv2': case 'iv2[]': return gl.INT_VEC2;
        case 'iv3': case 'iv3[]': return gl.INT_VEC3;
        case 'iv4': case 'iv4[]': return gl.INT_VEC4;
        case 'm3': case 'm3[]': return gl.FLOAT_MAT3;
        case 'm4': case 'm4[]': return gl.FLOAT_MAT4;
        default: console.error(`unknown uniform kind '${kind}'`);
    }
}

export function isArrayUniform(kind: UniformKind) {
    return kind.endsWith('[]');
}

export type UniformSetter = (gl: GLRenderingContext, location: number, value: any) => void
export type UniformSetters = { [k: string]: UniformSetter }

function uniform1f(gl: GLRenderingContext, location: number, value: any) { gl.uniform1f(location, value); }
function uniform1fv(gl: GLRenderingContext, location: number, value: any) { gl.uniform1fv(location, value); }
function uniform1i(gl: GLRenderingContext, location: number, value: any) { gl.uniform1i(location, value); }
function uniform1iv(gl: GLRenderingContext, location: number, value: any) { gl.uniform1iv(location, value); }
function uniform2fv(gl: GLRenderingContext, location: number, value: any) { gl.uniform2fv(location, value); }
function uniform3fv(gl: GLRenderingContext, location: number, value: any) { gl.uniform3fv(location, value); }
function uniform4fv(gl: GLRenderingContext, location: number, value: any) { gl.uniform4fv(location, value); }
function uniform2iv(gl: GLRenderingContext, location: number, value: any) { gl.uniform2iv(location, value); }
function uniform3iv(gl: GLRenderingContext, location: number, value: any) { gl.uniform3iv(location, value); }
function uniform4iv(gl: GLRenderingContext, location: number, value: any) { gl.uniform4iv(location, value); }
function uniformMatrix3fv(gl: GLRenderingContext, location: number, value: any) { gl.uniformMatrix3fv(location, false, value); }
function uniformMatrix4fv(gl: GLRenderingContext, location: number, value: any) { gl.uniformMatrix4fv(location, false, value); }

function getUniformSetter(kind: UniformKind): UniformSetter {
    switch (kind) {
        case 'f': return uniform1f;
        case 'f[]': return uniform1fv;
        case 'i': case 't': case 'b': return uniform1i;
        case 'i[]': case 't[]': case 'b[]': return uniform1iv;
        case 'v2': case 'v2[]': return uniform2fv;
        case 'v3': case 'v3[]': return uniform3fv;
        case 'v4': case 'v4[]': return uniform4fv;
        case 'iv2': case 'iv2[]': return uniform2iv;
        case 'iv3': case 'iv3[]': return uniform3iv;
        case 'iv4': case 'iv4[]': return uniform4iv;
        case 'm3': case 'm3[]': return uniformMatrix3fv;
        case 'm4': case 'm4[]': return uniformMatrix4fv;
    }
}

export function getUniformSetters(schema: RenderableSchema) {
    const setters: UniformSetters = {};
    Object.keys(schema).forEach(k => {
        const spec = schema[k];
        if (spec.type === 'uniform') {
            setters[k] = getUniformSetter(spec.kind);
        } else if (spec.type === 'texture') {
            setters[k] = getUniformSetter('t');
        }
    });
    return setters;
}

export function getUniformGlslType(kind: UniformKind): string {
    switch (kind) {
        case 'f': return 'float';
        case 'i': return 'int';
        case 't': return 'sampler2D';
        case 'b': return 'bool';
        case 'v2': return 'vec2';
        case 'v3': return 'vec3';
        case 'v4': return 'vec4';
        case 'm3': return 'mat3';
        case 'm4': return 'mat4';
    }
    throw new Error(`${kind} has no primitive GLSL type.`);
}

export function isUniformValueScalar(kind: UniformKind): boolean {
    switch (kind) {
        case 'f':
        case 'i':
        case 'b':
            return true;
        default:
            return false;
    }
}

export function cloneUniformValues(uniformValues: UniformValues): UniformValues {
    const clonedValues: UniformValues = {};
    Object.keys(uniformValues).forEach(k => {
        clonedValues[k] = ValueCell.create(deepClone(uniformValues[k].ref.value));
    });
    return clonedValues;
}