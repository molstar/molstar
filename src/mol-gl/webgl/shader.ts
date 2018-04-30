/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createReferenceCache, ReferenceCache } from 'mol-util/reference-cache';
import { Context } from './context';

function addLineNumbers(source: string) {
    const lines = source.split('\n')
    for (let i = 0; i < lines.length; ++i) {
        lines[i] = (i + 1) + ': ' + lines[i]
    }
    return lines.join('\n')
}

export type ShaderType = 'vert' | 'frag'
export interface ShaderProps { type: ShaderType, source: string }
export interface Shader {
    attach: (program: WebGLProgram) => void
    destroy: () => void
}

function createShader(ctx: Context, props: ShaderProps): Shader {
    const { gl } = ctx
    const { type, source } = props

    const shader = gl.createShader(type === 'vert' ? gl.VERTEX_SHADER : gl.FRAGMENT_SHADER)
    if (shader === null) {
        throw new Error(`Error creating ${type} shader`)
    }

    gl.shaderSource(shader, source)
    gl.compileShader(shader)

    if (gl.getShaderParameter(shader, gl.COMPILE_STATUS) === false) {
        console.warn(`'${type}' shader info log '${gl.getShaderInfoLog(shader)}'\n${addLineNumbers(source)}`)
        throw new Error(`Error compiling ${type} shader`)
    }

    return {
        attach: (program: WebGLProgram) => {
            gl.attachShader(program, shader)
        },
        destroy: () => {
            gl.deleteShader(shader)
        }
    }
}

export type ShaderCache = ReferenceCache<Shader, ShaderProps, Context>

export function createShaderCache(): ShaderCache {
    return createReferenceCache(
        (props: ShaderProps) => JSON.stringify(props),
        (ctx: Context, props: ShaderProps) => createShader(ctx, props),
        (shader: Shader) => { shader.destroy() }
    )
}