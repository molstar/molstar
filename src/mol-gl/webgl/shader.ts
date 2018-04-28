/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Context } from './context'

function addLineNumbers(source: string) {
    const lines = source.split('\n')
    for (let i = 0; i < lines.length; ++i) {
        lines[i] = (i + 1) + ': ' + lines[i]
    }
    return lines.join('\n')
}

type ShaderType = 'vert' | 'frag'

function createShader(gl: WebGLRenderingContext, type: ShaderType, source: string) {
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

    return shader
}

export function getShader(ctx: Context, type: ShaderType, source: string) {
    let shaderRef = ctx.shaderCache.get(source)
    if (!shaderRef) {
        shaderRef = { usageCount: 0, value: createShader(ctx.gl, type, source) }
        ctx.shaderCache.set(source, shaderRef)
    }
    shaderRef.usageCount += 1
    return {
        free: () => {
            if (shaderRef) {
                shaderRef.usageCount -= 1
                shaderRef = undefined
            }
        },
        value: shaderRef.value
    }
}