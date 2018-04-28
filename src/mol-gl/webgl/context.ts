/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

interface Reference<T> { usageCount: number, value: T }

export interface Context {
    gl: WebGLRenderingContext
    shaderCache: Map<string, Reference<WebGLShader>>
    extensions: {
        angleInstancedArrays: ANGLE_instanced_arrays
        oesElementIndexUint: OES_element_index_uint
    }
}

export function createContext(gl: WebGLRenderingContext): Context {
    const angleInstancedArrays = gl.getExtension('ANGLE_instanced_arrays')
    if (angleInstancedArrays === null) {
        throw new Error('Could not get "ANGLE_instanced_arrays" extension')
    }
    const oesElementIndexUint = gl.getExtension('OES_element_index_uint')
    if (oesElementIndexUint === null) {
        throw new Error('Could not get "OES_element_index_uint" extension')
    }
    return {
        gl,
        shaderCache: new Map(),
        extensions: { angleInstancedArrays, oesElementIndexUint }
    }
}