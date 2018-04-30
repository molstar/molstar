/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createProgramCache, ProgramCache } from './program'
import { createShaderCache, ShaderCache } from './shader'

// const extensions = [
//     'OES_element_index_uint',
//     'ANGLE_instanced_arrays'
// ]
// const optionalExtensions = [
//     'EXT_disjoint_timer_query'
// ]

function unbindResources (gl: WebGLRenderingContext) {
    // bind null to all texture units
    const maxTextureImageUnits = gl.getParameter(gl.MAX_TEXTURE_IMAGE_UNITS)
    for (let i = 0; i < maxTextureImageUnits; ++i) {
        gl.activeTexture(gl.TEXTURE0 + i)
        gl.bindTexture(gl.TEXTURE_2D, null)
        gl.bindTexture(gl.TEXTURE_CUBE_MAP, null)
    }

    // assign the smallest possible buffer to all attributes
    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    const maxVertexAttribs = gl.getParameter(gl.MAX_VERTEX_ATTRIBS);
    for (let i = 0; i < maxVertexAttribs; ++i) {
        gl.vertexAttribPointer(i, 1, gl.FLOAT, false, 0, 0);
    }

    // bind null to all buffers
    gl.bindBuffer(gl.ARRAY_BUFFER, null)
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null)
    gl.bindRenderbuffer(gl.RENDERBUFFER, null)
    gl.bindFramebuffer(gl.FRAMEBUFFER, null)
}

type RequiredExtensions = {
    angleInstancedArrays: ANGLE_instanced_arrays
    oesElementIndexUint: OES_element_index_uint
}

export interface Context {
    gl: WebGLRenderingContext
    extensions: RequiredExtensions
    shaderCache: ShaderCache
    programCache: ProgramCache
    destroy: () => void
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
        extensions: { angleInstancedArrays, oesElementIndexUint },
        shaderCache: createShaderCache(),
        programCache: createProgramCache(),
        destroy: () => {
            unbindResources(gl)
            // TODO destroy buffers and textures
        }
    }
}