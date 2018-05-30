/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createProgramCache, ProgramCache } from './program'
import { createShaderCache, ShaderCache } from './shader'

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
    unbindFramebuffer(gl)
}

function unbindFramebuffer(gl: WebGLRenderingContext) {
    gl.bindFramebuffer(gl.FRAMEBUFFER, null)
}

export function createImageData(buffer: Uint8Array, width: number, height: number) {
    const w = width * 4
    const h = height
    const data = new Uint8ClampedArray(width * height * 4)
    for (let i = 0, maxI = h / 2; i < maxI; ++i) {
        for (let j = 0, maxJ = w; j < maxJ; ++j) {
            const index1 = i * w + j;
            const index2 = (h-i-1) * w + j;
            data[index1] = buffer[index2];
            data[index2] = buffer[index1];
        }
    }
    return new ImageData(data, width, height);
}

type Extensions = {
    angleInstancedArrays: ANGLE_instanced_arrays
    standardDerivatives: OES_standard_derivatives
    oesElementIndexUint: OES_element_index_uint | null
    oesVertexArrayObject: OES_vertex_array_object | null
}

/** A WebGL context object, including the rendering context, resource caches and counts */
export interface Context {
    gl: WebGLRenderingContext
    extensions: Extensions
    shaderCache: ShaderCache
    programCache: ProgramCache
    bufferCount: number
    textureCount: number
    vaoCount: number
    readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => void
    unbindFramebuffer: () => void
    destroy: () => void
}

export function createContext(gl: WebGLRenderingContext): Context {
    const angleInstancedArrays = gl.getExtension('ANGLE_instanced_arrays')
    if (angleInstancedArrays === null) {
        throw new Error('Could not get "ANGLE_instanced_arrays" extension')
    }
    const standardDerivatives = gl.getExtension('OES_standard_derivatives')
    if (standardDerivatives === null) {
        throw new Error('Could not get "OES_standard_derivatives" extension')
    }
    const oesElementIndexUint = gl.getExtension('OES_element_index_uint')
    if (oesElementIndexUint === null) {
        console.warn('Could not get "OES_element_index_uint" extension')
    }
    const oesVertexArrayObject = gl.getExtension('OES_vertex_array_object')
    if (oesVertexArrayObject === null) {
        console.log('Could not get "OES_vertex_array_object" extension')
    }

    const shaderCache = createShaderCache()
    const programCache = createProgramCache()

    return {
        gl,
        extensions: { angleInstancedArrays, standardDerivatives, oesElementIndexUint, oesVertexArrayObject },
        shaderCache,
        programCache,
        bufferCount: 0,
        textureCount: 0,
        vaoCount: 0,
        readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => {
            if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE) {
                gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, buffer)
            } else {
                console.error('Reading pixels failed. Framebuffer not complete.')
            }
        },
        unbindFramebuffer: () => unbindFramebuffer(gl),
        destroy: () => {
            unbindResources(gl)
            programCache.dispose()
            shaderCache.dispose()
            // TODO destroy buffers and textures
        }
    }
}