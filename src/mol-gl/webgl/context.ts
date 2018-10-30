/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createProgramCache, ProgramCache } from './program'
import { createShaderCache, ShaderCache } from './shader'
import { GLRenderingContext, COMPAT_instanced_arrays, COMPAT_standard_derivatives, COMPAT_vertex_array_object, getInstancedArrays, getStandardDerivatives, getVertexArrayObject, isWebGL2, COMPAT_element_index_uint, getElementIndexUint, COMPAT_texture_float, getTextureFloat, COMPAT_texture_float_linear, getTextureFloatLinear, COMPAT_blend_minmax, getBlendMinMax, getFragDepth, COMPAT_frag_depth } from './compat';

export function getGLContext(canvas: HTMLCanvasElement, contextAttributes?: WebGLContextAttributes): GLRenderingContext | null {
    function getContext(contextId: 'webgl' | 'experimental-webgl' | 'webgl2') {
        try {
           return canvas.getContext(contextId, contextAttributes) as GLRenderingContext | null
        } catch (e) {
            return null
        }
    }
    return getContext('webgl2') ||  getContext('webgl') || getContext('experimental-webgl')
}

function getPixelRatio() {
    return (typeof window !== 'undefined') ? window.devicePixelRatio : 1
}

function unbindResources (gl: GLRenderingContext) {
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

function unbindFramebuffer(gl: GLRenderingContext) {
    gl.bindFramebuffer(gl.FRAMEBUFFER, null)
}

const tmpPixel = new Uint8Array(1 * 4);

function fence(gl: WebGL2RenderingContext) {
    return new Promise(resolve => {
        const sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0)
        if (!sync) {
            console.warn('could not create a WebGL2 sync object')
            gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, tmpPixel)
            resolve()
        } else {
            gl.flush(); // Ensure the fence is submitted.
            const check = () => {
                const status = gl.getSyncParameter(sync, gl.SYNC_STATUS)
                if (status == gl.SIGNALED) {
                    gl.deleteSync(sync);
                    resolve();
                } else {
                    setTimeout(check, 0)
                }
            }
            setTimeout(check, 0)
        }
    })
}

async function waitForGpuCommandsComplete(gl: GLRenderingContext) {
    if (isWebGL2(gl)) {
        await fence(gl)
    } else {
        console.info('webgl sync object not supported in webgl 1')
        gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, tmpPixel)
    }
}

export function createImageData(buffer: ArrayLike<number>, width: number, height: number) {
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

//

type Extensions = {
    instancedArrays: COMPAT_instanced_arrays
    standardDerivatives: COMPAT_standard_derivatives
    blendMinMax: COMPAT_blend_minmax
    textureFloat: COMPAT_texture_float
    textureFloatLinear: COMPAT_texture_float_linear
    elementIndexUint: COMPAT_element_index_uint | null
    vertexArrayObject: COMPAT_vertex_array_object | null
    fragDepth: COMPAT_frag_depth | null
}

/** A WebGL context object, including the rendering context, resource caches and counts */
export interface Context {
    readonly gl: GLRenderingContext
    readonly isWebGL2: boolean
    readonly extensions: Extensions
    readonly pixelRatio: number

    readonly shaderCache: ShaderCache
    readonly programCache: ProgramCache

    bufferCount: number
    framebufferCount: number
    renderbufferCount: number
    textureCount: number
    vaoCount: number

    drawCount: number
    instanceCount: number
    instancedDrawCount: number

    readonly maxTextureSize: number
    readonly maxDrawBuffers: number

    unbindFramebuffer: () => void
    readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => void
    readPixelsAsync: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => Promise<void>
    waitForGpuCommandsComplete: () => Promise<void>
    destroy: () => void
}

export function createContext(gl: GLRenderingContext): Context {
    const instancedArrays = getInstancedArrays(gl)
    if (instancedArrays === null) {
        throw new Error('Could not find support for "instanced_arrays"')
    }
    const standardDerivatives = getStandardDerivatives(gl)
    if (standardDerivatives === null) {
        throw new Error('Could not find support for "standard_derivatives"')
    }
    const blendMinMax = getBlendMinMax(gl)
    if (blendMinMax === null) {
        throw new Error('Could not find support for "blend_minmax"')
    }
    const textureFloat = getTextureFloat(gl)
    if (textureFloat === null) {
        throw new Error('Could not find support for "texture_float"')
    }
    const textureFloatLinear = getTextureFloatLinear(gl)
    if (textureFloatLinear === null) {
        throw new Error('Could not find support for "texture_float_linear"')
    }
    const elementIndexUint = getElementIndexUint(gl)
    if (elementIndexUint === null) {
        console.warn('Could not find support for "element_index_uint"')
    }
    const vertexArrayObject = getVertexArrayObject(gl)
    if (vertexArrayObject === null) {
        console.log('Could not find support for "vertex_array_object"')
    }
    const fragDepth = getFragDepth(gl)
    if (fragDepth === null) {
        console.log('Could not find support for "frag_depth"')
    }

    const shaderCache = createShaderCache()
    const programCache = createProgramCache()

    const parameters = {
        maxTextureSize: gl.getParameter(gl.MAX_TEXTURE_SIZE),
        maxDrawBuffers: isWebGL2(gl) ? gl.getParameter(gl.MAX_DRAW_BUFFERS) : 0,
        maxVertexTextureImageUnits: gl.getParameter(gl.MAX_VERTEX_TEXTURE_IMAGE_UNITS),
    }

    if (parameters.maxVertexTextureImageUnits < 4) {
        throw new Error('Need "MAX_VERTEX_TEXTURE_IMAGE_UNITS" >= 4')
    }

    let readPixelsAsync: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => Promise<void>
    if (isWebGL2(gl)) {
        const pbo = gl.createBuffer()
        readPixelsAsync = async (x: number, y: number, width: number, height: number, buffer: Uint8Array) => {
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, pbo)
            gl.bufferData(gl.PIXEL_PACK_BUFFER, width * height * 4, gl.STATIC_COPY)
            gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, 0)
            await fence(gl)
            gl.getBufferSubData(gl.PIXEL_PACK_BUFFER, 0, buffer);
        }
    } else {
        readPixelsAsync = async (x: number, y: number, width: number, height: number, buffer: Uint8Array) => {
            gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, buffer)
        }
    }

    return {
        gl,
        isWebGL2: isWebGL2(gl),
        extensions: {
            instancedArrays,
            standardDerivatives,
            blendMinMax,
            textureFloat,
            textureFloatLinear,
            elementIndexUint,
            vertexArrayObject,
            fragDepth
        },
        pixelRatio: getPixelRatio(),

        shaderCache,
        programCache,

        bufferCount: 0,
        framebufferCount: 0,
        renderbufferCount: 0,
        textureCount: 0,
        vaoCount: 0,

        drawCount: 0,
        instanceCount: 0,
        instancedDrawCount: 0,

        get maxTextureSize () { return parameters.maxTextureSize },
        get maxDrawBuffers () { return parameters.maxDrawBuffers },

        unbindFramebuffer: () => unbindFramebuffer(gl),
        readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => {
            gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, buffer)
            // TODO check is very expensive
            // if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE) {
            //     gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, buffer)
            // } else {
            //     console.error('Reading pixels failed. Framebuffer not complete.')
            // }
        },
        readPixelsAsync,
        waitForGpuCommandsComplete: () => waitForGpuCommandsComplete(gl),

        destroy: () => {
            unbindResources(gl)
            programCache.dispose()
            shaderCache.dispose()
            // TODO destroy buffers and textures
        }
    }
}