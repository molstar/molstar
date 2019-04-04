/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createProgramCache, ProgramCache } from './program'
import { createShaderCache, ShaderCache } from './shader'
import { GLRenderingContext, COMPAT_instanced_arrays, COMPAT_standard_derivatives, COMPAT_vertex_array_object, getInstancedArrays, getStandardDerivatives, getVertexArrayObject, isWebGL2, COMPAT_element_index_uint, getElementIndexUint, COMPAT_texture_float, getTextureFloat, COMPAT_texture_float_linear, getTextureFloatLinear, COMPAT_blend_minmax, getBlendMinMax, getFragDepth, COMPAT_frag_depth, COMPAT_color_buffer_float, getColorBufferFloat, COMPAT_draw_buffers, getDrawBuffers } from './compat';
import { createFramebufferCache, FramebufferCache } from './framebuffer';
import { Scheduler } from 'mol-task';
import { isProductionMode } from 'mol-util/debug';

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

function getErrorDescription(gl: GLRenderingContext, error: number) {
    switch (error) {
        case gl.NO_ERROR: return 'no error'
        case gl.INVALID_ENUM: return 'invalid enum'
        case gl.INVALID_VALUE: return 'invalid value'
        case gl.INVALID_OPERATION: return 'invalid operation'
        case gl.INVALID_FRAMEBUFFER_OPERATION: return 'invalid framebuffer operation'
        case gl.OUT_OF_MEMORY: return 'out of memory'
        case gl.CONTEXT_LOST_WEBGL: return 'context lost'
    }
    return 'unknown error'
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

function checkSync(gl: WebGL2RenderingContext, sync: WebGLSync, resolve: () => void) {
    if (gl.getSyncParameter(sync, gl.SYNC_STATUS) === gl.SIGNALED) {
        gl.deleteSync(sync)
        resolve()
    } else {
        Scheduler.setImmediate(checkSync, gl, sync, resolve)
    }
}

function fence(gl: WebGL2RenderingContext, resolve: () => void) {
    const sync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0)
    if (!sync) {
        console.warn('Could not create a WebGLSync object')
        gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, tmpPixel)
        resolve()
    } else {
        Scheduler.setImmediate(checkSync, gl, sync, resolve)
    }
}

let SentWebglSyncObjectNotSupportedInWebglMessage = false
function waitForGpuCommandsComplete(gl: GLRenderingContext): Promise<void> {
    return new Promise(resolve => {
        if (isWebGL2(gl)) {
            // TODO seems quite slow
            fence(gl, resolve)
        } else {
            if (!SentWebglSyncObjectNotSupportedInWebglMessage) {
                console.info('Sync object not supported in WebGL')
                SentWebglSyncObjectNotSupportedInWebglMessage = true
            }
            waitForGpuCommandsCompleteSync(gl)
            resolve()
        }
    })
}

function waitForGpuCommandsCompleteSync(gl: GLRenderingContext): void {
    gl.bindFramebuffer(gl.FRAMEBUFFER, null)
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, tmpPixel)
}

function readPixels(gl: GLRenderingContext, x: number, y: number, width: number, height: number, buffer: Uint8Array | Float32Array) {
    if (!isProductionMode && gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) {
        console.error('Reading pixels failed. Framebuffer not complete.')
        return
    }
    if (buffer instanceof Uint8Array) {
        gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, buffer)
    } else {
        gl.readPixels(x, y, width, height, gl.RGBA, gl.FLOAT, buffer)
    }
    if (!isProductionMode) {
        const error = gl.getError()
        if (error) console.log(`Error reading pixels: '${getErrorDescription(gl, error)}'`)
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

export type WebGLExtensions = {
    instancedArrays: COMPAT_instanced_arrays
    standardDerivatives: COMPAT_standard_derivatives
    blendMinMax: COMPAT_blend_minmax
    textureFloat: COMPAT_texture_float
    textureFloatLinear: COMPAT_texture_float_linear
    elementIndexUint: COMPAT_element_index_uint | null
    vertexArrayObject: COMPAT_vertex_array_object | null
    fragDepth: COMPAT_frag_depth | null
    colorBufferFloat: COMPAT_color_buffer_float | null
    drawBuffers: COMPAT_draw_buffers | null
}

export type WebGLStats = {
    bufferCount: number
    framebufferCount: number
    renderbufferCount: number
    textureCount: number
    vaoCount: number

    drawCount: number
    instanceCount: number
    instancedDrawCount: number
}

export type WebGLState = {
    currentProgramId: number
    currentMaterialId: number
}

/** A WebGL context object, including the rendering context, resource caches and counts */
export interface WebGLContext {
    readonly gl: GLRenderingContext
    readonly isWebGL2: boolean
    readonly pixelRatio: number

    readonly extensions: WebGLExtensions
    readonly state: WebGLState
    readonly stats: WebGLStats

    readonly shaderCache: ShaderCache
    readonly programCache: ProgramCache
    readonly framebufferCache: FramebufferCache

    readonly maxTextureSize: number
    readonly maxDrawBuffers: number

    unbindFramebuffer: () => void
    readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array | Float32Array) => void
    readPixelsAsync: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => Promise<void>
    waitForGpuCommandsComplete: () => Promise<void>
    waitForGpuCommandsCompleteSync: () => void
    destroy: () => void
}

export function createContext(gl: GLRenderingContext): WebGLContext {
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
    const colorBufferFloat = getColorBufferFloat(gl)
    if (colorBufferFloat === null) {
        console.log('Could not find support for "color_buffer_float"')
    }
    const drawBuffers = getDrawBuffers(gl)
    if (drawBuffers === null) {
        console.log('Could not find support for "draw_buffers"')
    }

    const state: WebGLState = {
        currentProgramId: -1,
        currentMaterialId: -1,
    }

    const stats: WebGLStats = {
        bufferCount: 0,
        framebufferCount: 0,
        renderbufferCount: 0,
        textureCount: 0,
        vaoCount: 0,

        drawCount: 0,
        instanceCount: 0,
        instancedDrawCount: 0,
    }

    const extensions: WebGLExtensions = {
        instancedArrays,
        standardDerivatives,
        blendMinMax,
        textureFloat,
        textureFloatLinear,
        elementIndexUint,
        vertexArrayObject,
        fragDepth,
        colorBufferFloat,
        drawBuffers
    }

    const shaderCache: ShaderCache = createShaderCache(gl)
    const programCache: ProgramCache = createProgramCache(gl, state, extensions, shaderCache)
    const framebufferCache: FramebufferCache = createFramebufferCache(gl, stats)

    const parameters = {
        maxTextureSize: gl.getParameter(gl.MAX_TEXTURE_SIZE) as number,
        maxDrawBuffers: isWebGL2(gl) ? gl.getParameter(gl.MAX_DRAW_BUFFERS) as number : 0,
        maxVertexTextureImageUnits: gl.getParameter(gl.MAX_VERTEX_TEXTURE_IMAGE_UNITS) as number,
    }

    if (parameters.maxVertexTextureImageUnits < 8) {
        throw new Error('Need "MAX_VERTEX_TEXTURE_IMAGE_UNITS" >= 8')
    }

    let readPixelsAsync: (x: number, y: number, width: number, height: number, buffer: Uint8Array) => Promise<void>
    if (isWebGL2(gl)) {
        const pbo = gl.createBuffer()
        let _buffer: Uint8Array | undefined = void 0
        let _resolve: (() => void) | undefined = void 0
        let _reading = false

        const bindPBO = () => {
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, pbo)
            gl.getBufferSubData(gl.PIXEL_PACK_BUFFER, 0, _buffer!)
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null)
            _reading = false
            _resolve!()
            _resolve = void 0
            _buffer = void 0
        }
        readPixelsAsync = (x: number, y: number, width: number, height: number, buffer: Uint8Array): Promise<void> => new Promise<void>((resolve, reject) => {
            if (_reading) {
                reject('Can not call multiple readPixelsAsync at the same time')
                return
            }
            _reading = true;
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, pbo)
            gl.bufferData(gl.PIXEL_PACK_BUFFER, width * height * 4, gl.STREAM_READ)
            gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, 0)
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null)
            // need to unbind/bind PBO before/after async awaiting the fence
            _resolve = resolve
            _buffer = buffer
            fence(gl, bindPBO)
        })
    } else {
        readPixelsAsync = async (x: number, y: number, width: number, height: number, buffer: Uint8Array) => {
            readPixels(gl, x, y, width, height, buffer)
        }
    }

    return {
        gl,
        isWebGL2: isWebGL2(gl),
        get pixelRatio () {
            // this can change during the lifetime of a rendering context, so need to re-obtain on access
            return getPixelRatio()
        },

        extensions,
        state,
        stats,

        shaderCache,
        programCache,
        framebufferCache,

        get maxTextureSize () { return parameters.maxTextureSize },
        get maxDrawBuffers () { return parameters.maxDrawBuffers },

        unbindFramebuffer: () => unbindFramebuffer(gl),
        readPixels: (x: number, y: number, width: number, height: number, buffer: Uint8Array | Float32Array) => {
            readPixels(gl, x, y, width, height, buffer)
        },
        readPixelsAsync,
        waitForGpuCommandsComplete: () => waitForGpuCommandsComplete(gl),
        waitForGpuCommandsCompleteSync: () => waitForGpuCommandsCompleteSync(gl),

        destroy: () => {
            unbindResources(gl)
            programCache.dispose()
            shaderCache.dispose()
            framebufferCache.dispose()
            // TODO destroy buffers and textures
        }
    }
}