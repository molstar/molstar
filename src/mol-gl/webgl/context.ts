/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createProgramCache, ProgramCache } from './program'
import { createShaderCache, ShaderCache } from './shader'
import { GLRenderingContext, COMPAT_instanced_arrays, COMPAT_standard_derivatives, COMPAT_vertex_array_object, getInstancedArrays, getStandardDerivatives, getVertexArrayObject, isWebGL2, COMPAT_element_index_uint, getElementIndexUint, COMPAT_texture_float, getTextureFloat, COMPAT_texture_float_linear, getTextureFloatLinear, COMPAT_blend_minmax, getBlendMinMax, getFragDepth, COMPAT_frag_depth, COMPAT_color_buffer_float, getColorBufferFloat, COMPAT_draw_buffers, getDrawBuffers, getShaderTextureLod, COMPAT_shader_texture_lod, getDepthTexture, COMPAT_depth_texture } from './compat';
import { createFramebufferCache, FramebufferCache, checkFramebufferStatus } from './framebuffer';
import { Scheduler } from 'mol-task';
import { isDebugMode } from 'mol-util/debug';

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

export function checkError(gl: GLRenderingContext) {
    const error = gl.getError()
    if (error) throw new Error(`WebGL error: '${getErrorDescription(gl, error)}'`)
}

function unbindResources (gl: GLRenderingContext) {
    // bind null to all texture units
    const maxTextureImageUnits = gl.getParameter(gl.MAX_TEXTURE_IMAGE_UNITS)
    for (let i = 0; i < maxTextureImageUnits; ++i) {
        gl.activeTexture(gl.TEXTURE0 + i)
        gl.bindTexture(gl.TEXTURE_2D, null)
        gl.bindTexture(gl.TEXTURE_CUBE_MAP, null)
        if (isWebGL2(gl)) {
            gl.bindTexture(gl.TEXTURE_2D_ARRAY, null)
            gl.bindTexture(gl.TEXTURE_3D, null)
        }
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
    if (isDebugMode) checkFramebufferStatus(gl)
    if (buffer instanceof Uint8Array) {
        gl.readPixels(x, y, width, height, gl.RGBA, gl.UNSIGNED_BYTE, buffer)
    } else if (buffer instanceof Float32Array) {
        gl.readPixels(x, y, width, height, gl.RGBA, gl.FLOAT, buffer)
    } else {
        throw new Error('unsupported readPixels buffer type')
    }
    if (isDebugMode) checkError(gl)
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
    elementIndexUint: COMPAT_element_index_uint
    depthTexture: COMPAT_depth_texture

    vertexArrayObject: COMPAT_vertex_array_object | null
    fragDepth: COMPAT_frag_depth | null
    colorBufferFloat: COMPAT_color_buffer_float | null
    drawBuffers: COMPAT_draw_buffers | null
    shaderTextureLod: COMPAT_shader_texture_lod | null
}

function createExtensions(gl: GLRenderingContext): WebGLExtensions {
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
        throw new Error('Could not find support for "element_index_uint"')
    }
    const depthTexture = getDepthTexture(gl)
    if (depthTexture === null) {
        throw new Error('Could not find support for "depth_texture"')
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
    const shaderTextureLod = getShaderTextureLod(gl)
    if (shaderTextureLod === null) {
        console.log('Could not find support for "shader_texture_lod"')
    }
    
    return {
        instancedArrays,
        standardDerivatives,
        blendMinMax,
        textureFloat,
        textureFloatLinear,
        elementIndexUint,
        depthTexture,

        vertexArrayObject,
        fragDepth,
        colorBufferFloat,
        drawBuffers,
        shaderTextureLod,
    }
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

function createStats(): WebGLStats {
    return {
        bufferCount: 0,
        framebufferCount: 0,
        renderbufferCount: 0,
        textureCount: 0,
        vaoCount: 0,

        drawCount: 0,
        instanceCount: 0,
        instancedDrawCount: 0,
    }
}

export type WebGLState = {
    currentProgramId: number
    currentMaterialId: number
    currentRenderItemId: number

    /**
     * specifies which WebGL capability to enable
     * - `gl.BLEND`: blending of the computed fragment color values
     * - `gl.CULL_FACE`: culling of polygons
     * - `gl.DEPTH_TEST`: depth comparisons and updates to the depth buffer
     * - `gl.DITHER`: dithering of color components before they get written to the color buffer
     * - `gl.POLYGON_OFFSET_FILL`: adding an offset to depth values of polygon's fragments
     * - `gl.SAMPLE_ALPHA_TO_COVERAGE`: computation of a temporary coverage value determined by the alpha value
     * - `gl.SAMPLE_COVERAGE`: ANDing the fragment's coverage with the temporary coverage value
     * - `gl.SCISSOR_TEST`: scissor test that discards fragments that are outside of the scissor rectangle
     * - `gl.STENCIL_TEST`: stencil testing and updates to the stencil buffer
     */
    enable: (cap: number) => void
    /**
     * specifies which WebGL capability to disable
     * - `gl.BLEND`: blending of the computed fragment color values
     * - `gl.CULL_FACE`: culling of polygons
     * - `gl.DEPTH_TEST`: depth comparisons and updates to the depth buffer
     * - `gl.DITHER`: dithering of color components before they get written to the color buffer
     * - `gl.POLYGON_OFFSET_FILL`: adding an offset to depth values of polygon's fragments
     * - `gl.SAMPLE_ALPHA_TO_COVERAGE`: computation of a temporary coverage value determined by the alpha value
     * - `gl.SAMPLE_COVERAGE`: ANDing the fragment's coverage with the temporary coverage value
     * - `gl.SCISSOR_TEST`: scissor test that discards fragments that are outside of the scissor rectangle
     * - `gl.STENCIL_TEST`: stencil testing and updates to the stencil buffer
     */
    disable: (cap: number) => void

    /** specifies whether polygons are front- or back-facing by setting a winding orientation */
    frontFace: (mode: number) => void
    /** specifies whether or not front- and/or back-facing polygons can be culled */
    cullFace: (mode: number) => void
    /** sets whether writing into the depth buffer is enabled or disabled */
    depthMask: (flag: boolean) => void
    /** sets which color components to enable or to disable */
    colorMask: (red: boolean, green: boolean, blue: boolean, alpha: boolean) => void
    /** specifies the color values used when clearing color buffers, used when calling `gl.clear`, clamped to [0, 1] */
    clearColor: (red: number, green: number, blue: number, alpha: number) => void

    /** defines which function is used for blending pixel arithmetic */
    blendFunc: (src: number, dst: number) => void
    /** defines which function is used for blending pixel arithmetic for RGB and alpha components separately */
    blendFuncSeparate: (srcRGB: number, dstRGB: number, srcAlpha: number, dstAlpha: number) => void

    /** set both the RGB blend equation and alpha blend equation to a single equation, determines how a new pixel is combined with an existing */
    blendEquation: (mode: number) => void
    /** set the RGB blend equation and alpha blend equation separately, determines how a new pixel is combined with an existing */
    blendEquationSeparate: (modeRGB: number, modeAlpha: number) => void
}

function createState(gl: GLRenderingContext): WebGLState {
    const enabledCapabilities: { [k: number]: boolean } = {}

    let currentFrontFace = gl.getParameter(gl.FRONT_FACE)
    let currentCullFace = gl.getParameter(gl.CULL_FACE_MODE)
    let currentDepthMask = gl.getParameter(gl.DEPTH_WRITEMASK)
    let currentColorMask = gl.getParameter(gl.COLOR_WRITEMASK)
    let currentClearColor = gl.getParameter(gl.COLOR_CLEAR_VALUE)

    let currentBlendSrcRGB = gl.getParameter(gl.BLEND_SRC_RGB)
    let currentBlendDstRGB = gl.getParameter(gl.BLEND_DST_RGB)
    let currentBlendSrcAlpha = gl.getParameter(gl.BLEND_SRC_ALPHA)
    let currentBlendDstAlpha = gl.getParameter(gl.BLEND_DST_ALPHA)

    let currentBlendEqRGB = gl.getParameter(gl.BLEND_EQUATION_RGB)
    let currentBlendEqAlpha = gl.getParameter(gl.BLEND_EQUATION_ALPHA)

    return {
        currentProgramId: -1,
        currentMaterialId: -1,
        currentRenderItemId: -1,

        enable: (cap: number) => {
            if (enabledCapabilities[cap] !== true ) {
                gl.enable(cap)
                enabledCapabilities[cap] = true
            }
        },
        disable: (cap: number) => {
            if (enabledCapabilities[cap] !== false) {
                gl.disable(cap)
                enabledCapabilities[cap] = false
            }
        },

        frontFace: (mode: number) => {
            if (mode !== currentFrontFace) {
                gl.frontFace(mode)
                currentFrontFace = mode
            }
        },
        cullFace: (mode: number) => {
            if (mode !== currentCullFace) {
                gl.cullFace(mode)
                currentCullFace = mode
            }
        },
        depthMask: (flag: boolean) => {
            if (flag !== currentDepthMask) {
                gl.depthMask(flag)
                currentDepthMask = flag
            }
        },
        colorMask: (red: boolean, green: boolean, blue: boolean, alpha: boolean) => {
            if (red !== currentColorMask[0] || green !== currentColorMask[1] || blue !== currentColorMask[2] || alpha !== currentColorMask[3])
            gl.colorMask(red, green, blue, alpha)
            currentColorMask[0] = red
            currentColorMask[1] = green
            currentColorMask[2] = blue
            currentColorMask[3] = alpha
        },
        clearColor: (red: number, green: number, blue: number, alpha: number) => {
            if (red !== currentClearColor[0] || green !== currentClearColor[1] || blue !== currentClearColor[2] || alpha !== currentClearColor[3])
            gl.clearColor(red, green, blue, alpha)
            currentClearColor[0] = red
            currentClearColor[1] = green
            currentClearColor[2] = blue
            currentClearColor[3] = alpha
        },

        blendFunc: (src: number, dst: number) => {
            if (src !== currentBlendSrcRGB || dst !== currentBlendDstRGB || src !== currentBlendSrcAlpha || dst !== currentBlendDstAlpha) {
                gl.blendFunc(src, dst)
                currentBlendSrcRGB = src
                currentBlendDstRGB = dst
                currentBlendSrcAlpha = src
                currentBlendDstAlpha = dst
            }
        },
        blendFuncSeparate: (srcRGB: number, dstRGB: number, srcAlpha: number, dstAlpha: number) => {
            if (srcRGB !== currentBlendSrcRGB || dstRGB !== currentBlendDstRGB || srcAlpha !== currentBlendSrcAlpha || dstAlpha !== currentBlendDstAlpha) {
                gl.blendFuncSeparate(srcRGB, dstRGB, srcAlpha, dstAlpha)
                currentBlendSrcRGB = srcRGB
                currentBlendDstRGB = dstRGB
                currentBlendSrcAlpha = srcAlpha
                currentBlendDstAlpha = dstAlpha
            }
        },

        blendEquation: (mode: number) => {
            if (mode !== currentBlendEqRGB || mode !== currentBlendEqAlpha) {
                gl.blendEquation(mode)
                currentBlendEqRGB = mode
                currentBlendEqAlpha = mode
            }
        },
        blendEquationSeparate: (modeRGB: number, modeAlpha: number) => {
            if (modeRGB !== currentBlendEqRGB || modeAlpha !== currentBlendEqAlpha) {
                gl.blendEquationSeparate(modeRGB, modeAlpha)
                currentBlendEqRGB = modeRGB
                currentBlendEqAlpha = modeAlpha
            }
        }
    }
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
    const extensions = createExtensions(gl)
    const state = createState(gl)
    const stats = createStats()

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