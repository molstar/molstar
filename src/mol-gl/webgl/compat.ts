/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { isDebugMode } from '../../mol-util/debug';
import { getErrorDescription } from './context';
import { getProgram } from './program';
import { getShader } from './shader';

export type GLRenderingContext = WebGLRenderingContext | WebGL2RenderingContext

export function isWebGL(gl: any): gl is WebGLRenderingContext {
    return typeof WebGLRenderingContext !== 'undefined' && gl instanceof WebGLRenderingContext;
}

export function isWebGL2(gl: any): gl is WebGL2RenderingContext {
    return typeof WebGL2RenderingContext !== 'undefined' && gl instanceof WebGL2RenderingContext;
}

/**
 * See https://registry.khronos.org/webgl/extensions/ANGLE_instanced_arrays/
 */
export interface COMPAT_instanced_arrays {
    /**
     * Renders primitives from array data like the `drawArrays` method. In addition, it can execute multiple instances of the range of elements.
     * @param mode the type primitive to render.
     * @param first the starting index in the array of vector points.
     * @param count the number of indices to be rendered.
     * @param primcount the number of instances of the range of elements to execute.
     */
    drawArraysInstanced(mode: number, first: number, count: number, primcount: number): void;
    /**
     * Renders primitives from array data like the `drawElements` method. In addition, it can execute multiple instances of a set of elements.
     * @param mode the type primitive to render.
     * @param count the number of elements to be rendered.
     * @param type the type of the values in the element array buffer.
     * @param offset an offset in the element array buffer. Must be a valid multiple of the size of the given `type`.
     * @param primcount the number of instances of the set of elements to execute.
     */
    drawElementsInstanced(mode: number, count: number, type: number, offset: number, primcount: number): void;
    /**
     * Modifies the rate at which generic vertex attributes advance when rendering multiple instances of primitives with `drawArraysInstanced` and `drawElementsInstanced`
     * @param index the index of the generic vertex attributes.
     * @param divisor the number of instances that will pass between updates of the generic attribute.
     */
    vertexAttribDivisor(index: number, divisor: number): void;
    readonly VERTEX_ATTRIB_ARRAY_DIVISOR: number;
}

export function getInstancedArrays(gl: GLRenderingContext): COMPAT_instanced_arrays | null {
    if (isWebGL2(gl)) {
        return {
            drawArraysInstanced: gl.drawArraysInstanced.bind(gl),
            drawElementsInstanced: gl.drawElementsInstanced.bind(gl),
            vertexAttribDivisor: gl.vertexAttribDivisor.bind(gl),
            VERTEX_ATTRIB_ARRAY_DIVISOR: gl.VERTEX_ATTRIB_ARRAY_DIVISOR
        };
    } else {
        const ext = gl.getExtension('ANGLE_instanced_arrays');
        if (ext === null) return null;
        return {
            drawArraysInstanced: ext.drawArraysInstancedANGLE.bind(ext),
            drawElementsInstanced: ext.drawElementsInstancedANGLE.bind(ext),
            vertexAttribDivisor: ext.vertexAttribDivisorANGLE.bind(ext),
            VERTEX_ATTRIB_ARRAY_DIVISOR: ext.VERTEX_ATTRIB_ARRAY_DIVISOR_ANGLE
        };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_standard_derivatives/
 */
export interface COMPAT_standard_derivatives {
    readonly FRAGMENT_SHADER_DERIVATIVE_HINT: number;
}

export function getStandardDerivatives(gl: GLRenderingContext): COMPAT_standard_derivatives | null {
    if (isWebGL2(gl)) {
        return { FRAGMENT_SHADER_DERIVATIVE_HINT: gl.FRAGMENT_SHADER_DERIVATIVE_HINT };
    } else {
        const ext = gl.getExtension('OES_standard_derivatives');
        if (ext === null) return null;
        return { FRAGMENT_SHADER_DERIVATIVE_HINT: ext.FRAGMENT_SHADER_DERIVATIVE_HINT_OES };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_element_index_uint/
 */
export interface COMPAT_element_index_uint {
}

export function getElementIndexUint(gl: GLRenderingContext): COMPAT_element_index_uint | null {
    return isWebGL2(gl) ? {} : gl.getExtension('OES_element_index_uint');
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_vertex_array_object/
 */
export interface COMPAT_vertex_array_object {
    readonly VERTEX_ARRAY_BINDING: number;
    bindVertexArray(arrayObject: WebGLVertexArrayObject | null): void;
    createVertexArray(): WebGLVertexArrayObject | null;
    deleteVertexArray(arrayObject: WebGLVertexArrayObject): void;
    isVertexArray(value: any): value is WebGLVertexArrayObject;
}

export function getVertexArrayObject(gl: GLRenderingContext): COMPAT_vertex_array_object | null {
    if (isWebGL2(gl)) {
        return {
            VERTEX_ARRAY_BINDING: gl.VERTEX_ARRAY_BINDING,
            bindVertexArray: gl.bindVertexArray.bind(gl),
            createVertexArray: gl.createVertexArray.bind(gl),
            deleteVertexArray: gl.deleteVertexArray.bind(gl),
            isVertexArray: gl.isVertexArray.bind(gl)
        };
    } else {
        const ext = gl.getExtension('OES_vertex_array_object');
        if (ext === null) return null;
        return {
            VERTEX_ARRAY_BINDING: ext.VERTEX_ARRAY_BINDING_OES,
            bindVertexArray: ext.bindVertexArrayOES.bind(ext),
            createVertexArray: ext.createVertexArrayOES.bind(ext),
            deleteVertexArray: ext.deleteVertexArrayOES.bind(ext),
            isVertexArray: ext.isVertexArrayOES.bind(ext)
        };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_texture_float/
 */
export interface COMPAT_texture_float {
}

export function getTextureFloat(gl: GLRenderingContext): COMPAT_texture_float | null {
    return isWebGL2(gl) ? {} : gl.getExtension('OES_texture_float');
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_texture_float_linear/
 */
export interface COMPAT_texture_float_linear {
}

export function getTextureFloatLinear(gl: GLRenderingContext): COMPAT_texture_float_linear | null {
    return gl.getExtension('OES_texture_float_linear');
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_texture_half_float/
 */
export interface COMPAT_texture_half_float {
    readonly HALF_FLOAT: number
}

export function getTextureHalfFloat(gl: GLRenderingContext): COMPAT_texture_half_float | null {
    if (isWebGL2(gl)) {
        return { HALF_FLOAT: gl.HALF_FLOAT };
    } else {
        const ext = gl.getExtension('OES_texture_half_float');
        if (ext === null) return null;
        return { HALF_FLOAT: ext.HALF_FLOAT_OES };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_texture_half_float_linear/
 */
export interface COMPAT_texture_half_float_linear {
}

export function getTextureHalfFloatLinear(gl: GLRenderingContext): COMPAT_texture_half_float_linear | null {
    return gl.getExtension('OES_texture_half_float_linear');
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_blend_minmax/
 */
export interface COMPAT_blend_minmax {
    readonly MIN: number
    readonly MAX: number
}

export function getBlendMinMax(gl: GLRenderingContext): COMPAT_blend_minmax | null {
    if (isWebGL2(gl)) {
        return { MIN: gl.MIN, MAX: gl.MAX };
    } else {
        const ext = gl.getExtension('EXT_blend_minmax');
        if (ext === null) return null;
        return { MIN: ext.MIN_EXT, MAX: ext.MAX_EXT };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_frag_depth/
 */
export interface COMPAT_frag_depth {
}

export function getFragDepth(gl: GLRenderingContext): COMPAT_frag_depth | null {
    return isWebGL2(gl) ? {} : gl.getExtension('EXT_frag_depth');
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_color_buffer_float/
 */
export interface COMPAT_color_buffer_float {
    readonly RGBA32F: number;
}

export function getColorBufferFloat(gl: GLRenderingContext): COMPAT_color_buffer_float | null {
    if (isWebGL2(gl)) {
        if (gl.getExtension('EXT_color_buffer_float') === null) return null;
        gl.getExtension('EXT_float_blend');
        return { RGBA32F: gl.RGBA32F };
    } else {
        const ext = gl.getExtension('WEBGL_color_buffer_float');
        if (ext === null) {
            // test as support may not be advertised by browsers
            gl.getExtension('OES_texture_float');
            return testColorBuffer(gl, gl.FLOAT) ? { RGBA32F: 0x8814 } : null;
        }
        gl.getExtension('EXT_float_blend');
        return { RGBA32F: ext.RGBA32F_EXT };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_color_buffer_half_float/
 */
export interface COMPAT_color_buffer_half_float {
    readonly RGBA16F: number;
}

export function getColorBufferHalfFloat(gl: GLRenderingContext): COMPAT_color_buffer_half_float | null {
    if (isWebGL2(gl)) {
        if (gl.getExtension('EXT_color_buffer_half_float') === null) return null;
        gl.getExtension('EXT_float_blend');
        return { RGBA16F: gl.RGBA16F };
    } else {
        const ext = gl.getExtension('EXT_color_buffer_half_float');
        if (ext === null) {
            // test as support may not be advertised by browsers
            gl.getExtension('OES_texture_half_float');
            return testColorBuffer(gl, 0x8D61) ? { RGBA16F: 0x881A } : null;
        }
        gl.getExtension('EXT_float_blend');
        return { RGBA16F: ext.RGBA16F_EXT };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_draw_buffers/
 */
export interface COMPAT_draw_buffers {
    drawBuffers(buffers: number[]): void;
    readonly COLOR_ATTACHMENT0: number;
    readonly COLOR_ATTACHMENT1: number;
    readonly COLOR_ATTACHMENT2: number;
    readonly COLOR_ATTACHMENT3: number;
    readonly COLOR_ATTACHMENT4: number;
    readonly COLOR_ATTACHMENT5: number;
    readonly COLOR_ATTACHMENT6: number;
    readonly COLOR_ATTACHMENT7: number;
    readonly DRAW_BUFFER0: number;
    readonly DRAW_BUFFER1: number;
    readonly DRAW_BUFFER2: number;
    readonly DRAW_BUFFER3: number;
    readonly DRAW_BUFFER4: number;
    readonly DRAW_BUFFER5: number;
    readonly DRAW_BUFFER6: number;
    readonly DRAW_BUFFER7: number;
    readonly MAX_COLOR_ATTACHMENTS: number;
    readonly MAX_DRAW_BUFFERS: number;
}

export function getDrawBuffers(gl: GLRenderingContext): COMPAT_draw_buffers | null {
    if (isWebGL2(gl)) {
        return {
            drawBuffers: gl.drawBuffers.bind(gl),
            COLOR_ATTACHMENT0: gl.COLOR_ATTACHMENT0,
            COLOR_ATTACHMENT1: gl.COLOR_ATTACHMENT1,
            COLOR_ATTACHMENT2: gl.COLOR_ATTACHMENT2,
            COLOR_ATTACHMENT3: gl.COLOR_ATTACHMENT3,
            COLOR_ATTACHMENT4: gl.COLOR_ATTACHMENT4,
            COLOR_ATTACHMENT5: gl.COLOR_ATTACHMENT5,
            COLOR_ATTACHMENT6: gl.COLOR_ATTACHMENT6,
            COLOR_ATTACHMENT7: gl.COLOR_ATTACHMENT7,
            DRAW_BUFFER0: gl.DRAW_BUFFER0,
            DRAW_BUFFER1: gl.DRAW_BUFFER1,
            DRAW_BUFFER2: gl.DRAW_BUFFER2,
            DRAW_BUFFER3: gl.DRAW_BUFFER3,
            DRAW_BUFFER4: gl.DRAW_BUFFER4,
            DRAW_BUFFER5: gl.DRAW_BUFFER5,
            DRAW_BUFFER6: gl.DRAW_BUFFER6,
            DRAW_BUFFER7: gl.DRAW_BUFFER7,
            MAX_COLOR_ATTACHMENTS: gl.MAX_COLOR_ATTACHMENTS,
            MAX_DRAW_BUFFERS: gl.MAX_DRAW_BUFFERS,
        };
    } else {
        const ext = gl.getExtension('WEBGL_draw_buffers');
        if (ext === null) return null;
        return {
            drawBuffers: ext.drawBuffersWEBGL.bind(ext),
            COLOR_ATTACHMENT0: ext.COLOR_ATTACHMENT0_WEBGL,
            COLOR_ATTACHMENT1: ext.COLOR_ATTACHMENT1_WEBGL,
            COLOR_ATTACHMENT2: ext.COLOR_ATTACHMENT2_WEBGL,
            COLOR_ATTACHMENT3: ext.COLOR_ATTACHMENT3_WEBGL,
            COLOR_ATTACHMENT4: ext.COLOR_ATTACHMENT4_WEBGL,
            COLOR_ATTACHMENT5: ext.COLOR_ATTACHMENT5_WEBGL,
            COLOR_ATTACHMENT6: ext.COLOR_ATTACHMENT6_WEBGL,
            COLOR_ATTACHMENT7: ext.COLOR_ATTACHMENT7_WEBGL,
            DRAW_BUFFER0: ext.DRAW_BUFFER0_WEBGL,
            DRAW_BUFFER1: ext.DRAW_BUFFER1_WEBGL,
            DRAW_BUFFER2: ext.DRAW_BUFFER2_WEBGL,
            DRAW_BUFFER3: ext.DRAW_BUFFER3_WEBGL,
            DRAW_BUFFER4: ext.DRAW_BUFFER4_WEBGL,
            DRAW_BUFFER5: ext.DRAW_BUFFER5_WEBGL,
            DRAW_BUFFER6: ext.DRAW_BUFFER6_WEBGL,
            DRAW_BUFFER7: ext.DRAW_BUFFER7_WEBGL,
            MAX_COLOR_ATTACHMENTS: ext.MAX_COLOR_ATTACHMENTS_WEBGL,
            MAX_DRAW_BUFFERS: ext.MAX_DRAW_BUFFERS_WEBGL,
        };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_draw_buffers_indexed/
 */
export interface COMPAT_draw_buffers_indexed {
    /**
     * Enables blending for an individual draw buffer.
     *
     * @param target must be BLEND.
     * @param index is an integer i specifying the draw buffer associated with the symbolic constant DRAW_BUFFERi.
     */
    enablei: (target: number, index: number) => void;
    /**
     * Disables  blending for an individual draw buffer.
     *
     * @param target must be BLEND.
     * @param index is an integer i specifying the draw buffer associated with the symbolic constant DRAW_BUFFERi.
     */
    disablei: (buf: number, mode: number) => void;
    /**
     * The buf argument is an integer i that indicates that the blend equations should be modified for DRAW_BUFFERi.
     *
     * mode accepts the same tokens as mode in blendEquation.
     */
    blendEquationi: (target: number, index: number) => void;
    /**
     * The buf argument is an integer i that indicates that the blend equations should be modified for DRAW_BUFFERi.
     *
     * modeRGB and modeAlpha accept the same tokens as modeRGB and modeAlpha in blendEquationSeparate.
     */
    blendEquationSeparatei: (buf: number, modeRGB: number, modeAlpha: number) => void;
    /**
     * The buf argument is an integer i that indicates that the blend functions should be modified for DRAW_BUFFERi.
     *
     * src and dst accept the same tokens as src and dst in blendFunc.
     */
    blendFunci: (buf: number, src: number, dst: number) => void;
    /**
     * The buf argument is an integer i that indicates that the blend functions should be modified for DRAW_BUFFERi.
     *
     * srcRGB, dstRGB, srcAlpha, and dstAlpha accept the same tokens as srcRGB, dstRGB, srcAlpha, and dstAlpha parameters in blendEquationSeparate.
     */
    blendFuncSeparatei: (buf: number, srcRGB: number, dstRGB: number, srcAlpha: number, dstAlpha: number) => void;
    /**
     * The buf argument is an integer i that indicates that the write mask should be modified for DRAW_BUFFERi.
     *
     * r, g, b, and a indicate whether R, G, B, or A values, respectively, are written or not (a value of TRUE means that the corresponding value is written).
     */
    colorMaski: (buf: number, r: boolean, g: boolean, b: boolean, a: boolean) => void;
}

export function getDrawBuffersIndexed(gl: GLRenderingContext): COMPAT_draw_buffers_indexed | null {
    const ext = gl.getExtension('OES_draw_buffers_indexed');
    if (ext === null) return null;
    return {
        enablei: ext.enableiOES.bind(ext),
        disablei: ext.disableiOES.bind(ext),
        blendEquationi: ext.blendEquationiOES.bind(ext),
        blendEquationSeparatei: ext.blendEquationSeparateiOES.bind(ext),
        blendFunci: ext.blendFunciOES.bind(ext),
        blendFuncSeparatei: ext.blendFuncSeparateiOES.bind(ext),
        colorMaski: ext.colorMaskiOES.bind(ext),
    };
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_shader_texture_lod/
 */
export interface COMPAT_shader_texture_lod {
}

export function getShaderTextureLod(gl: GLRenderingContext): COMPAT_shader_texture_lod | null {
    return isWebGL2(gl) ? {} : gl.getExtension('EXT_shader_texture_lod');
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_depth_texture/
 */
export interface COMPAT_depth_texture {
    readonly UNSIGNED_INT_24_8: number;
}

export function getDepthTexture(gl: GLRenderingContext): COMPAT_depth_texture | null {
    if (isWebGL2(gl)) {
        return {
            UNSIGNED_INT_24_8: gl.UNSIGNED_INT_24_8
        };
    } else {
        const ext = gl.getExtension('WEBGL_depth_texture');
        if (ext === null) return null;
        return {
            UNSIGNED_INT_24_8: ext.UNSIGNED_INT_24_8_WEBGL
        };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_sRGB/
 */
export interface COMPAT_sRGB {
    readonly FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING: number;
    readonly SRGB8_ALPHA8: number;
    readonly SRGB8: number;
    readonly SRGB: number;
}

export function getSRGB(gl: GLRenderingContext): COMPAT_sRGB | null {
    if (isWebGL2(gl)) {
        return {
            FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING: gl.FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING,
            SRGB8_ALPHA8: gl.SRGB8_ALPHA8,
            SRGB8: gl.SRGB8,
            SRGB: gl.SRGB
        };
    } else {
        const ext = gl.getExtension('EXT_sRGB');
        if (ext === null) return null;
        return {
            FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING: ext.FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING_EXT,
            SRGB8_ALPHA8: ext.SRGB8_ALPHA8_EXT,
            SRGB8: ext.SRGB_ALPHA_EXT,
            SRGB: ext.SRGB_EXT
        };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_disjoint_timer_query/ and https://registry.khronos.org/webgl/extensions/EXT_disjoint_timer_query_webgl2/
 */
export interface COMPAT_disjoint_timer_query {
    /** A GLint indicating the number of bits used to hold the query result for the given target. */
    QUERY_COUNTER_BITS: number
    /** A WebGLQuery object, which is the currently active query for the given target. */
    CURRENT_QUERY: number
    /** A GLuint64EXT containing the query result. */
    QUERY_RESULT: number
    /** A GLboolean indicating whether or not a query result is available. */
    QUERY_RESULT_AVAILABLE: number
    /** Elapsed time (in nanoseconds). */
    TIME_ELAPSED: number
    /** The current time. */
    TIMESTAMP: number
    /** A GLboolean indicating whether or not the GPU performed any disjoint operation. */
    GPU_DISJOINT: number

    /** Creates a new WebGLTimerQueryEXT. */
    createQuery: () => WebGLQuery
    /** Deletes a given WebGLTimerQueryEXT. */
    deleteQuery: (query: WebGLQuery) => void
    /** Returns true if a given object is a valid WebGLTimerQueryEXT. */
    isQuery: (query: WebGLQuery) => boolean
    /** The timer starts when all commands prior to beginQueryEXT have been fully executed. */
    beginQuery: (target: number, query: WebGLQuery) => void
    /** The timer stops when all commands prior to endQueryEXT have been fully executed. */
    endQuery: (target: number) => void
    /** Records the current time into the corresponding query object. */
    queryCounter: (query: WebGLQuery, target: number) => void
    /** Returns information about a query target. */
    getQuery: (target: number, pname: number) => WebGLQuery | number
    /** Return the state of a query object. */
    getQueryParameter: (query: WebGLQuery, pname: number) => number | boolean
}

export function getDisjointTimerQuery(gl: GLRenderingContext): COMPAT_disjoint_timer_query | null {
    if (isWebGL2(gl)) {
        // Firefox has EXT_disjoint_timer_query in webgl2
        const ext = gl.getExtension('EXT_disjoint_timer_query_webgl2') || gl.getExtension('EXT_disjoint_timer_query');
        if (ext === null) return null;
        return {
            QUERY_COUNTER_BITS: ext.QUERY_COUNTER_BITS_EXT,
            CURRENT_QUERY: gl.CURRENT_QUERY,
            QUERY_RESULT: gl.QUERY_RESULT,
            QUERY_RESULT_AVAILABLE: gl.QUERY_RESULT_AVAILABLE,
            TIME_ELAPSED: ext.TIME_ELAPSED_EXT,
            TIMESTAMP: ext.TIMESTAMP_EXT,
            GPU_DISJOINT: ext.GPU_DISJOINT_EXT,

            createQuery: gl.createQuery.bind(gl),
            deleteQuery: gl.deleteQuery.bind(gl),
            isQuery: gl.isQuery.bind(gl),
            beginQuery: gl.beginQuery.bind(gl),
            endQuery: gl.endQuery.bind(gl),
            queryCounter: ext.queryCounterEXT.bind(ext),
            getQuery: gl.getQuery.bind(gl),
            getQueryParameter: gl.getQueryParameter.bind(gl),
        };
    } else {
        const ext = gl.getExtension('EXT_disjoint_timer_query');
        if (ext === null) return null;
        return {
            QUERY_COUNTER_BITS: ext.QUERY_COUNTER_BITS_EXT,
            CURRENT_QUERY: ext.CURRENT_QUERY_EXT,
            QUERY_RESULT: ext.QUERY_RESULT_EXT,
            QUERY_RESULT_AVAILABLE: ext.QUERY_RESULT_AVAILABLE_EXT,
            TIME_ELAPSED: ext.TIME_ELAPSED_EXT,
            TIMESTAMP: ext.TIMESTAMP_EXT,
            GPU_DISJOINT: ext.GPU_DISJOINT_EXT,

            createQuery: ext.createQueryEXT.bind(ext),
            deleteQuery: ext.deleteQueryEXT.bind(ext),
            isQuery: ext.isQueryEXT.bind(ext),
            beginQuery: ext.beginQueryEXT.bind(ext),
            endQuery: ext.endQueryEXT.bind(ext),
            queryCounter: ext.queryCounterEXT.bind(ext),
            getQuery: ext.getQueryEXT.bind(ext),
            getQueryParameter: ext.getQueryObjectEXT.bind(ext),
        };
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_multi_draw/
 */
export interface COMPAT_multi_draw {
    /**
     * Renders multiple primitives from array data. It is identical to multiple calls to the `drawArrays` method.
     */
    readonly multiDrawArrays: (mode: number, firstsList: Int32Array, firstsOffset: number, countsList: Int32Array, countsOffset: number, drawcount: number) => void;
    /**
     * Renders multiple primitives from array data. It is identical to multiple calls to the `drawElements` method.
     */
    readonly multiDrawElements: (mode: number, countsList: Int32Array, countsOffset: number, type: number, offsetsList: Int32Array, offsetsOffset: number, drawcount: number) => void;
    /**
     * Renders multiple primitives from array data. It is identical to multiple calls to the `drawArraysInstanced` method.
     */
    readonly multiDrawArraysInstanced: (mode: number, firstsList: Int32Array, firstsOffset: number, countsList: Int32Array, countsOffset: number, instanceCountsList: Int32Array, instanceCountsOffset: number, drawcount: number) => void;
    /**
     * Renders multiple primitives from array data. It is identical to multiple calls to the `drawElementsInstanced` method.
     */
    readonly multiDrawElementsInstanced: (mode: number, countsList: Int32Array, countsOffset: number, type: number, offsetsList: Int32Array, offsetsOffset: number, instanceCountsList: Int32Array, instanceCountsOffset: number, drawcount: number) => void;
}

export function getMultiDraw(gl: GLRenderingContext): COMPAT_multi_draw | null {
    const ext = gl.getExtension('WEBGL_multi_draw');
    if (ext) {
        return {
            multiDrawArrays: ext.multiDrawArraysWEBGL.bind(ext),
            multiDrawElements: ext.multiDrawElementsWEBGL.bind(ext),
            multiDrawArraysInstanced: ext.multiDrawArraysInstancedWEBGL.bind(ext),
            multiDrawElementsInstanced: ext.multiDrawElementsInstancedWEBGL.bind(ext),
        };
    } else {
        return null;
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_draw_instanced_base_vertex_base_instance/
 */
export interface COMPAT_draw_instanced_base_vertex_base_instance {
    /**
     * Behaves identically to DrawArraysInstanced except that `baseInstance` is passed down to DrawArraysOneInstance instead of zero.
     */
    readonly drawArraysInstancedBaseInstance: (mode: number, first: number, count: number,
        instanceCount: number, baseInstance: number) => void;
    /**
     * Behaves identically to DrawElementsInstanced except that `baseVertex` and `baseInstance` are passed down to DrawElementsOneInstance instead of zero.
     */
    readonly drawElementsInstancedBaseVertexBaseInstance: (mode: number, count: number, type: number, offset: number, instanceCount: number, baseVertex: number, baseInstance: number) => void;
}

export function getDrawInstancedBaseVertexBaseInstance(gl: GLRenderingContext): COMPAT_draw_instanced_base_vertex_base_instance | null {
    const ext = gl.getExtension('WEBGL_draw_instanced_base_vertex_base_instance');
    if (ext) {
        return {
            drawArraysInstancedBaseInstance: ext.drawArraysInstancedBaseInstanceWEBGL.bind(ext),
            drawElementsInstancedBaseVertexBaseInstance: ext.drawElementsInstancedBaseVertexBaseInstanceWEBGL.bind(ext),
        };
    } else {
        return null;
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_multi_draw_instanced_base_vertex_base_instance/
 */
export interface COMPAT_multi_draw_instanced_base_vertex_base_instance {
    /**
     * Behaves identically to DrawArraysInstancedBaseInstance except that a list of arrays is specified instead. The number of lists is specified in the `drawcount` parameter.
     */
    readonly multiDrawArraysInstancedBaseInstance: (mode: number, firstsList: Int32Array, firstsOffset: number, countsList: Int32Array, countsOffset: number, instanceCountsList: Int32Array, instanceCountsOffset: number, baseInstancesList: Uint32Array, baseInstancesOffset: number, drawcount: number) => void;
    /**
     * Behaves identically to DrawElementsInstancedBaseVertexBaseInstance except that a list of arrays is specified instead. The number of lists is specified in the `drawcount` parameter.
     */
    readonly multiDrawElementsInstancedBaseVertexBaseInstance: (mode: number, countsList: Int32Array, countsOffset: number, type: number, offsetsList: Int32Array, offsetsOffset: number, instanceCountsList: Int32Array, instanceCountsOffset: number, baseVerticesList: Int32Array, baseVerticesOffset: number, baseInstancesList: Uint32Array, baseInstancesOffset: number, drawcount: number) => void;
}

export function getMultiDrawInstancedBaseVertexBaseInstance(gl: GLRenderingContext): COMPAT_multi_draw_instanced_base_vertex_base_instance | null {
    const ext = gl.getExtension('WEBGL_multi_draw_instanced_base_vertex_base_instance');
    if (ext) {
        return {
            multiDrawArraysInstancedBaseInstance: ext.multiDrawArraysInstancedBaseInstanceWEBGL.bind(ext),
            multiDrawElementsInstancedBaseVertexBaseInstance: ext.multiDrawElementsInstancedBaseVertexBaseInstanceWEBGL.bind(ext),
        };
    } else {
        return null;
    }
}

/**
 * See https://registry.khronos.org/webgl/extensions/KHR_parallel_shader_compile/
 */
export interface COMPAT_parallel_shader_compile {
    readonly COMPLETION_STATUS: number;
}

export function getParallelShaderCompile(gl: GLRenderingContext): COMPAT_parallel_shader_compile | null {
    const ext = gl.getExtension('KHR_parallel_shader_compile');
    if (ext === null) return null;
    return {
        COMPLETION_STATUS: ext.COMPLETION_STATUS_KHR,
    };
}

/**
 * See https://registry.khronos.org/webgl/extensions/OES_fbo_render_mipmap/
 */
export interface COMPAT_fboRenderMipmap {
}

export function getFboRenderMipmap(gl: GLRenderingContext): COMPAT_fboRenderMipmap | null {
    return isWebGL2(gl) ? {} : gl.getExtension('OES_fbo_render_mipmap');
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_provoking_vertex/
 */
export interface COMPAT_provoking_vertex {
    readonly FIRST_VERTEX_CONVENTION: number;
    readonly LAST_VERTEX_CONVENTION: number;
    readonly PROVOKING_VERTEX: number;
    provokingVertex(provokeMode: number): void;
}

export function getProvokingVertex(gl: GLRenderingContext): COMPAT_provoking_vertex | null {
    if (isWebGL2(gl)) {
        const ext = gl.getExtension('WEBGL_provoking_vertex');
        if (ext) {
            return {
                FIRST_VERTEX_CONVENTION: ext.FIRST_VERTEX_CONVENTION_WEBGL,
                LAST_VERTEX_CONVENTION: ext.LAST_VERTEX_CONVENTION_WEBGL,
                PROVOKING_VERTEX: ext.PROVOKING_VERTEX_WEBGL,
                provokingVertex: ext.provokingVertexWEBGL.bind(ext)
            };
        }
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_clip_cull_distance/
 */
export interface COMPAT_clip_cull_distance {
    readonly MAX_CLIP_DISTANCES: number;
    readonly MAX_CULL_DISTANCES: number;
    readonly MAX_COMBINED_CLIP_AND_CULL_DISTANCES: number;

    readonly CLIP_DISTANCE0: number;
    readonly CLIP_DISTANCE1: number;
    readonly CLIP_DISTANCE2: number;
    readonly CLIP_DISTANCE3: number;
    readonly CLIP_DISTANCE4: number;
    readonly CLIP_DISTANCE5: number;
    readonly CLIP_DISTANCE6: number;
    readonly CLIP_DISTANCE7: number;
}

export function getClipCullDistance(gl: GLRenderingContext): COMPAT_clip_cull_distance | null {
    if (isWebGL2(gl)) {
        const ext = gl.getExtension('WEBGL_clip_cull_distance');
        if (ext) {
            return {
                MAX_CLIP_DISTANCES: ext.MAX_CLIP_DISTANCES_WEBGL,
                MAX_CULL_DISTANCES: ext.MAX_CULL_DISTANCES_WEBGL,
                MAX_COMBINED_CLIP_AND_CULL_DISTANCES: ext.MAX_COMBINED_CLIP_AND_CULL_DISTANCES_WEBGL,

                CLIP_DISTANCE0: ext.CLIP_DISTANCE0_WEBGL,
                CLIP_DISTANCE1: ext.CLIP_DISTANCE1_WEBGL,
                CLIP_DISTANCE2: ext.CLIP_DISTANCE2_WEBGL,
                CLIP_DISTANCE3: ext.CLIP_DISTANCE3_WEBGL,
                CLIP_DISTANCE4: ext.CLIP_DISTANCE4_WEBGL,
                CLIP_DISTANCE5: ext.CLIP_DISTANCE5_WEBGL,
                CLIP_DISTANCE6: ext.CLIP_DISTANCE6_WEBGL,
                CLIP_DISTANCE7: ext.CLIP_DISTANCE7_WEBGL
            };
        }
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_conservative_depth/
 */
export interface COMPAT_conservative_depth {
}

export function getConservativeDepth(gl: GLRenderingContext): COMPAT_conservative_depth | null {
    if (isWebGL2(gl)) {
        const ext = gl.getExtension('EXT_conservative_depth');
        if (ext) {
            return {};
        }
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_stencil_texturing/
 */
export interface COMPAT_stencil_texturing {
    readonly DEPTH_STENCIL_TEXTURE_MODE: number;
    readonly STENCIL_INDEX: number;
}

export function getStencilTexturing(gl: GLRenderingContext): COMPAT_stencil_texturing | null {
    if (isWebGL2(gl)) {
        const ext = gl.getExtension('WEBGL_stencil_texturing');
        if (ext) {
            return {
                DEPTH_STENCIL_TEXTURE_MODE: ext.DEPTH_STENCIL_TEXTURE_MODE_WEBGL,
                STENCIL_INDEX: ext.STENCIL_INDEX_WEBGL
            };
        }
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_clip_control/
 */
export interface COMPAT_clip_control {
    readonly LOWER_LEFT: number;
    readonly UPPER_LEFT: number;

    readonly NEGATIVE_ONE_TO_ONE: number;
    readonly ZERO_TO_ONE: number;

    readonly CLIP_ORIGIN: number;
    readonly CLIP_DEPTH_MODE: number;

    /**
     * @param origin must be LOWER_LEFT (default) or UPPER_LEFT.
     * @param depth must be NEGATIVE_ONE_TO_ONE (default) or ZERO_TO_ONE.
     */
    clipControl(origin: number, depth: number): void
}

export function getClipControl(gl: GLRenderingContext): COMPAT_clip_control | null {
    const ext = gl.getExtension('EXT_clip_control');
    if (ext) {
        return {
            LOWER_LEFT: ext.LOWER_LEFT_EXT,
            UPPER_LEFT: ext.UPPER_LEFT_EXT,

            NEGATIVE_ONE_TO_ONE: ext.NEGATIVE_ONE_TO_ONE_EXT,
            ZERO_TO_ONE: ext.ZERO_TO_ONE_EXT,

            CLIP_ORIGIN: ext.CLIP_ORIGIN_EXT,
            CLIP_DEPTH_MODE: ext.CLIP_DEPTH_MODE_EXT,

            clipControl: ext.clipControlEXT.bind(ext)
        };
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_render_snorm/
 */
export interface COMPAT_render_snorm {
}

export function getRenderSnorm(gl: GLRenderingContext): COMPAT_render_snorm | null {
    if (isWebGL2(gl)) {
        const ext = gl.getExtension('EXT_render_snorm');
        if (ext) {
            return {};
        }
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/WEBGL_render_shared_exponent/
 */
export interface COMPAT_render_shared_exponent {
}

export function getRenderSharedExponent(gl: GLRenderingContext): COMPAT_render_shared_exponent | null {
    if (isWebGL2(gl)) {
        const ext = gl.getExtension('WEBGL_render_shared_exponent');
        if (ext) {
            return {};
        }
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_texture_norm16/
 */
export interface COMPAT_texture_norm16 {
    readonly R16: number;
    readonly RG16: number;
    readonly RGB16: number;
    readonly RGBA16: number;
    readonly R16_SNORM: number;
    readonly RG16_SNORM: number;
    readonly RGB16_SNORM: number;
    readonly RGBA16_SNORM: number;
}

export function getTextureNorm16(gl: GLRenderingContext): COMPAT_texture_norm16 | null {
    const ext = gl.getExtension('EXT_texture_norm16');
    if (ext) {
        return {
            R16: ext.R16_EXT,
            RG16: ext.RG16_EXT,
            RGB16: ext.RGB16_EXT,
            RGBA16: ext.RGBA16_EXT,
            R16_SNORM: ext.R16_SNORM_EXT,
            RG16_SNORM: ext.RG16_SNORM_EXT,
            RGB16_SNORM: ext.RGB16_SNORM_EXT,
            RGBA16_SNORM: ext.RGBA16_SNORM_EXT
        };
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/EXT_depth_clamp/
 */
export interface COMPAT_depth_clamp {
    readonly DEPTH_CLAMP: number;
}

export function getDepthClamp(gl: GLRenderingContext): COMPAT_depth_clamp | null {
    const ext = gl.getExtension('EXT_depth_clamp');
    if (ext) {
        return {
            DEPTH_CLAMP: ext.DEPTH_CLAMP_EXT
        };
    }
    return null;
}

/**
 * See https://registry.khronos.org/webgl/extensions/OVR_multiview2/
 */
export interface COMPAT_multiview2 {
    readonly FRAMEBUFFER_ATTACHMENT_TEXTURE_NUM_VIEWS: number;
    readonly FRAMEBUFFER_ATTACHMENT_TEXTURE_BASE_VIEW_INDEX: number;
    readonly MAX_VIEWS: number;
    readonly FRAMEBUFFER_INCOMPLETE_VIEW_TARGETS: number;

    /**
     * Binds a texture to a framebuffer attachment point for multiview rendering.
     *
     * @param target must be FRAMEBUFFER.
     * @param attachment must be FRAMEBUFFER_ATTACHMENT_TEXTURE_LAYER_OVR.
     * @param texture is the texture to bind.
     * @param level is the mipmap level of the texture to bind.
     * @param baseViewIndex is the index of the first view in the multiview texture.
     * @param numViews is the number of views in the multiview texture.
     */
    framebufferTextureMultiview(target: number, attachment: number, texture: number, level: number, baseViewIndex: number, numViews: number): void;
}

export function getMultiview2(gl: GLRenderingContext): COMPAT_multiview2 | null {
    if (isWebGL2(gl)) {
        const ext = gl.getExtension('OVR_multiview2');
        if (ext) {
            return {
                FRAMEBUFFER_ATTACHMENT_TEXTURE_NUM_VIEWS: ext.FRAMEBUFFER_ATTACHMENT_TEXTURE_NUM_VIEWS_OVR,
                FRAMEBUFFER_ATTACHMENT_TEXTURE_BASE_VIEW_INDEX: ext.FRAMEBUFFER_ATTACHMENT_TEXTURE_BASE_VIEW_INDEX_OVR,
                MAX_VIEWS: ext.MAX_VIEWS_OVR,
                FRAMEBUFFER_INCOMPLETE_VIEW_TARGETS: ext.FRAMEBUFFER_INCOMPLETE_VIEW_TARGETS_OVR,

                framebufferTextureMultiview: ext.framebufferTextureMultiviewOVR.bind(ext)
            };
        }
    }
    return null;
}

export function getNoNonInstancedActiveAttribs(gl: GLRenderingContext): boolean {
    if (!isWebGL2(gl)) return false;

    if (typeof navigator !== 'undefined') {
        const ffMatch = window.navigator.userAgent.match(/Firefox\/([0-9]+)\./);
        if (!ffMatch) return true;

        const ffVersion = parseInt(ffMatch[1]);
        // supported since FF 85 (https://bugzilla.mozilla.org/show_bug.cgi?id=1679693)
        return ffVersion >= 85;
    }
    return false;
}

//

const TextureTestVertShader = `
attribute vec4 aPosition;

void main() {
    gl_Position = aPosition;
}`;

const TextureTestFragShader = `
precision mediump float;
uniform vec4 uColor;
uniform sampler2D uTexture;

void main() {
    gl_FragColor = texture2D(uTexture, vec2(0.5, 0.5)) * uColor;
}`;

const TextureTestTexCoords = new Float32Array([
    -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0
]);

// adapted from https://stackoverflow.com/questions/28827511/
export function testColorBuffer(gl: GLRenderingContext, type: number) {
    // setup shaders
    const vertShader = getShader(gl, { type: 'vert', source: TextureTestVertShader });
    const fragShader = getShader(gl, { type: 'frag', source: TextureTestFragShader });
    if (!vertShader || !fragShader) return false;

    // setup program
    const program = getProgram(gl);
    gl.attachShader(program, vertShader);
    gl.attachShader(program, fragShader);
    gl.linkProgram(program);
    gl.useProgram(program);

    // look up where the vertex data needs to go.
    const positionLocation = gl.getAttribLocation(program, 'aPosition');
    const colorLoc = gl.getUniformLocation(program, 'uColor');
    if (!colorLoc) {
        if (isDebugMode) {
            console.log(`error getting 'uColor' uniform location`);
        }
        return false;
    }

    // provide texture coordinates for the rectangle.
    const positionBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, TextureTestTexCoords, gl.STATIC_DRAW);
    gl.enableVertexAttribArray(positionLocation);
    gl.vertexAttribPointer(positionLocation, 2, gl.FLOAT, false, 0, 0);

    const whiteTex = gl.createTexture();
    const whiteData = new Uint8Array([255, 255, 255, 255]);
    gl.bindTexture(gl.TEXTURE_2D, whiteTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, whiteData);

    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, type, null);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

    const fb = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    if (status !== gl.FRAMEBUFFER_COMPLETE) {
        if (isDebugMode) {
            console.log(`error creating framebuffer for '${type}'`);
        }
        return false;
    }

    // Draw the rectangle.
    gl.bindTexture(gl.TEXTURE_2D, whiteTex);
    gl.uniform4fv(colorLoc, [0, 10, 20, 1]);
    gl.drawArrays(gl.TRIANGLES, 0, 6);

    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.clearColor(1, 0, 0, 1);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.uniform4fv(colorLoc, [0, 1 / 10, 1 / 20, 1]);
    gl.drawArrays(gl.TRIANGLES, 0, 6);

    // Check if rendered correctly
    const pixel = new Uint8Array(4);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, pixel);
    if (pixel[0] !== 0 || pixel[1] < 248 || pixel[2] < 248 || pixel[3] < 254) {
        if (isDebugMode) {
            console.log(`not able to actually render to '${type}' texture`);
        }
        return false;
    }

    // Check reading from float texture
    if (type === gl.FLOAT) {
        gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
        const floatPixel = new Float32Array(4);
        gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, floatPixel);
        const error = gl.getError();
        if (error) {
            if (isDebugMode) {
                console.log(`error reading float pixels: '${getErrorDescription(gl, error)}'`);
            }
            return false;
        }
    }

    return true;
}
