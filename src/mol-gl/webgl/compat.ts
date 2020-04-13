/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

export type GLRenderingContext = WebGLRenderingContext | WebGL2RenderingContext

export function isWebGL(gl: any): gl is WebGLRenderingContext {
    return typeof WebGLRenderingContext !== 'undefined' && gl instanceof WebGLRenderingContext;
}

export function isWebGL2(gl: any): gl is WebGL2RenderingContext {
    return typeof WebGL2RenderingContext !== 'undefined' && gl instanceof WebGL2RenderingContext;
}

export interface COMPAT_instanced_arrays {
    drawArraysInstanced(mode: number, first: number, count: number, primcount: number): void;
    drawElementsInstanced(mode: number, count: number, type: number, offset: number, primcount: number): void;
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

export interface COMPAT_element_index_uint {
}

export function getElementIndexUint(gl: GLRenderingContext): COMPAT_element_index_uint | null {
    return isWebGL2(gl) ? {} : gl.getExtension('OES_element_index_uint');
}

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

export interface COMPAT_texture_float {
}

export function getTextureFloat(gl: GLRenderingContext): COMPAT_texture_float | null {
    return isWebGL2(gl) ? {} : gl.getExtension('OES_texture_float');
}

export interface COMPAT_texture_float_linear {
}

export function getTextureFloatLinear(gl: GLRenderingContext): COMPAT_texture_float_linear | null {
    return gl.getExtension('OES_texture_float_linear');
}

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

export interface COMPAT_frag_depth {
}

export function getFragDepth(gl: GLRenderingContext): COMPAT_frag_depth | null {
    return isWebGL2(gl) ? {} : gl.getExtension('EXT_frag_depth');
}

export interface COMPAT_color_buffer_float {
    readonly RGBA32F: number;
}

export function getColorBufferFloat(gl: GLRenderingContext): COMPAT_color_buffer_float | null {
    if (isWebGL2(gl)) {
        if (gl.getExtension('EXT_color_buffer_float') === null) return null;
        return { RGBA32F: gl.RGBA32F };
    } else {
        const ext = gl.getExtension('WEBGL_color_buffer_float');
        if (ext === null) return null;
        return { RGBA32F: ext.RGBA32F_EXT };
    }
}

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

export interface COMPAT_shader_texture_lod {
}

export function getShaderTextureLod(gl: GLRenderingContext): COMPAT_shader_texture_lod | null {
    return isWebGL2(gl) ? {} : gl.getExtension('EXT_shader_texture_lod');
}

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