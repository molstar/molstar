/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

export interface COMPAT_texture_half_float_linear {
}

export function getTextureHalfFloatLinear(gl: GLRenderingContext): COMPAT_texture_half_float_linear | null {
    return gl.getExtension('OES_texture_half_float_linear');
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
