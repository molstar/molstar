/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GLRenderingContext, COMPAT_instanced_arrays, COMPAT_standard_derivatives, COMPAT_vertex_array_object, getInstancedArrays, getStandardDerivatives, COMPAT_element_index_uint, getElementIndexUint, COMPAT_texture_float, getTextureFloat, COMPAT_texture_float_linear, getTextureFloatLinear, COMPAT_blend_minmax, getBlendMinMax, getFragDepth, COMPAT_frag_depth, COMPAT_color_buffer_float, getColorBufferFloat, COMPAT_draw_buffers, getDrawBuffers, getShaderTextureLod, COMPAT_shader_texture_lod, getDepthTexture, COMPAT_depth_texture, COMPAT_sRGB, getSRGB, getTextureHalfFloat, getTextureHalfFloatLinear, COMPAT_texture_half_float, COMPAT_texture_half_float_linear, COMPAT_color_buffer_half_float, getColorBufferHalfFloat, getVertexArrayObject } from './compat';
import { isDebugMode } from '../../mol-util/debug';

export type WebGLExtensions = {
    instancedArrays: COMPAT_instanced_arrays
    elementIndexUint: COMPAT_element_index_uint

    standardDerivatives: COMPAT_standard_derivatives | null
    textureFloat: COMPAT_texture_float | null
    textureFloatLinear: COMPAT_texture_float_linear | null
    textureHalfFloat: COMPAT_texture_half_float | null
    textureHalfFloatLinear: COMPAT_texture_half_float_linear | null
    depthTexture: COMPAT_depth_texture | null
    blendMinMax: COMPAT_blend_minmax | null
    vertexArrayObject: COMPAT_vertex_array_object | null
    fragDepth: COMPAT_frag_depth | null
    colorBufferFloat: COMPAT_color_buffer_float | null
    colorBufferHalfFloat: COMPAT_color_buffer_half_float | null
    drawBuffers: COMPAT_draw_buffers | null
    shaderTextureLod: COMPAT_shader_texture_lod | null
    sRGB: COMPAT_sRGB | null
}

export function createExtensions(gl: GLRenderingContext): WebGLExtensions {
    const instancedArrays = getInstancedArrays(gl);
    if (instancedArrays === null) {
        throw new Error('Could not find support for "instanced_arrays"');
    }
    const elementIndexUint = getElementIndexUint(gl);
    if (elementIndexUint === null) {
        throw new Error('Could not find support for "element_index_uint"');
    }
    const standardDerivatives = getStandardDerivatives(gl);
    if (standardDerivatives === null) {
        throw new Error('Could not find support for "standard_derivatives"');
    }

    const textureFloat = getTextureFloat(gl);
    if (isDebugMode && textureFloat === null) {
        console.log('Could not find support for "texture_float"');
    }
    const textureFloatLinear = getTextureFloatLinear(gl);
    if (isDebugMode && textureFloatLinear === null) {
        // TODO handle non-support downstream (no gpu gaussian calc, no gpu mc???)
        // - can't be a required extension because it is not supported by `headless-gl`
        console.log('Could not find support for "texture_float_linear"');
    }
    const textureHalfFloat = getTextureHalfFloat(gl);
    if (isDebugMode && textureHalfFloat === null) {
        console.log('Could not find support for "texture_half_float"');
    }
    const textureHalfFloatLinear = getTextureHalfFloatLinear(gl);
    if (isDebugMode && textureHalfFloatLinear === null) {
        // TODO handle non-support downstream (no gpu gaussian calc, no gpu mc???)
        // - can't be a required extension because it is not supported by `headless-gl`
        console.log('Could not find support for "texture_half_float_linear"');
    }
    const depthTexture = getDepthTexture(gl);
    if (isDebugMode && depthTexture === null) {
        console.log('Could not find support for "depth_texture"');
    }
    const blendMinMax = getBlendMinMax(gl);
    if (isDebugMode && blendMinMax === null) {
        // TODO handle non-support downstream (e.g. no gpu gaussian calc)
        // - can't be a required extension because it is not supported by `headless-gl`
        console.log('Could not find support for "blend_minmax"');
    }
    const vertexArrayObject = getVertexArrayObject(gl);
    if (isDebugMode && vertexArrayObject === null) {
        console.log('Could not find support for "vertex_array_object"');
    }
    const fragDepth = getFragDepth(gl);
    if (isDebugMode && fragDepth === null) {
        console.log('Could not find support for "frag_depth"');
    }
    const colorBufferFloat = getColorBufferFloat(gl);
    if (isDebugMode && colorBufferFloat === null) {
        console.log('Could not find support for "color_buffer_float"');
    }
    const colorBufferHalfFloat = getColorBufferHalfFloat(gl);
    if (isDebugMode && colorBufferHalfFloat === null) {
        console.log('Could not find support for "color_buffer_half_float"');
    }
    const drawBuffers = getDrawBuffers(gl);
    if (isDebugMode && drawBuffers === null) {
        console.log('Could not find support for "draw_buffers"');
    }
    const shaderTextureLod = getShaderTextureLod(gl);
    if (isDebugMode && shaderTextureLod === null) {
        console.log('Could not find support for "shader_texture_lod"');
    }
    const sRGB = getSRGB(gl);
    if (isDebugMode && sRGB === null) {
        console.log('Could not find support for "sRGB"');
    }

    return {
        instancedArrays,
        standardDerivatives,
        textureFloat,
        textureFloatLinear,
        textureHalfFloat,
        textureHalfFloatLinear,
        elementIndexUint,
        depthTexture,

        blendMinMax,
        vertexArrayObject,
        fragDepth,
        colorBufferFloat,
        colorBufferHalfFloat,
        drawBuffers,
        shaderTextureLod,
        sRGB,
    };
}