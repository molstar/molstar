/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GLRenderingContext, COMPAT_instanced_arrays, COMPAT_standard_derivatives, COMPAT_vertex_array_object, getInstancedArrays, getStandardDerivatives, getVertexArrayObject, COMPAT_element_index_uint, getElementIndexUint, COMPAT_texture_float, getTextureFloat, COMPAT_texture_float_linear, getTextureFloatLinear, COMPAT_blend_minmax, getBlendMinMax, getFragDepth, COMPAT_frag_depth, COMPAT_color_buffer_float, getColorBufferFloat, COMPAT_draw_buffers, getDrawBuffers, getShaderTextureLod, COMPAT_shader_texture_lod, getDepthTexture, COMPAT_depth_texture } from './compat';

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

export function createExtensions(gl: GLRenderingContext): WebGLExtensions {
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