/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GLRenderingContext, COMPAT_instanced_arrays, COMPAT_standard_derivatives, COMPAT_vertex_array_object, getInstancedArrays, getStandardDerivatives, getVertexArrayObject, COMPAT_element_index_uint, getElementIndexUint, COMPAT_texture_float, getTextureFloat, COMPAT_texture_float_linear, getTextureFloatLinear, COMPAT_blend_minmax, getBlendMinMax, getFragDepth, COMPAT_frag_depth, COMPAT_color_buffer_float, getColorBufferFloat, COMPAT_draw_buffers, getDrawBuffers, getShaderTextureLod, COMPAT_shader_texture_lod, getDepthTexture, COMPAT_depth_texture } from './compat';

export type WebGLExtensions = {
    instancedArrays: COMPAT_instanced_arrays
    textureFloat: COMPAT_texture_float
    elementIndexUint: COMPAT_element_index_uint

    standardDerivatives: COMPAT_standard_derivatives | null
    textureFloatLinear: COMPAT_texture_float_linear | null
    depthTexture: COMPAT_depth_texture | null
    blendMinMax: COMPAT_blend_minmax | null
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
    const textureFloat = getTextureFloat(gl)
    if (textureFloat === null) {
        throw new Error('Could not find support for "texture_float"')
    }
    const elementIndexUint = getElementIndexUint(gl)
    if (elementIndexUint === null) {
        throw new Error('Could not find support for "element_index_uint"')
    }

    const standardDerivatives = getStandardDerivatives(gl)
    if (standardDerivatives === null) {
        // TODO handle non-support downstream (e.g. no flat shading)
        // throw new Error('Could not find support for "standard_derivatives"')
        console.log('Could not find support for "standard_derivatives"')
    }
    const textureFloatLinear = getTextureFloatLinear(gl)
    if (textureFloatLinear === null) {
        // TODO handle non-support downstream (no gpu gaussian calc, no gpu mc???)
        // throw new Error('Could not find support for "texture_float_linear"')
        console.log('Could not find support for "texture_float_linear"')
    }
    const depthTexture = getDepthTexture(gl)
    if (depthTexture === null) {
        console.log('Could not find support for "depth_texture"')
    }
    const blendMinMax = getBlendMinMax(gl)
    if (blendMinMax === null) {
        console.log('Could not find support for "blend_minmax"')
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
        textureFloat,
        textureFloatLinear,
        elementIndexUint,
        depthTexture,

        blendMinMax,
        vertexArrayObject,
        fragDepth,
        colorBufferFloat,
        drawBuffers,
        shaderTextureLod,
    }
}