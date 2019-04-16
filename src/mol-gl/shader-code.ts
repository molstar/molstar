/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { idFactory } from 'mol-util/id-factory';
import { WebGLExtensions } from './webgl/context';
import { isWebGL2, GLRenderingContext } from './webgl/compat';

export type DefineKind = 'boolean' | 'string' | 'number'
export type DefineType = boolean | string
export type DefineValues = { [k: string]: ValueCell<DefineType> }

const shaderCodeId = idFactory()

export interface ShaderExtensions {
    readonly standardDerivatives?: boolean
    readonly fragDepth?: boolean
    readonly drawBuffers?: boolean
    readonly shaderTextureLod?: boolean
}

export interface ShaderCode {
    readonly id: number
    readonly vert: string
    readonly frag: string
    readonly extensions: ShaderExtensions
}

export function ShaderCode(vert: string, frag: string, extensions: ShaderExtensions = {}): ShaderCode {
    return { id: shaderCodeId(), vert, frag, extensions }
}

export const PointsShaderCode = ShaderCode(
    require('mol-gl/shader/points.vert').default,
    require('mol-gl/shader/points.frag').default
)

export const SpheresShaderCode = ShaderCode(
    require('mol-gl/shader/spheres.vert').default,
    require('mol-gl/shader/spheres.frag').default,
    { fragDepth: true }
)

export const TextShaderCode = ShaderCode(
    require('mol-gl/shader/text.vert').default,
    require('mol-gl/shader/text.frag').default,
    { standardDerivatives: true }
)

export const LinesShaderCode = ShaderCode(
    require('mol-gl/shader/lines.vert').default,
    require('mol-gl/shader/lines.frag').default
)

export const MeshShaderCode = ShaderCode(
    require('mol-gl/shader/mesh.vert').default,
    require('mol-gl/shader/mesh.frag').default,
    { standardDerivatives: true }
)

export const DirectVolumeShaderCode = ShaderCode(
    require('mol-gl/shader/direct-volume.vert').default,
    require('mol-gl/shader/direct-volume.frag').default,
    { fragDepth: true }
)



export type ShaderDefines = {
    [k: string]: ValueCell<DefineType>
}

function getDefinesCode (defines: ShaderDefines) {
    if (defines === undefined) return ''
    const lines = []
    for (const name in defines) {
        const define = defines[name]
        const v = define.ref.value
        if (v !== undefined) {
            if (typeof v === 'string') {
                lines.push(`#define ${name}_${v}`)
            } else if (typeof v === 'number') {
                lines.push(`#define ${name} ${v}`)
            } else if (typeof v === 'boolean') {
                if (v) lines.push(`#define ${name}`)
            } else {
                throw new Error('unknown define type')
            }
        }
    }
    return lines.join('\n') + '\n'
}

function getGlsl100FragPrefix(extensions: WebGLExtensions, shaderExtensions: ShaderExtensions) {
    const prefix: string[] = []
    if (shaderExtensions.standardDerivatives) {
        prefix.push('#extension GL_OES_standard_derivatives : enable')
        prefix.push('#define enabledStandardDerivatives')
    }
    if (shaderExtensions.fragDepth) {
        if (extensions.fragDepth) {
            prefix.push('#extension GL_EXT_frag_depth : enable')
            prefix.push('#define enabledFragDepth')
        } else {
            throw new Error(`requested 'GL_EXT_frag_depth' extension is unavailable`)
        }
    }
    if (shaderExtensions.drawBuffers) {
        if (extensions.drawBuffers) {
            prefix.push('#extension GL_EXT_draw_buffers : require')
            prefix.push('#define requiredDrawBuffers')
        } else {
            throw new Error(`requested 'GL_EXT_draw_buffers' extension is unavailable`)
        }
    }
    if (shaderExtensions.shaderTextureLod) {
        if (extensions.shaderTextureLod) {
            prefix.push('#extension GL_EXT_shader_texture_lod : enable')
            prefix.push('#define enabledShaderTextureLod')
        } else {
            throw new Error(`requested 'GL_EXT_shader_texture_lod' extension is unavailable`)
        }
    }
    return prefix.join('\n') + '\n'
}

const glsl300VertPrefix = `#version 300 es
#define attribute in
#define varying out
#define texture2D texture
`

const glsl300FragPrefix = `#version 300 es
layout(location = 0) out highp vec4 out_FragData0;
layout(location = 1) out highp vec4 out_FragData1;
layout(location = 2) out highp vec4 out_FragData2;
layout(location = 3) out highp vec4 out_FragData3;
layout(location = 4) out highp vec4 out_FragData4;
layout(location = 5) out highp vec4 out_FragData5;
layout(location = 6) out highp vec4 out_FragData6;
layout(location = 7) out highp vec4 out_FragData7;

#define varying in
#define texture2D texture
#define texture2DLodEXT textureLod

#define gl_FragColor out_FragData0
#define gl_FragDepthEXT gl_FragDepth

#define enabledStandardDerivatives
#define enabledFragDepth
#define requiredDrawBuffers
`

function transformGlsl300Frag(frag: string) {
    return frag.replace(/gl_FragData\[([0-7])\]/g, 'out_FragData$1')
}

export function addShaderDefines(gl: GLRenderingContext, extensions: WebGLExtensions, defines: ShaderDefines, shaders: ShaderCode): ShaderCode {
    const webgl2 = isWebGL2(gl)
    const header = getDefinesCode(defines)
    const vertPrefix = webgl2 ? glsl300VertPrefix : ''
    const fragPrefix = webgl2 ? glsl300FragPrefix : getGlsl100FragPrefix(extensions, shaders.extensions)
    const frag = webgl2 ? transformGlsl300Frag(shaders.frag) : shaders.frag
    return {
        id: shaderCodeId(),
        vert: `${vertPrefix}${header}${shaders.vert}`,
        frag: `${fragPrefix}${header}${frag}`,
        extensions: shaders.extensions
    }
}