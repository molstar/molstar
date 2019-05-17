/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { idFactory } from 'mol-util/id-factory';
import { WebGLExtensions } from './webgl/extensions';
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

import apply_fog from 'mol-gl/shader/chunks/apply-fog.glsl'
import apply_light_color from 'mol-gl/shader/chunks/apply-light-color.glsl'
import apply_marker_color from 'mol-gl/shader/chunks/apply-marker-color.glsl'
import assign_color_varying from 'mol-gl/shader/chunks/assign-color-varying.glsl'
import assign_group from 'mol-gl/shader/chunks/assign-group.glsl'
import assign_marker_varying from 'mol-gl/shader/chunks/assign-marker-varying.glsl'
import assign_material_color from 'mol-gl/shader/chunks/assign-material-color.glsl'
import assign_normal from 'mol-gl/shader/chunks/assign-normal.glsl'
import assign_position from 'mol-gl/shader/chunks/assign-position.glsl'
import assign_size from 'mol-gl/shader/chunks/assign-size.glsl'
import color_frag_params from 'mol-gl/shader/chunks/color-frag-params.glsl'
import color_vert_params from 'mol-gl/shader/chunks/color-vert-params.glsl'
import common_frag_params from 'mol-gl/shader/chunks/common-frag-params.glsl'
import common_vert_params from 'mol-gl/shader/chunks/common-vert-params.glsl'
import common from 'mol-gl/shader/chunks/common.glsl'
import light_frag_params from 'mol-gl/shader/chunks/light-frag-params.glsl'
import matrix_scale from 'mol-gl/shader/chunks/matrix-scale.glsl'
import normal_frag_params from 'mol-gl/shader/chunks/normal-frag-params.glsl'
import read_from_texture from 'mol-gl/shader/chunks/read-from-texture.glsl'
import size_vert_params from 'mol-gl/shader/chunks/size-vert-params.glsl'
import texture3d_from_2d_linear from 'mol-gl/shader/chunks/texture3d-from-2d-nearest.glsl'
import texture3d_from_2d_nearest from 'mol-gl/shader/chunks/texture3d-from-2d-nearest.glsl'

const ShaderChunks: { [k: string]: string } = {
    apply_fog,
    apply_light_color,
    apply_marker_color,
    assign_color_varying,
    assign_group,
    assign_marker_varying,
    assign_material_color,
    assign_normal,
    assign_position,
    assign_size,
    color_frag_params,
    color_vert_params,
    common_frag_params,
    common_vert_params,
    common,
    light_frag_params,
    matrix_scale,
    normal_frag_params,
    read_from_texture,
    size_vert_params,
    texture3d_from_2d_linear,
    texture3d_from_2d_nearest
}

const reInclude = /^(?!\/\/)\s*#include\s+(\S+)/gmi
const reSingleLineComment = /[ \t]*\/\/.*\n/g
const reMultiLineComment = /[ \t]*\/\*[\s\S]*?\*\//g
const reMultipleLinebreaks = /\n{2,}/g

function addIncludes(text: string) {
    return text
        .replace(reInclude, (_, p1) => {
            const chunk = ShaderChunks[p1]
            if (!chunk) throw new Error(`empty chunk, '${p1}'`)
            return chunk
        })
        .trim()
        .replace(reSingleLineComment, '\n')
        .replace(reMultiLineComment, '\n')
        .replace(reMultipleLinebreaks, '\n')
}

export function ShaderCode(vert: string, frag: string, extensions: ShaderExtensions = {}): ShaderCode {
    return { id: shaderCodeId(), vert: addIncludes(vert), frag: addIncludes(frag), extensions }
}

import points_vert from 'mol-gl/shader/points.vert'
import points_frag from 'mol-gl/shader/points.frag'
export const PointsShaderCode = ShaderCode(points_vert, points_frag)

import spheres_vert from 'mol-gl/shader/spheres.vert'
import spheres_frag from 'mol-gl/shader/spheres.frag'
export const SpheresShaderCode = ShaderCode(spheres_vert, spheres_frag, { fragDepth: true })

import text_vert from 'mol-gl/shader/text.vert'
import text_frag from 'mol-gl/shader/text.frag'
export const TextShaderCode = ShaderCode(text_vert, text_frag, { standardDerivatives: true })

import lines_vert from 'mol-gl/shader/lines.vert'
import lines_frag from 'mol-gl/shader/lines.frag'
export const LinesShaderCode = ShaderCode(lines_vert, lines_frag)

import mesh_vert from 'mol-gl/shader/mesh.vert'
import mesh_frag from 'mol-gl/shader/mesh.frag'
export const MeshShaderCode = ShaderCode(mesh_vert, mesh_frag, { standardDerivatives: true })

import direct_volume_vert from 'mol-gl/shader/direct-volume.vert'
import direct_volume_frag from 'mol-gl/shader/direct-volume.frag'
export const DirectVolumeShaderCode = ShaderCode(direct_volume_vert, direct_volume_frag, { fragDepth: true })

//

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