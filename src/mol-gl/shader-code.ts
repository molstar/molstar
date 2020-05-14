/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../mol-util';
import { idFactory } from '../mol-util/id-factory';
import { WebGLExtensions } from './webgl/extensions';
import { isWebGL2, GLRenderingContext } from './webgl/compat';

export type DefineKind = 'boolean' | 'string' | 'number'
export type DefineType = boolean | string
export type DefineValues = { [k: string]: ValueCell<DefineType> }

const shaderCodeId = idFactory();

export interface ShaderExtensions {
    readonly standardDerivatives?: boolean
    readonly fragDepth?: boolean
    readonly drawBuffers?: boolean
    readonly shaderTextureLod?: boolean
}

export interface ShaderCode {
    readonly id: number
    readonly name: string
    readonly vert: string
    readonly frag: string
    readonly extensions: ShaderExtensions
}

import apply_fog from './shader/chunks/apply-fog.glsl';
import apply_interior_color from './shader/chunks/apply-interior-color.glsl';
import apply_light_color from './shader/chunks/apply-light-color.glsl';
import apply_marker_color from './shader/chunks/apply-marker-color.glsl';
import assign_clipping_varying from './shader/chunks/assign-clipping-varying.glsl';
import assign_color_varying from './shader/chunks/assign-color-varying.glsl';
import assign_group from './shader/chunks/assign-group.glsl';
import assign_marker_varying from './shader/chunks/assign-marker-varying.glsl';
import assign_material_color from './shader/chunks/assign-material-color.glsl';
import assign_position from './shader/chunks/assign-position.glsl';
import assign_size from './shader/chunks/assign-size.glsl';
import check_picking_alpha from './shader/chunks/check-picking-alpha.glsl';
import clip_instance from './shader/chunks/clip-instance.glsl';
import clip_pixel from './shader/chunks/clip-pixel.glsl';
import color_frag_params from './shader/chunks/color-frag-params.glsl';
import color_vert_params from './shader/chunks/color-vert-params.glsl';
import common_clip from './shader/chunks/common-clip.glsl';
import common_frag_params from './shader/chunks/common-frag-params.glsl';
import common_vert_params from './shader/chunks/common-vert-params.glsl';
import common from './shader/chunks/common.glsl';
import light_frag_params from './shader/chunks/light-frag-params.glsl';
import matrix_scale from './shader/chunks/matrix-scale.glsl';
import normal_frag_params from './shader/chunks/normal-frag-params.glsl';
import read_from_texture from './shader/chunks/read-from-texture.glsl';
import size_vert_params from './shader/chunks/size-vert-params.glsl';
import texture3d_from_2d_linear from './shader/chunks/texture3d-from-2d-nearest.glsl';
import texture3d_from_2d_nearest from './shader/chunks/texture3d-from-2d-nearest.glsl';

const ShaderChunks: { [k: string]: string } = {
    apply_fog,
    apply_interior_color,
    apply_light_color,
    apply_marker_color,
    assign_clipping_varying,
    assign_color_varying,
    assign_group,
    assign_marker_varying,
    assign_material_color,
    assign_position,
    assign_size,
    check_picking_alpha,
    clip_instance,
    clip_pixel,
    color_frag_params,
    color_vert_params,
    common_clip,
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
};

const reInclude = /^(?!\/\/)\s*#include\s+(\S+)/gmi;
const reSingleLineComment = /[ \t]*\/\/.*\n/g;
const reMultiLineComment = /[ \t]*\/\*[\s\S]*?\*\//g;
const reMultipleLinebreaks = /\n{2,}/g;

function addIncludes(text: string) {
    return text
        .replace(reInclude, (_, p1) => {
            const chunk = ShaderChunks[p1];
            if (!chunk) throw new Error(`empty chunk, '${p1}'`);
            return chunk;
        })
        .trim()
        .replace(reSingleLineComment, '\n')
        .replace(reMultiLineComment, '\n')
        .replace(reMultipleLinebreaks, '\n');
}

export function ShaderCode(name: string, vert: string, frag: string, extensions: ShaderExtensions = {}): ShaderCode {
    return { id: shaderCodeId(), name, vert: addIncludes(vert), frag: addIncludes(frag), extensions };
}

import points_vert from './shader/points.vert';
import points_frag from './shader/points.frag';
export const PointsShaderCode = ShaderCode('points', points_vert, points_frag);

import spheres_vert from './shader/spheres.vert';
import spheres_frag from './shader/spheres.frag';
export const SpheresShaderCode = ShaderCode('spheres', spheres_vert, spheres_frag, { fragDepth: true });

import text_vert from './shader/text.vert';
import text_frag from './shader/text.frag';
export const TextShaderCode = ShaderCode('text', text_vert, text_frag, { standardDerivatives: true });

import lines_vert from './shader/lines.vert';
import lines_frag from './shader/lines.frag';
export const LinesShaderCode = ShaderCode('lines', lines_vert, lines_frag);

import mesh_vert from './shader/mesh.vert';
import mesh_frag from './shader/mesh.frag';
export const MeshShaderCode = ShaderCode('mesh', mesh_vert, mesh_frag, { standardDerivatives: true });

import direct_volume_vert from './shader/direct-volume.vert';
import direct_volume_frag from './shader/direct-volume.frag';
export const DirectVolumeShaderCode = ShaderCode('direct-volume', direct_volume_vert, direct_volume_frag, { fragDepth: true });

import image_vert from './shader/image.vert';
import image_frag from './shader/image.frag';
export const ImageShaderCode = ShaderCode('image', image_vert, image_frag, { fragDepth: true });

//

export type ShaderDefines = {
    [k: string]: ValueCell<DefineType>
}

function getDefinesCode (defines: ShaderDefines) {
    if (defines === undefined) return '';
    const lines = [];
    for (const name in defines) {
        const define = defines[name];
        const v = define.ref.value;
        if (v !== undefined) {
            if (typeof v === 'string') {
                lines.push(`#define ${name}_${v}`);
            } else if (typeof v === 'number') {
                lines.push(`#define ${name} ${v}`);
            } else if (typeof v === 'boolean') {
                if (v) lines.push(`#define ${name}`);
            } else {
                throw new Error('unknown define type');
            }
        }
    }
    return lines.join('\n') + '\n';
}

function getGlsl100FragPrefix(extensions: WebGLExtensions, shaderExtensions: ShaderExtensions) {
    const prefix: string[] = [];
    if (shaderExtensions.standardDerivatives) {
        prefix.push('#extension GL_OES_standard_derivatives : enable');
        prefix.push('#define enabledStandardDerivatives');
    }
    if (shaderExtensions.fragDepth) {
        if (extensions.fragDepth) {
            prefix.push('#extension GL_EXT_frag_depth : enable');
            prefix.push('#define enabledFragDepth');
        } else {
            throw new Error(`requested 'GL_EXT_frag_depth' extension is unavailable`);
        }
    }
    if (shaderExtensions.drawBuffers) {
        if (extensions.drawBuffers) {
            prefix.push('#extension GL_EXT_draw_buffers : require');
            prefix.push('#define requiredDrawBuffers');
        } else {
            throw new Error(`requested 'GL_EXT_draw_buffers' extension is unavailable`);
        }
    }
    if (shaderExtensions.shaderTextureLod) {
        if (extensions.shaderTextureLod) {
            prefix.push('#extension GL_EXT_shader_texture_lod : enable');
            prefix.push('#define enabledShaderTextureLod');
        } else {
            throw new Error(`requested 'GL_EXT_shader_texture_lod' extension is unavailable`);
        }
    }
    return prefix.join('\n') + '\n';
}

const glsl300VertPrefix = `#version 300 es
#define attribute in
#define varying out
#define texture2D texture
`;

const glsl300FragPrefixCommon = `
#define varying in
#define texture2D texture
#define texture2DLodEXT textureLod

#define gl_FragColor out_FragData0
#define gl_FragDepthEXT gl_FragDepth

#define requiredDrawBuffers
`;

function getGlsl300FragPrefix(gl: WebGL2RenderingContext, extensions: WebGLExtensions, shaderExtensions: ShaderExtensions) {
    const prefix = [ '#version 300 es' ];
    if (shaderExtensions.standardDerivatives) {
        prefix.push('#define enabledStandardDerivatives');
    }
    if (shaderExtensions.fragDepth) {
        prefix.push('#define enabledFragDepth');
    }
    if (extensions.drawBuffers) {
        const maxDrawBuffers = gl.getParameter(gl.MAX_DRAW_BUFFERS) as number;
        for (let i = 0, il = maxDrawBuffers; i < il; ++i) {
            prefix.push(`layout(location = ${i}) out highp vec4 out_FragData${i};`);
        }
    }
    prefix.push(glsl300FragPrefixCommon);
    return prefix.join('\n') + '\n';
}

function transformGlsl300Frag(frag: string) {
    return frag.replace(/gl_FragData\[([0-9]+)\]/g, 'out_FragData$1');
}

export function addShaderDefines(gl: GLRenderingContext, extensions: WebGLExtensions, defines: ShaderDefines, shaders: ShaderCode): ShaderCode {
    const header = getDefinesCode(defines);
    const vertPrefix = isWebGL2(gl) ? glsl300VertPrefix : '';
    const fragPrefix = isWebGL2(gl)
        ? getGlsl300FragPrefix(gl, extensions, shaders.extensions)
        : getGlsl100FragPrefix(extensions, shaders.extensions);
    const frag = isWebGL2(gl) ? transformGlsl300Frag(shaders.frag) : shaders.frag;
    return {
        id: shaderCodeId(),
        name: shaders.name,
        vert: `${vertPrefix}${header}${shaders.vert}`,
        frag: `${fragPrefix}${header}${frag}`,
        extensions: shaders.extensions
    };
}