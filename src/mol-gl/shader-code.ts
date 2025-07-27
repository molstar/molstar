/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ValueCell } from '../mol-util';
import { idFactory } from '../mol-util/id-factory';
import { WebGLExtensions } from './webgl/extensions';
import { isWebGL2, GLRenderingContext } from './webgl/compat';
import { assertUnreachable } from '../mol-util/type-helpers';

export type DefineKind = 'boolean' | 'string' | 'number'
export type DefineType = boolean | string
export type DefineValues = { [k: string]: ValueCell<DefineType> }

const shaderCodeId = idFactory();

type ShaderExtensionsValue = 'required' | 'optional'
export interface ShaderExtensions {
    readonly fragDepth?: ShaderExtensionsValue
    readonly drawBuffers?: ShaderExtensionsValue
    readonly shaderTextureLod?: ShaderExtensionsValue
    /** Needed to enable the `gl_DrawID` built-in */
    readonly multiDraw?: ShaderExtensionsValue
    readonly clipCullDistance?: ShaderExtensionsValue
    readonly conservativeDepth?: ShaderExtensionsValue
    /** Needed to enable the `gl_ViewID_OVR` built-in */
    readonly multiview2?: ShaderExtensionsValue
}

type FragOutTypes = { [k in number]: 'vec4' | 'ivec4' }
type IgnoreDefine = (name: string, variant: string, defines: ShaderDefines) => boolean

export interface ShaderCode {
    readonly id: number
    readonly name: string
    readonly vert: string
    readonly frag: string
    readonly extensions: ShaderExtensions
    /** Fragment shader output type only applicable for webgl2 */
    readonly outTypes: FragOutTypes
    readonly ignoreDefine?: IgnoreDefine
}

import { apply_fog } from './shader/chunks/apply-fog.glsl';
import { apply_interior_color } from './shader/chunks/apply-interior-color.glsl';
import { apply_light_color } from './shader/chunks/apply-light-color.glsl';
import { apply_marker_color } from './shader/chunks/apply-marker-color.glsl';
import { assign_clipping_varying } from './shader/chunks/assign-clipping-varying.glsl';
import { assign_color_varying } from './shader/chunks/assign-color-varying.glsl';
import { assign_group } from './shader/chunks/assign-group.glsl';
import { assign_marker_varying } from './shader/chunks/assign-marker-varying.glsl';
import { assign_material_color } from './shader/chunks/assign-material-color.glsl';
import { assign_position } from './shader/chunks/assign-position.glsl';
import { assign_size } from './shader/chunks/assign-size.glsl';
import { check_picking_alpha } from './shader/chunks/check-picking-alpha.glsl';
import { check_transparency } from './shader/chunks/check-transparency.glsl';
import { clip_instance } from './shader/chunks/clip-instance.glsl';
import { clip_pixel } from './shader/chunks/clip-pixel.glsl';
import { color_frag_params } from './shader/chunks/color-frag-params.glsl';
import { color_vert_params } from './shader/chunks/color-vert-params.glsl';
import { common_clip } from './shader/chunks/common-clip.glsl';
import { common_frag_params } from './shader/chunks/common-frag-params.glsl';
import { common_vert_params } from './shader/chunks/common-vert-params.glsl';
import { common } from './shader/chunks/common.glsl';
import { fade_lod } from './shader/chunks/fade-lod.glsl';
import { float_to_rgba } from './shader/chunks/float-to-rgba.glsl';
import { light_frag_params } from './shader/chunks/light-frag-params.glsl';
import { normal_frag_params } from './shader/chunks/normal-frag-params.glsl';
import { read_from_texture } from './shader/chunks/read-from-texture.glsl';
import { rgba_to_float } from './shader/chunks/rgba-to-float.glsl';
import { size_vert_params } from './shader/chunks/size-vert-params.glsl';
import { texture3d_from_1d_trilinear } from './shader/chunks/texture3d-from-1d-trilinear.glsl';
import { texture3d_from_2d_linear } from './shader/chunks/texture3d-from-2d-linear.glsl';
import { texture3d_from_2d_nearest } from './shader/chunks/texture3d-from-2d-nearest.glsl';
import { wboit_write } from './shader/chunks/wboit-write.glsl';
import { dpoit_write } from './shader/chunks/dpoit-write.glsl';

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
    check_transparency,
    clip_instance,
    clip_pixel,
    color_frag_params,
    color_vert_params,
    common_clip,
    common_frag_params,
    common_vert_params,
    common,
    fade_lod,
    float_to_rgba,
    light_frag_params,
    normal_frag_params,
    read_from_texture,
    rgba_to_float,
    size_vert_params,
    texture3d_from_1d_trilinear,
    texture3d_from_2d_linear,
    texture3d_from_2d_nearest,
    wboit_write,
    dpoit_write
};

const reInclude = /^(?!\/\/)\s*#include\s+(\S+)/gm;
const reUnrollLoop = /#pragma unroll_loop_start\s+for\s*\(\s*int\s+i\s*=\s*(\d+)\s*;\s*i\s*<\s*(\d+)\s*;\s*\+\+i\s*\s*\)\s*{([\s\S]+?)}\s+#pragma unroll_loop_end/g;
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

function unrollLoops(str: string) {
    return str.replace(reUnrollLoop, loopReplacer);
}

function loopReplacer(_match: string, start: string, end: string, snippet: string) {
    let out = '';
    for (let i = parseInt(start); i < parseInt(end); ++i) {
        out += snippet
            .replace(/\[\s*i\s*\]/g, `[${i}]`)
            .replace(/UNROLLED_LOOP_INDEX/g, `${i}`);
    }
    return out;
}

function replaceCounts(str: string, defines: ShaderDefines) {
    if (defines.dLightCount) str = str.replace(/dLightCount/g, `${defines.dLightCount.ref.value}`);
    if (defines.dClipObjectCount) str = str.replace(/dClipObjectCount/g, `${defines.dClipObjectCount.ref.value}`);
    return str;
}

function preprocess(str: string, defines: ShaderDefines) {
    return unrollLoops(replaceCounts(str, defines));
}

export function ShaderCode(name: string, vert: string, frag: string, extensions: ShaderExtensions = {}, outTypes: FragOutTypes = {}, ignoreDefine?: IgnoreDefine): ShaderCode {
    return { id: shaderCodeId(), name, vert: addIncludes(vert), frag: addIncludes(frag), extensions, outTypes, ignoreDefine };
}

// Note: `drawBuffers` need to be 'optional' for wboit

function ignoreDefine(name: string, variant: string, defines: ShaderDefines): boolean {
    if (variant.startsWith('color') || variant === 'tracing') {
        if (name === 'dLightCount') {
            return !!defines.dIgnoreLight?.ref.value;
        }
    } else {
        const ignore = [
            'dColorType', 'dUsePalette',
            'dOverpaintType', 'dOverpaint',
            'dSubstanceType', 'dSubstance',
            'dColorMarker', 'dCelShaded',
            'dLightCount',
        ];
        if (variant !== 'depth' && !variant.startsWith('pick')) {
            ignore.push('dXrayShaded');
        }
        if (variant !== 'emissive') {
            ignore.push('dEmissiveType', 'dEmissive');
        }
        return ignore.includes(name);
    }
    return false;
};

function ignoreDefineUnlit(name: string, variant: string, defines: ShaderDefines): boolean {
    if (name === 'dLightCount') return true;
    return ignoreDefine(name, variant, defines);
};

import { points_vert } from './shader/points.vert';
import { points_frag } from './shader/points.frag';
export const PointsShaderCode = ShaderCode('points', points_vert, points_frag, { drawBuffers: 'optional' }, {}, ignoreDefineUnlit);

import { spheres_vert } from './shader/spheres.vert';
import { spheres_frag } from './shader/spheres.frag';
export const SpheresShaderCode = ShaderCode('spheres', spheres_vert, spheres_frag, { fragDepth: 'required', drawBuffers: 'optional' }, {}, ignoreDefine);

import { cylinders_vert } from './shader/cylinders.vert';
import { cylinders_frag } from './shader/cylinders.frag';
export const CylindersShaderCode = ShaderCode('cylinders', cylinders_vert, cylinders_frag, { fragDepth: 'required', drawBuffers: 'optional' }, {}, ignoreDefine);

import { text_vert } from './shader/text.vert';
import { text_frag } from './shader/text.frag';
export const TextShaderCode = ShaderCode('text', text_vert, text_frag, { drawBuffers: 'optional' }, {}, ignoreDefineUnlit);

import { lines_vert } from './shader/lines.vert';
import { lines_frag } from './shader/lines.frag';
export const LinesShaderCode = ShaderCode('lines', lines_vert, lines_frag, { drawBuffers: 'optional' }, {}, ignoreDefineUnlit);

import { mesh_vert } from './shader/mesh.vert';
import { mesh_frag } from './shader/mesh.frag';
export const MeshShaderCode = ShaderCode('mesh', mesh_vert, mesh_frag, { drawBuffers: 'optional' }, {}, ignoreDefine);

import { directVolume_vert } from './shader/direct-volume.vert';
import { directVolume_frag } from './shader/direct-volume.frag';
export const DirectVolumeShaderCode = ShaderCode('direct-volume', directVolume_vert, directVolume_frag, { fragDepth: 'optional', drawBuffers: 'optional' }, {}, ignoreDefine);

import { image_vert } from './shader/image.vert';
import { image_frag } from './shader/image.frag';
export const ImageShaderCode = ShaderCode('image', image_vert, image_frag, { drawBuffers: 'optional' }, {}, ignoreDefineUnlit);

//

export type ShaderDefines = {
    [k: string]: ValueCell<DefineType>
}

function getDefinesCode(defines: ShaderDefines, ignore?: IgnoreDefine) {
    if (defines === undefined) return '';
    const variant = (defines.dRenderVariant?.ref.value || '') as string;

    const lines = [];
    for (const name in defines) {
        if (ignore?.(name, variant, defines)) continue;

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
                assertUnreachable(v);
            }
        }
    }
    return lines.join('\n') + '\n';
}

function getGlsl100VertPrefix(extensions: WebGLExtensions, shaderExtensions: ShaderExtensions) {
    const prefix: string[] = [];
    if (shaderExtensions.drawBuffers) {
        if (extensions.drawBuffers) {
            prefix.push('#define requiredDrawBuffers');
        } else if (shaderExtensions.drawBuffers === 'required') {
            throw new Error(`required 'GL_EXT_draw_buffers' extension not available`);
        }
    }
    if (shaderExtensions.multiDraw) {
        if (extensions.multiDraw) {
            prefix.push('#extension GL_ANGLE_multi_draw : require');
            prefix.push('#define enabledMultiDraw');
        } else if (shaderExtensions.multiDraw === 'required') {
            throw new Error(`required 'GL_ANGLE_multi_draw' extension not available`);
        }
    }
    return prefix.join('\n') + '\n';
}

function getGlsl100FragPrefix(extensions: WebGLExtensions, shaderExtensions: ShaderExtensions) {
    const prefix: string[] = [
        '#extension GL_OES_standard_derivatives : enable'
    ];
    if (shaderExtensions.fragDepth) {
        if (extensions.fragDepth) {
            prefix.push('#extension GL_EXT_frag_depth : enable');
            prefix.push('#define enabledFragDepth');
        } else if (shaderExtensions.fragDepth === 'required') {
            throw new Error(`required 'GL_EXT_frag_depth' extension not available`);
        }
    }
    if (shaderExtensions.drawBuffers) {
        if (extensions.drawBuffers) {
            prefix.push('#extension GL_EXT_draw_buffers : require');
            prefix.push('#define requiredDrawBuffers');
            prefix.push('#define gl_FragColor gl_FragData[0]');
        } else if (shaderExtensions.drawBuffers === 'required') {
            throw new Error(`required 'GL_EXT_draw_buffers' extension not available`);
        }
    }
    if (shaderExtensions.shaderTextureLod) {
        if (extensions.shaderTextureLod) {
            prefix.push('#extension GL_EXT_shader_texture_lod : enable');
            prefix.push('#define enabledShaderTextureLod');
        } else if (shaderExtensions.shaderTextureLod === 'required') {
            throw new Error(`required 'GL_EXT_shader_texture_lod' extension not available`);
        }
    }
    if (extensions.depthTexture) {
        prefix.push('#define depthTextureSupport');
    }
    return prefix.join('\n') + '\n';
}

const glsl300VertPrefixCommon = `
#define attribute in
#define varying out
#define texture2D texture
`;

const glsl300FragPrefixCommon = `
#define varying in
#define texture2D texture
#define textureCube texture
#define texture2DLodEXT textureLod
#define textureCubeLodEXT textureLod

#define gl_FragColor out_FragData0
#define gl_FragDepthEXT gl_FragDepth
`;

function getGlsl300VertPrefix(extensions: WebGLExtensions, shaderExtensions: ShaderExtensions) {
    const prefix = [
        '#version 300 es',
    ];
    if (shaderExtensions.drawBuffers) {
        if (extensions.drawBuffers) {
            prefix.push('#define requiredDrawBuffers');
        }
    }
    if (shaderExtensions.multiDraw) {
        if (extensions.multiDraw) {
            prefix.push('#extension GL_ANGLE_multi_draw : require');
            prefix.push('#define enabledMultiDraw');
        } else if (shaderExtensions.multiDraw === 'required') {
            throw new Error(`required 'GL_ANGLE_multi_draw' extension not available`);
        }
    }
    if (shaderExtensions.clipCullDistance) {
        if (extensions.clipCullDistance) {
            prefix.push('#extension GL_ANGLE_clip_cull_distance : enable');
            prefix.push('#define enabledClipCullDistance');
        } else if (shaderExtensions.clipCullDistance === 'required') {
            throw new Error(`required 'GL_ANGLE_clip_cull_distance' extension not available`);
        }
    }
    if (shaderExtensions.conservativeDepth) {
        if (extensions.conservativeDepth) {
            prefix.push('#extension GL_EXT_conservative_depth : enable');
            prefix.push('#define enabledConservativeDepth');
        } else if (shaderExtensions.conservativeDepth === 'required') {
            throw new Error(`required 'GL_EXT_conservative_depth' extension not available`);
        }
    }
    if (shaderExtensions.multiview2) {
        if (extensions.multiview2) {
            prefix.push('#extension GL_OVR_multiview2 : require');
            prefix.push('#define enabledMultiview2');
        } else if (shaderExtensions.multiview2 === 'required') {
            throw new Error(`required 'GL_OVR_multiview2' extension not available`);
        }
    }
    if (extensions.noNonInstancedActiveAttribs) {
        prefix.push('#define noNonInstancedActiveAttribs');
    }
    prefix.push(glsl300VertPrefixCommon);
    return prefix.join('\n') + '\n';
}

function getGlsl300FragPrefix(gl: WebGL2RenderingContext, extensions: WebGLExtensions, shaderExtensions: ShaderExtensions, outTypes: FragOutTypes) {
    const prefix = [
        '#version 300 es',
        `layout(location = 0) out highp ${outTypes[0] || 'vec4'} out_FragData0;`
    ];
    if (shaderExtensions.fragDepth) {
        if (extensions.fragDepth) {
            prefix.push('#define enabledFragDepth');
        }
    }
    if (shaderExtensions.drawBuffers) {
        if (extensions.drawBuffers) {
            prefix.push('#define requiredDrawBuffers');
            const maxDrawBuffers = gl.getParameter(gl.MAX_DRAW_BUFFERS) as number;
            for (let i = 1, il = maxDrawBuffers; i < il; ++i) {
                prefix.push(`layout(location = ${i}) out highp ${outTypes[i] || 'vec4'} out_FragData${i};`);
            }
        }
    }
    if (shaderExtensions.shaderTextureLod) {
        if (extensions.shaderTextureLod) {
            prefix.push('#define enabledShaderTextureLod');
        }
    }
    if (extensions.depthTexture) {
        prefix.push('#define depthTextureSupport');
    }
    prefix.push(glsl300FragPrefixCommon);
    return prefix.join('\n') + '\n';
}

function transformGlsl300Frag(frag: string) {
    return frag.replace(/gl_FragData\[([0-9]+)\]/g, 'out_FragData$1');
}

export function addShaderDefines(gl: GLRenderingContext, extensions: WebGLExtensions, defines: ShaderDefines, shaders: ShaderCode): ShaderCode {
    const vertHeader = getDefinesCode(defines, shaders.ignoreDefine);
    const fragHeader = getDefinesCode(defines, shaders.ignoreDefine);
    const vertPrefix = isWebGL2(gl)
        ? getGlsl300VertPrefix(extensions, shaders.extensions)
        : getGlsl100VertPrefix(extensions, shaders.extensions);
    const fragPrefix = isWebGL2(gl)
        ? getGlsl300FragPrefix(gl, extensions, shaders.extensions, shaders.outTypes)
        : getGlsl100FragPrefix(extensions, shaders.extensions);
    const frag = isWebGL2(gl) ? transformGlsl300Frag(shaders.frag) : shaders.frag;
    return {
        id: shaderCodeId(),
        name: shaders.name,
        vert: `${vertPrefix}${vertHeader}${preprocess(shaders.vert, defines)}`,
        frag: `${fragPrefix}${fragHeader}${preprocess(frag, defines)}`,
        extensions: shaders.extensions,
        outTypes: shaders.outTypes
    };
}