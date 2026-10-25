/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { assertUnreachable } from '@molstar/core/util/type-helpers';
import type { DefineLookup } from '../renderable/schema.js';
import type { IgnoreDefine, ShaderCode, ShaderDefines, ShaderExtensions } from '../shader-code.js';
import { type PreprocessorDefines, preprocessShader } from '../shader-preprocessor.js';
import { isWebGL2, type GLRenderingContext } from './compat.js';
import type { WebGLParameters } from './context.js';
import type { WebGLExtensions } from './extensions.js';
import { type FragOutTypes, getGlslProgram } from './programs.js';
import type { ShaderDeclarations } from './shader-declarations.js';

function replaceCounts(str: string, defines: ShaderDefines) {
  if (defines.dLightCount) str = str.replace(/dLightCount/g, `${defines.dLightCount.ref.value}`);
  return str;
}

/** Defines of `#define` lines (prefix and header), as the GLSL compiler sees them before the source */
function getPreprocessorDefines(code: string, version: '100' | '300'): PreprocessorDefines {
  const defines = new Map<string, string>([['__VERSION__', version]]);
  for (const m of code.matchAll(/^[ \t]*#define[ \t]+(\w+)(\()?[ \t]*(.*)$/gm)) {
    if (m[2]) throw new Error(`function-like macro '${m[1]}' in shader prefix or header`);
    defines.set(m[1], m[3].trim());
  }
  return defines;
}

/** Resolves the conditionals of `str` in TS, the active `#define`s are kept for the GLSL compiler */
function preprocess(str: string, shaderDefines: ShaderDefines, defines: PreprocessorDefines) {
  return preprocessShader(replaceCounts(str, shaderDefines), defines);
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
  // no `gl_VertexID` in GLSL ES 1.0
  const prefix: string[] = ['#define emulateVertexId'];
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
  if (extensions.clipControl) {
    // clip depth in [0, 1] (`EXT_clip_control`), see `WebGLContext.depthRangeConversion`
    prefix.push('#define clipZeroToOne');
  }
  return prefix.join('\n') + '\n';
}

function getGlsl100FragPrefix(extensions: WebGLExtensions, shaderExtensions: ShaderExtensions) {
  const prefix: string[] = ['#extension GL_OES_standard_derivatives : enable'];
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
  if (extensions.clipControl) {
    // clip depth in [0, 1] (`EXT_clip_control`), see `WebGLContext.depthRangeConversion`
    prefix.push('#define clipZeroToOne');
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
  const prefix = ['#version 300 es'];
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
  if (shaderExtensions.multiview2) {
    if (extensions.multiview2) {
      prefix.push('#extension GL_OVR_multiview2 : require');
      prefix.push('#define enabledMultiview2');
    } else if (shaderExtensions.multiview2 === 'required') {
      throw new Error(`required 'GL_OVR_multiview2' extension not available`);
    }
  }
  if (!extensions.noNonInstancedActiveAttribs) {
    // `gl_VertexID` is unusable, emulated with the `aVertex` attribute
    prefix.push('#define emulateVertexId');
  }
  if (extensions.clipControl) {
    // clip depth in [0, 1] (`EXT_clip_control`), see `WebGLContext.depthRangeConversion`
    prefix.push('#define clipZeroToOne');
  }
  prefix.push(glsl300VertPrefixCommon);
  return prefix.join('\n') + '\n';
}

/** prefix defines that are the same for all programs of a context, `undefined` for all other defines */
export function getContextDefineLookup(extensions: WebGLExtensions): DefineLookup {
  return (name) => {
    if (name === 'emulateVertexId') return !extensions.noNonInstancedActiveAttribs;
    if (name === 'clipZeroToOne') return !!extensions.clipControl;
    return undefined;
  };
}

function getGlsl300FragPrefix(
  extensions: WebGLExtensions,
  parameters: WebGLParameters,
  shaderExtensions: ShaderExtensions,
  outTypes: FragOutTypes,
) {
  const prefix = ['#version 300 es'];
  // extension directives come before any declaration
  if (shaderExtensions.conservativeDepth) {
    if (extensions.conservativeDepth) {
      prefix.push('#extension GL_EXT_conservative_depth : enable');
      prefix.push('#define enabledConservativeDepth');
    } else if (shaderExtensions.conservativeDepth === 'required') {
      throw new Error(`required 'GL_EXT_conservative_depth' extension not available`);
    }
  }
  prefix.push(`layout(location = 0) out highp ${outTypes[0] || 'vec4'} out_FragData0;`);
  if (shaderExtensions.fragDepth) {
    if (extensions.fragDepth) {
      prefix.push('#define enabledFragDepth');
    }
  }
  if (shaderExtensions.drawBuffers) {
    if (extensions.drawBuffers) {
      prefix.push('#define requiredDrawBuffers');
      for (let i = 1, il = parameters.maxDrawBuffers; i < il; ++i) {
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
  if (extensions.clipControl) {
    // clip depth in [0, 1] (`EXT_clip_control`), see `WebGLContext.depthRangeConversion`
    prefix.push('#define clipZeroToOne');
  }
  prefix.push(glsl300FragPrefixCommon);
  return prefix.join('\n') + '\n';
}

function transformGlsl300Frag(frag: string) {
  return frag.replace(/gl_FragData\[([0-9]+)\]/g, 'out_FragData$1');
}

export type ShaderUniformDeclarations = {
  /** GLSL uniform declarations generated from the schema, injected into both shader stages */
  code: string;
  /** GLSL (element) types of the declared uniforms, declarations in the shader code are removed */
  types: Map<string, string>;
};

const reUniformDeclaration = /^[ \t]*uniform\s+((?:highp\s+|mediump\s+|lowp\s+)?\w+)\s+([^;{}]+);/gm;

/** Remove declarations of uniforms that are generated, see `getShaderUniformDeclarations`. */
export function stripUniformDeclarations(source: string, types: ReadonlyMap<string, string>) {
  return source.replace(reUniformDeclaration, (match, type: string, names: string) => {
    if (type.includes('sampler')) return match;
    const baseType = type.replace(/^(highp|mediump|lowp)\s+/, '');
    const kept: string[] = [];
    for (const n of names.split(',')) {
      const name = n.trim();
      const generatedType = types.get(name.replace(/\[[^\]]*\]$/, '').trim());
      if (generatedType === undefined) {
        kept.push(name);
      } else if (generatedType !== baseType) {
        throw new Error(
          `uniform '${name}' is declared as '${baseType}' in the shader but as '${generatedType}' in the schema`,
        );
      }
    }
    if (kept.length === 0) return '';
    return `uniform ${type} ${kept.join(', ')};`;
  });
}

const reSamplerDeclaration = /^[ \t]*uniform\s+(?:(?:highp|mediump|lowp)\s+)?(\w*sampler\w+)\s+(\w+)\s*;/gm;

/** Remove declarations of samplers that are generated, see `getShaderDeclarations`. */
export function stripSamplerDeclarations(source: string, types: ReadonlyMap<string, readonly string[]>) {
  return source.replace(reSamplerDeclaration, (match, type: string, name: string) => {
    const expected = types.get(name);
    if (expected === undefined) return match;
    if (!expected.includes(type)) {
      throw new Error(`sampler '${name}' is declared as '${type}' in the shader, expected one of '${expected}'`);
    }
    return '';
  });
}

const reAttributeDeclaration = /^[ \t]*attribute\s+(?:(?:highp|mediump|lowp)\s+)?(\w+)\s+(\w+)\s*;/gm;

/** Remove declarations of attributes that are generated, see `getShaderDeclarations`. */
export function stripAttributeDeclarations(source: string, types: ReadonlyMap<string, string>) {
  return source.replace(reAttributeDeclaration, (match, type: string, name: string) => {
    const expected = types.get(name);
    if (expected === undefined) return match;
    if (expected !== type) {
      throw new Error(`attribute '${name}' is declared as '${type}' in the shader but as '${expected}' in the schema`);
    }
    return '';
  });
}

/** GLSL sources of a program, with prefixes, defines and generated declarations */
export interface GlslSources {
  readonly vert: string;
  readonly frag: string;
}

/** The GLSL of `shaderCode` (see `getGlslProgram`) with prefixes for the context, defines and declarations */
export function addShaderDefines(
  gl: GLRenderingContext,
  extensions: WebGLExtensions,
  parameters: WebGLParameters,
  defines: ShaderDefines,
  shaderCode: ShaderCode,
  declarations?: ShaderDeclarations,
): GlslSources {
  const source = getGlslProgram(shaderCode.name);
  const vertHeader = getDefinesCode(defines, shaderCode.ignoreDefine);
  const fragHeader = getDefinesCode(defines, shaderCode.ignoreDefine);
  const vertPrefix = isWebGL2(gl)
    ? getGlsl300VertPrefix(extensions, shaderCode.extensions)
    : getGlsl100VertPrefix(extensions, shaderCode.extensions);
  const fragPrefix = isWebGL2(gl)
    ? getGlsl300FragPrefix(extensions, parameters, shaderCode.extensions, source.outTypes)
    : getGlsl100FragPrefix(extensions, shaderCode.extensions);
  const version = isWebGL2(gl) ? '300' : '100';
  const vertDefines = getPreprocessorDefines(`${vertPrefix}${vertHeader}`, version);
  const fragDefines = getPreprocessorDefines(`${fragPrefix}${fragHeader}`, version);
  let vert = preprocess(source.vert, defines, vertDefines);
  let frag = preprocess(isWebGL2(gl) ? transformGlsl300Frag(source.frag) : source.frag, defines, fragDefines);
  if (declarations) {
    const { uniformTypes, samplerTypes, attributeTypes } = declarations;
    // counts may be ignored defines (e.g. `dLightCount`), see `replaceCounts`
    const vertDeclarations = preprocess(declarations.vert, defines, vertDefines);
    const fragDeclarations = preprocess(declarations.frag, defines, fragDefines);
    vert = stripUniformDeclarations(stripSamplerDeclarations(vert, samplerTypes), uniformTypes);
    vert = `${vertDeclarations}\n${stripAttributeDeclarations(vert, attributeTypes)}`;
    frag = `${fragDeclarations}\n${stripUniformDeclarations(stripSamplerDeclarations(frag, samplerTypes), uniformTypes)}`;
  }
  return {
    vert: `${vertPrefix}${vertHeader}${vert}`,
    frag: `${fragPrefix}${fragHeader}${frag}`,
  };
}
