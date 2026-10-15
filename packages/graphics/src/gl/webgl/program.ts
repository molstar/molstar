/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ShaderCode, DefineValues } from '../shader-code.js';
import { addShaderDefines, type GlslSources } from './shader-code.js';
import type { WebGLState } from './state.js';
import type { WebGLExtensions } from './extensions.js';
import { getUniformSetters, getUniformType, type UniformSetters } from './uniform.js';
import { isArrayUniform, type UniformsList, type UniformType } from '../uniform.js';
import { getUniformBlockLayout, UniformBlockInfo, UniformBlocks, getUniformBlock } from '../uniform-layout.js';
import { getShaderDeclarations } from './shader-declarations.js';
import { type AttributeBuffers, bindAttributeBuffers, getAttribType } from './buffer.js';
import { type AttributeLayout, getAttributeLayout } from '../attribute-layout.js';
import { asGLTexture, type GLTexture, type TextureId } from './texture.js';
import { isNullTexture, type Textures } from '../texture.js';
import { idFactory } from '@molstar/core/util/id-factory';
import {
  type DefineLookup,
  getDefineLookup,
  getDefineValuesLookup,
  getTextureDimension,
  isUsed,
  type RenderableSchema,
  type TextureDimension,
} from '../renderable/schema.js';
import { isDebugMode } from '@molstar/core/util/debug';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import type { ShaderType, Shader } from './shader.js';
import type { WebGLParameters } from './context.js';
import type { GraphicsProgram, ProgramVariant } from '../render-item.js';

const getNextProgramId = idFactory();

export interface Program extends GraphicsProgram {
  use: () => void;
  setUniforms: (uniformValues: UniformsList) => void;
  uniform: (k: string, v: UniformType) => void;
  bindAttributes: (attribueBuffers: AttributeBuffers) => void;
  offsetAttributes: (attributeBuffers: AttributeBuffers, offset: number) => void;
  bindTextures: (textures: Textures, startingTargetUnit: number) => void;

  reset: () => void;
  destroy: () => void;
}

export type Programs = { [k: string]: Program };

type Locations = { [k: string]: number };

function getLocations(gl: GLRenderingContext, program: WebGLProgram, schema: RenderableSchema, useBlocks: boolean) {
  const locations: Locations = {};
  Object.keys(schema).forEach((k) => {
    const spec = schema[k];
    if (spec.type === 'uniform') {
      // uniforms in blocks have no location, they are backed by a buffer
      if (useBlocks && getUniformBlock(spec) !== undefined) return;
      let loc = gl.getUniformLocation(program, k);
      // headless-gl requires a '[0]' suffix for array uniforms (https://github.com/stackgl/headless-gl/issues/170)
      if (loc === null && isArrayUniform(spec.kind)) loc = gl.getUniformLocation(program, k + '[0]');
      // unused uniforms will result in a `null` location which is usually fine
      // if (loc === null) console.info(`Could not get uniform location for '${k}'`);
      locations[k] = loc as number;
    } else if (spec.type === 'texture') {
      const loc = gl.getUniformLocation(program, k);
      // unused uniforms will result in a `null` location which is usually fine
      // if (loc === null) console.info(`Could not get uniform location for '${k}'`);
      locations[k] = loc as number;
    }
  });
  return locations;
}

function checkActiveAttributes(
  gl: GLRenderingContext,
  program: WebGLProgram,
  schema: RenderableSchema,
  defines: DefineLookup,
) {
  const attribCount = gl.getProgramParameter(program, gl.ACTIVE_ATTRIBUTES);
  for (let i = 0; i < attribCount; ++i) {
    const info = gl.getActiveAttrib(program, i);
    if (info) {
      const { name, type } = info;
      if (name.startsWith('__activeAttribute')) {
        // name assigned by `gl.shim.ts`, ignore for checks
        continue;
      }
      if (name === 'gl_InstanceID') continue; // WebGL2 built-in
      if (name === 'gl_VertexID') continue; // WebGL2 built-in
      if (name === 'gl_DrawID') continue; // WEBGL_multi_draw built-in
      if (name === 'gl_ViewID_OVR') continue; // OVR_multiview2 built-in
      const spec = schema[name];
      if (spec === undefined) {
        throw new Error(`missing 'uniform' or 'texture' with name '${name}' in schema`);
      }
      if (spec.type !== 'attribute') {
        throw new Error(`'${name}' must be of type 'attribute' but is '${spec.type}'`);
      }
      const attribType = getAttribType(gl, spec.kind, spec.itemSize);
      if (attribType !== type) {
        throw new Error(`unexpected attribute type '${attribType}' for ${name}, expected '${type}'`);
      }
      if (!isUsed(spec, defines)) {
        throw new Error(`attribute '${name}' is used, but its spec 'usedIf' condition does not hold`);
      }
    }
  }
}

function checkActiveUniforms(
  gl: GLRenderingContext,
  program: WebGLProgram,
  schema: RenderableSchema,
  defines: DefineLookup,
) {
  const attribCount = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < attribCount; ++i) {
    const info = gl.getActiveUniform(program, i);
    if (info) {
      const { name, type } = info;
      if (name.startsWith('__activeUniform')) {
        // name assigned by `gl.shim.ts`, ignore for checks
        continue;
      }
      if (name === 'gl_InstanceID') continue; // WebGL2 built-in
      if (name === 'gl_VertexID') continue; // WebGL2 built-in
      if (name === 'gl_DrawID') continue; // WEBGL_multi_draw built-in
      if (name === 'gl_ViewID_OVR') continue; // OVR_multiview2 built-in
      const baseName = name.replace(/[[0-9]+\]$/, ''); // 'array' uniforms
      const spec = schema[baseName];
      if (spec === undefined) {
        throw new Error(`missing 'uniform' or 'texture' with name '${name}' in schema`);
      }
      if (spec.type === 'uniform') {
        const uniformType = getUniformType(gl, spec.kind);
        if (uniformType !== type) {
          throw new Error(`unexpected uniform type for ${name}`);
        }
      } else if (spec.type === 'texture') {
        const target = getDimensionTarget(gl, getTextureDimension(spec, defines));
        if (getSamplerTarget(gl, type) !== target) {
          throw new Error(`sampler type of '${name}' does not match the spec dimension`);
        }
        const isInteger = spec.sampleType === 'sint' || spec.sampleType === 'uint';
        if (isWebGL2(gl) && isIntegerSampler(gl, type) !== isInteger) {
          throw new Error(`sampler type of '${name}' does not match the spec sample type '${spec.sampleType}'`);
        }
        if (!isUsed(spec, defines)) {
          throw new Error(`texture '${name}' is used, but its spec 'usedIf' condition does not hold`);
        }
      } else {
        throw new Error(`'${name}' must be of type 'uniform' or 'texture' but is '${spec.type}'`);
      }
    }
  }
}

/** sampler types of the active texture uniforms, as reported by the driver */
function getSamplerTypes(gl: GLRenderingContext, program: WebGLProgram, schema: RenderableSchema) {
  const samplerTypes: { [k: string]: number } = {};
  const count = gl.getProgramParameter(program, gl.ACTIVE_UNIFORMS);
  for (let i = 0; i < count; ++i) {
    const info = gl.getActiveUniform(program, i);
    if (info && schema[info.name]?.type === 'texture') samplerTypes[info.name] = info.type;
  }
  return samplerTypes;
}

function getSamplerTarget(gl: GLRenderingContext, samplerType: number) {
  if (samplerType === gl.SAMPLER_2D) return gl.TEXTURE_2D;
  if (samplerType === gl.SAMPLER_CUBE) return gl.TEXTURE_CUBE_MAP;
  if (isWebGL2(gl)) {
    switch (samplerType) {
      case gl.INT_SAMPLER_2D:
      case gl.UNSIGNED_INT_SAMPLER_2D:
      case gl.SAMPLER_2D_SHADOW:
        return gl.TEXTURE_2D;
      case gl.SAMPLER_3D:
      case gl.INT_SAMPLER_3D:
      case gl.UNSIGNED_INT_SAMPLER_3D:
        return gl.TEXTURE_3D;
      case gl.INT_SAMPLER_CUBE:
      case gl.UNSIGNED_INT_SAMPLER_CUBE:
      case gl.SAMPLER_CUBE_SHADOW:
        return gl.TEXTURE_CUBE_MAP;
      case gl.SAMPLER_2D_ARRAY:
      case gl.INT_SAMPLER_2D_ARRAY:
      case gl.UNSIGNED_INT_SAMPLER_2D_ARRAY:
      case gl.SAMPLER_2D_ARRAY_SHADOW:
        return gl.TEXTURE_2D_ARRAY;
    }
  }
  throw new Error(`unknown sampler type '${samplerType}'`);
}

function getDimensionTarget(gl: GLRenderingContext, dimension: TextureDimension) {
  switch (dimension) {
    case '2d':
      return gl.TEXTURE_2D;
    case '3d':
      return (gl as WebGL2RenderingContext).TEXTURE_3D;
    case 'cube':
      return gl.TEXTURE_CUBE_MAP;
  }
}

/** all `#define`d names, including local ones, approximates the defines a program was compiled with */
function getShaderDefineLookup(shaderCode: GlslSources): DefineLookup {
  const names = `${shaderCode.vert}\n${shaderCode.frag}`.matchAll(/^[ \t]*#define[ \t]+(\w+)/gm);
  return getDefineLookup(Array.from(names, (m) => m[1]));
}

function isIntegerSampler(gl: GLRenderingContext, samplerType: number) {
  if (!isWebGL2(gl)) return false;
  switch (samplerType) {
    case gl.INT_SAMPLER_2D:
    case gl.INT_SAMPLER_3D:
    case gl.INT_SAMPLER_CUBE:
    case gl.INT_SAMPLER_2D_ARRAY:
    case gl.UNSIGNED_INT_SAMPLER_2D:
    case gl.UNSIGNED_INT_SAMPLER_3D:
    case gl.UNSIGNED_INT_SAMPLER_CUBE:
    case gl.UNSIGNED_INT_SAMPLER_2D_ARRAY:
      return true;
  }
  return false;
}

function isIntegerFormat(gl: GLRenderingContext, format: number) {
  if (!isWebGL2(gl)) return false;
  return (
    format === gl.RED_INTEGER || format === gl.RG_INTEGER || format === gl.RGB_INTEGER || format === gl.RGBA_INTEGER
  );
}

/** the bound texture must match the sampler declared in the shader, otherwise sampling is undefined */
function checkTextureBinding(gl: GLRenderingContext, name: string, texture: GLTexture, samplerType: number) {
  const target = getSamplerTarget(gl, samplerType);
  if (texture.target !== target) {
    throw new Error(`texture '${name}' has target ${texture.target}, but the sampler requires ${target}`);
  }
  if (isIntegerSampler(gl, samplerType) !== isIntegerFormat(gl, texture.format)) {
    throw new Error(`texture '${name}' format ${texture.format} does not match the sampler (integer vs. float)`);
  }
}

/** active attributes must be at their fixed location, see `getAttributeLayout` */
function checkAttributeLocations(gl: GLRenderingContext, program: WebGLProgram, layout: AttributeLayout) {
  for (const { name, location } of layout.entries) {
    const actual = gl.getAttribLocation(program, name);
    if (actual !== -1 && actual !== location) {
      throw new Error(`attribute '${name}' is at location ${actual}, expected ${location}`);
    }
  }
}

function checkProgram(gl: GLRenderingContext, program: WebGLProgram) {
  // no-op in FF on Mac, see https://bugzilla.mozilla.org/show_bug.cgi?id=1284425
  // gl.validateProgram(program)
  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    throw new Error(`Could not compile WebGL program. \n\n${gl.getProgramInfoLog(program)}`);
  }
}

/** `undefined` (unlike an inactive `null`) location: not in schema or backed by a uniform buffer */
function throwMissingLocation(name: string): never {
  throw new Error(`uniform '${name}' has no location, it is either not in the schema or a uniform block member`);
}

const warnedClassicUniforms = new Set<string>();
/** only 'draw' uniforms are expected outside of uniform blocks, see `UniformSpec` */
function checkClassicUniforms(schema: RenderableSchema) {
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'uniform' || spec.kind === 't' || spec.block !== undefined) continue;
    if (warnedClassicUniforms.has(k)) continue;
    warnedClassicUniforms.add(k);
    console.warn(`uniform '${k}' is not in a uniform block, array uniforms need a 'length'`);
  }
}

/** checks that the CPU-computed std140 layouts match the layouts reported by the driver */
function checkUniformBlockLayouts(gl: WebGL2RenderingContext, program: WebGLProgram, schema: RenderableSchema) {
  for (const block of UniformBlocks) {
    const layout = getUniformBlockLayout(block, schema);
    if (!layout) continue;

    const { name } = UniformBlockInfo[block];
    const blockIndex = gl.getUniformBlockIndex(program, name);
    if (blockIndex === gl.INVALID_INDEX) continue;

    const dataSize: number = gl.getActiveUniformBlockParameter(program, blockIndex, gl.UNIFORM_BLOCK_DATA_SIZE);
    if (dataSize > layout.size) {
      throw new Error(`uniform block '${name}' data size mismatch, expected ${layout.size} got ${dataSize}`);
    }

    const indices = gl.getUniformIndices(
      program,
      layout.members.map((m) => (m.length ? `${m.name}[0]` : m.name)),
    );
    if (indices === null) continue;
    const validIndices: number[] = [];
    const validMembers: number[] = [];
    for (let i = 0; i < indices.length; ++i) {
      if (indices[i] !== gl.INVALID_INDEX) {
        validIndices.push(indices[i]);
        validMembers.push(i);
      }
    }
    if (validIndices.length === 0) continue;
    const offsets: number[] = gl.getActiveUniforms(program, validIndices, gl.UNIFORM_OFFSET);
    const strides: number[] = gl.getActiveUniforms(program, validIndices, gl.UNIFORM_ARRAY_STRIDE);
    for (let i = 0; i < validIndices.length; ++i) {
      const m = layout.members[validMembers[i]];
      if (offsets[i] !== m.offset) {
        throw new Error(
          `uniform block '${name}' offset mismatch for '${m.name}', expected ${m.offset} got ${offsets[i]}`,
        );
      }
      if (m.length && strides[i] !== m.stride) {
        throw new Error(
          `uniform block '${name}' array stride mismatch for '${m.name}', expected ${m.stride} got ${strides[i]}`,
        );
      }
    }
  }
}

export interface ProgramProps {
  defineValues: DefineValues;
  shaderCode: ShaderCode;
  schema: RenderableSchema;
}

export function getProgram(gl: GLRenderingContext) {
  const program = gl.createProgram();
  if (program === null) {
    throw new Error('Could not create WebGL program');
  }
  return program;
}

type ShaderGetter = (type: ShaderType, source: string) => Shader;

function normalizeVariant(variant: any): ProgramVariant {
  if (typeof variant !== 'string') throw new Error(`unknown program variant: ${variant}`);
  if (variant.startsWith('color')) return 'color';
  if (variant.startsWith('pick')) return 'pick';
  return variant as ProgramVariant;
}

export function createProgram(
  gl: GLRenderingContext,
  state: WebGLState,
  extensions: WebGLExtensions,
  parameters: WebGLParameters,
  getShader: ShaderGetter,
  props: ProgramProps,
): Program {
  const { defineValues, shaderCode: _shaderCode, schema } = props;

  let program = getProgram(gl);
  const programId = getNextProgramId();
  const variant = normalizeVariant(defineValues.dRenderVariant.ref.value);

  const useBlocks = isWebGL2(gl);
  const shaderCode = addShaderDefines(
    gl,
    extensions,
    parameters,
    defineValues,
    _shaderCode,
    getShaderDeclarations(schema, useBlocks, getDefineValuesLookup(defineValues)),
  );
  const vertShader = getShader('vert', shaderCode.vert);
  const fragShader = getShader('frag', shaderCode.frag);

  const attributeLayout = getAttributeLayout(schema);
  if (attributeLayout.size > parameters.maxVertexAttribs) {
    throw new Error(
      `'${_shaderCode.name}' needs ${attributeLayout.size} attribute locations, only ${parameters.maxVertexAttribs} are available`,
    );
  }

  let locations: Locations;
  let uniformSetters: UniformSetters;
  /** only set in debug mode */
  let samplerTypes: { [k: string]: number } | undefined;

  let linked = false;
  let finalized = false;
  let destroyed = false;

  function link() {
    vertShader.attach(program);
    fragShader.attach(program);
    for (const { name, location } of attributeLayout.entries) {
      gl.bindAttribLocation(program, location, name);
    }
    gl.linkProgram(program);
    if (isDebugMode) {
      checkProgram(gl, program);
    }
    linked = true;
  }
  if (variant === 'compute') link();

  function finalize() {
    locations = getLocations(gl, program, schema, useBlocks);
    uniformSetters = getUniformSetters(schema);
    if (useBlocks && isWebGL2(gl)) {
      for (const block of UniformBlocks) {
        const info = UniformBlockInfo[block];
        const blockIndex = gl.getUniformBlockIndex(program, info.name);
        if (blockIndex !== gl.INVALID_INDEX) {
          gl.uniformBlockBinding(program, blockIndex, info.bindingIndex);
        }
      }
    }
    if (isDebugMode) {
      const defines = getShaderDefineLookup(shaderCode);
      checkActiveAttributes(gl, program, schema, defines);
      checkAttributeLocations(gl, program, attributeLayout);
      checkActiveUniforms(gl, program, schema, defines);
      samplerTypes = getSamplerTypes(gl, program, schema);
      if (useBlocks && isWebGL2(gl)) {
        checkUniformBlockLayouts(gl, program, schema);
        checkClassicUniforms(schema);
      }
    }
    finalized = true;
  }

  return {
    id: programId,
    variant,

    isReady: () => {
      return finalized;
    },
    link: () => {
      if (!linked) link();
    },
    finalize(force?: boolean): boolean {
      if (!linked) link();
      if (
        !finalized &&
        (force ||
          !extensions.parallelShaderCompile ||
          gl.getProgramParameter(program, extensions.parallelShaderCompile.COMPLETION_STATUS))
      ) {
        finalize();
      }
      return finalized;
    },

    use: () => {
      // console.log('use', programId)
      if (isDebugMode && !finalized) throw new Error(`program not finalized: ${variant}`);
      state.currentProgramId = programId;
      gl.useProgram(program);
    },
    setUniforms: (uniformValues: UniformsList) => {
      for (let i = 0, il = uniformValues.length; i < il; ++i) {
        const [k, v] = uniformValues[i];
        if (v) {
          const l = locations[k];
          if (l === undefined) {
            if (isDebugMode) throwMissingLocation(k);
          } else if (l !== null) {
            uniformSetters[k](gl, l, v.ref.value);
          }
        }
      }
    },
    uniform: (k: string, v: UniformType) => {
      const l = locations[k];
      if (l === undefined) {
        if (isDebugMode) throwMissingLocation(k);
      } else if (l !== null) {
        uniformSetters[k](gl, l, v);
      }
    },
    bindAttributes: (attributeBuffers: AttributeBuffers) => {
      bindAttributeBuffers(state, attributeLayout, attributeBuffers);
    },
    offsetAttributes: (attributeBuffers: AttributeBuffers, offset) => {
      for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
        const [k, buffer] = attributeBuffers[i];
        buffer.changeOffset(attributeLayout.byName[k].location, offset);
      }
    },
    bindTextures: (textures: Textures, startingTargetUnit: number) => {
      for (let i = 0, il = textures.length; i < il; ++i) {
        const [k, texture] = textures[i];
        const l = locations[k];
        if (l !== null && l !== undefined) {
          const unit = (i + startingTargetUnit) as TextureId;
          if (isNullTexture(texture)) {
            gl.activeTexture(gl.TEXTURE0 + unit);
            gl.bindTexture(gl.TEXTURE_2D, null);
          } else {
            const glTexture = asGLTexture(texture);
            if (samplerTypes?.[k] !== undefined) checkTextureBinding(gl, k, glTexture, samplerTypes[k]);
            glTexture.bind(unit);
          }
          uniformSetters[k](gl, l, unit);
        }
      }
    },

    reset: () => {
      program = getProgram(gl);
      if (linked) link();
      if (finalized) finalize();
    },
    destroy: () => {
      if (destroyed) return;
      vertShader.destroy();
      fragShader.destroy();
      gl.deleteProgram(program);
      destroyed = true;
    },
  };
}
