/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ShaderCode, DefineValues } from '../shader-code.js';
import { addShaderDefines, type GlslSources } from './shader-code.js';
import type { WebGLState } from './state.js';
import type { WebGLExtensions } from './extensions.js';
import { getUniformSetters, getUniformType, type UniformSetters } from './uniform.js';
import { isArrayUniform, type UniformsList, type UniformType } from '../uniform.js';
import {
  getBindingLayout,
  getUniformBlockLayout,
  UniformBlockInfo,
  UniformBlocks,
  getUniformBlock,
} from '../uniform-layout.js';
import { getShaderDeclarations } from './shader-declarations.js';
import { type AttributeBuffers, bindAttributeBuffers, getAttribType } from './buffer.js';
import { type AttributeLayout, getAttributeLayout } from '../attribute-layout.js';
import { asGLTexture, type GLTexture, type TextureId } from './texture.js';
import { isDataPageTexture } from './data-pages.js';
import { isNullTexture, type Texture, type Textures } from '../texture.js';
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
import { arrayEqual } from '@molstar/core/util/array';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import type { ShaderType, Shader } from './shader.js';
import type { WebGLParameters } from './context.js';
import type { GraphicsProgram, ProgramVariant } from '../render-item.js';

const getNextProgramId = idFactory();

export interface Program extends GraphicsProgram {
  /**
   * Whether linking started by `link` has finished, without blocking (`KHR_parallel_shader_compile`); true without the
   * extension (linking then blocks in `finalize`) and once destroyed
   */
  isCompiled: () => boolean;
  use: () => void;
  setUniforms: (uniformValues: UniformsList) => void;
  uniform: (k: string, v: UniformType) => void;
  bindAttributes: (attribueBuffers: AttributeBuffers) => void;
  offsetAttributes: (attributeBuffers: AttributeBuffers, offset: number) => void;
  /** Whether the program reads any of `attributeBuffers` (once finalized): else `offsetAttributes` changes nothing */
  readsAnyAttribute: (attributeBuffers: AttributeBuffers) => boolean;
  /** Binds `textures` at their units, fixed per program (see `getTextureUnits`) */
  bindTextures: (textures: Textures) => void;
  /** Binds no texture at the units of the samplers `names` the program reads */
  unbindTextures: (names: readonly string[]) => void;

  destroy: () => void;
}

export type Programs = { [k: string]: Program };

type Locations = { [k: string]: number };

type TextureBindings = {
  readonly units: (TextureId | undefined)[];
  readonly textures: (Texture | undefined)[];
  /** the GL texture, `null` for a null texture */
  readonly resolved: (GLTexture | null)[];
  /** whether a data page (`data-pages.ts`), whose binding is skipped where the unit has it */
  readonly paged: boolean[];
};

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
/**
 * Texture units of a program: a running index over the active samplers in the global, material and object tiers of
 * the binding layout (the order of the WebGPU bind groups), so units are fixed per program, its sampler uniforms are set
 * once, and the units stay as low as the samplers allow (devices may have only 8)
 */
function getTextureUnits(schema: RenderableSchema, defines: DefineLookup, isActive: (name: string) => boolean) {
  const layout = getBindingLayout(schema, defines);
  const units: { [k: string]: TextureId } = {};
  let unit = 0;
  for (const tier of UniformBlocks) {
    for (const { name } of layout[tier].textures) {
      if (isActive(name)) units[name] = unit++ as TextureId;
    }
  }
  // active samplers the layout leaves out (none unless a `usedIf` is wrong, see `checkActiveUniforms`) after it
  for (const k of Object.keys(schema)) {
    if (schema[k].type === 'texture' && isActive(k) && units[k] === undefined) units[k] = unit++ as TextureId;
  }
  return units;
}

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
  useBlocks: boolean,
): Program {
  const { defineValues, shaderCode: _shaderCode, schema } = props;

  const program = getProgram(gl);
  const programId = getNextProgramId();
  const variant = normalizeVariant(defineValues.dRenderVariant.ref.value);

  const shaderCode = addShaderDefines(
    gl,
    extensions,
    parameters,
    defineValues,
    _shaderCode,
    getShaderDeclarations(schema, isWebGL2(gl), getDefineValuesLookup(defineValues), useBlocks),
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
  /**
   * the layout's attributes the program reads: the others are not offset per draw (e.g. the instanced attributes of
   * merged render items reading their instance data from a texture, `dInstanceData`)
   */
  let activeAttributes: Set<string>;
  let uniformSetters: UniformSetters;
  /**
   * The value last set per classic uniform (a GL program keeps them): an unchanged one is not set again, e.g. the
   * draw-scope `uLod`, the same for most objects (`uniform4fv` was 14% of the frames of 1000 unmerged objects)
   */
  const uploaded: { [k: string]: unknown } = Object.create(null);
  function setUniform(k: string, l: number, v: UniformType) {
    const prev = uploaded[k];
    if (typeof v === 'object') {
      if (prev !== undefined && arrayEqual(prev as ArrayLike<number>, v as ArrayLike<number>)) return;
      uploaded[k] = (v as unknown as { slice(): unknown }).slice();
    } else {
      if (prev === v) return;
      uploaded[k] = v;
    }
    uniformSetters[k](gl, l, v);
  }
  let textureUnits: { [k: string]: TextureId };
  /**
   * per textures list the units of its entries and their textures resolved (`isNullTexture`, `asGLTexture` are set
   * lookups), checked by identity per draw: the textures of a list may be replaced in place
   */
  let textureBindings = new WeakMap<Textures, TextureBindings>();
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

  /** Assigns the texture units and sets the sampler uniforms to them, leaving the current program as it was */
  function setTextureUnits(defines: DefineLookup) {
    // textures of inactive samplers get no unit and are not bound
    textureUnits = getTextureUnits(schema, defines, (k) => locations[k] !== null && locations[k] !== undefined);
    textureBindings = new WeakMap();
    const current = gl.getParameter(gl.CURRENT_PROGRAM) as WebGLProgram | null;
    gl.useProgram(program);
    for (const k of Object.keys(textureUnits)) gl.uniform1i(locations[k], textureUnits[k]);
    gl.useProgram(current);
  }

  function finalize() {
    locations = getLocations(gl, program, schema, useBlocks);
    activeAttributes = new Set();
    for (const { name } of attributeLayout.entries) {
      if (gl.getAttribLocation(program, name) !== -1) activeAttributes.add(name);
    }
    uniformSetters = getUniformSetters(schema);
    setTextureUnits(getShaderDefineLookup(shaderCode));
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
    isCompiled: () => {
      if (finalized || destroyed || !extensions.parallelShaderCompile) return true;
      return linked && !!gl.getProgramParameter(program, extensions.parallelShaderCompile.COMPLETION_STATUS);
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
            setUniform(k, l, v.ref.value);
          }
        }
      }
    },
    uniform: (k: string, v: UniformType) => {
      const l = locations[k];
      if (l === undefined) {
        if (isDebugMode) throwMissingLocation(k);
      } else if (l !== null) {
        setUniform(k, l, v);
      }
    },
    bindAttributes: (attributeBuffers: AttributeBuffers) => {
      bindAttributeBuffers(state, attributeLayout, attributeBuffers);
    },
    offsetAttributes: (attributeBuffers: AttributeBuffers, offset) => {
      for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
        const [k, buffer] = attributeBuffers[i];
        if (activeAttributes.has(k)) buffer.changeOffset(attributeLayout.byName[k].location, offset);
      }
    },
    readsAnyAttribute: (attributeBuffers: AttributeBuffers) => {
      for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
        if (activeAttributes.has(attributeBuffers[i][0])) return true;
      }
      return false;
    },
    unbindTextures: (names: readonly string[]) => {
      for (const k of names) {
        const unit = textureUnits[k];
        if (unit === undefined) continue;
        state.boundTextures[unit] = undefined;
        gl.activeTexture(gl.TEXTURE0 + unit);
        gl.bindTexture(gl.TEXTURE_2D, null);
      }
    },
    bindTextures: (textures: Textures) => {
      let b = textureBindings.get(textures);
      if (!b) {
        b = { units: textures.map(([k]) => textureUnits[k]), textures: [], resolved: [], paged: [] };
        textureBindings.set(textures, b);
      }
      const { units, resolved, paged } = b;
      const bound = state.boundTextures;
      for (let i = 0, il = textures.length; i < il; ++i) {
        const unit = units[i];
        if (unit === undefined) continue;
        const texture = textures[i][1];
        if (b.textures[i] !== texture) {
          b.textures[i] = texture;
          resolved[i] = isNullTexture(texture) ? null : asGLTexture(texture);
          paged[i] = isDataPageTexture(texture);
        }
        const glTexture = resolved[i];
        if (glTexture === null) {
          bound[unit] = undefined;
          gl.activeTexture(gl.TEXTURE0 + unit);
          gl.bindTexture(gl.TEXTURE_2D, null);
        } else {
          const k = textures[i][0];
          if (samplerTypes?.[k] !== undefined) checkTextureBinding(gl, k, glTexture, samplerTypes[k]);
          // data pages only, which are never render targets: skipping others left a target bound while drawn into
          // (a feedback loop, `transparency/dpoit/overlap`)
          if (paged[i]) {
            if (bound[unit] === glTexture) continue;
            bound[unit] = glTexture;
          } else {
            bound[unit] = undefined;
          }
          glTexture.bind(unit);
        }
      }
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
