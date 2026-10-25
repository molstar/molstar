/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import {
  type DefineLookup,
  type RenderableSchema,
  type TextureDimension,
  TextureDimensions,
  type TextureSampleType,
} from '../renderable/schema.js';
import type { ShaderUniformDeclarations } from './shader-code.js';
import {
  getElementKind,
  getUniformBlock,
  getUniformBlockLayout,
  UniformBlockInfo,
  type UniformBlockLayout,
  UniformBlocks,
} from '../uniform-layout.js';
import { getAttributeLayout } from '../attribute-layout.js';
import type { AttributeItemSize } from '../buffer.js';
import { getUniformGlslType } from './uniform.js';
import { getDataTextureAccessors } from './data-textures.js';

export function getUniformBlockDeclaration(layout: UniformBlockLayout) {
  const lines = [`layout(std140) uniform ${UniformBlockInfo[layout.block].name} {`];
  for (const m of layout.members) {
    const type = getUniformGlslType(getElementKind(m.kind));
    lines.push(m.length ? `    ${type} ${m.name}[${m.length}];` : `    ${type} ${m.name};`);
  }
  lines.push('};');
  return lines.join('\n');
}

/**
 * GLSL declarations of all uniforms of `schema`, the source of truth for shaders.
 * With `useBlocks`, blocks for block members and classic declarations for the rest,
 * otherwise classic declarations only. Textures are not included. Classic arrays are sized
 * by their `lengthDefine` unless it is a uniform (or `length`), arrays without either are left to the shader code.
 */
export function getShaderUniformDeclarations(
  schema: RenderableSchema,
  useBlocks: boolean,
): ShaderUniformDeclarations | undefined {
  const cache = useBlocks ? blockDeclarationsCache : classicDeclarationsCache;
  if (cache.has(schema)) return cache.get(schema);

  const lines: string[] = [
    // precision must match across stages and in GLSL ES 1.0 fragment shaders
    // it must be set before declarations, the shader body's statements come later
    'precision highp float;',
    'precision highp int;',
  ];
  const types = new Map<string, string>();
  if (useBlocks) {
    for (const block of UniformBlocks) {
      const layout = getUniformBlockLayout(block, schema);
      if (!layout) continue;
      lines.push(getUniformBlockDeclaration(layout));
      for (const m of layout.members) types.set(m.name, getUniformGlslType(getElementKind(m.kind)));
    }
  }
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'uniform' || spec.kind === 't') continue;
    if (useBlocks && getUniformBlock(spec) !== undefined) continue;
    const type = getUniformGlslType(getElementKind(spec.kind));
    if (!spec.kind.endsWith('[]')) {
      lines.push(`uniform ${type} ${k};`);
    } else if (spec.lengthDefine !== undefined && schema[spec.lengthDefine]?.type !== 'uniform') {
      // zero-length arrays don't compile
      lines.push(`#if ${spec.lengthDefine} > 0`, `uniform ${type} ${k}[${spec.lengthDefine}];`, '#endif');
    } else if (spec.length !== undefined) {
      // the maximum, where the length is a uniform (a loop bound, a define on WebGL1 only)
      lines.push(`uniform ${type} ${k}[${spec.length}];`);
    } else {
      continue;
    }
    types.set(k, type);
  }
  const declarations = types.size > 0 ? { code: lines.join('\n'), types } : undefined;
  cache.set(schema, declarations);
  return declarations;
}

const blockDeclarationsCache = new WeakMap<RenderableSchema, ShaderUniformDeclarations | undefined>();
const classicDeclarationsCache = new WeakMap<RenderableSchema, ShaderUniformDeclarations | undefined>();

//

/** `undefined` if not available (3D textures on WebGL1); integer samples are read as float on WebGL1 */
export function getSamplerGlslType(
  sampleType: TextureSampleType,
  dimension: TextureDimension,
  isWebGL2: boolean,
): string | undefined {
  if (dimension === '3d' && !isWebGL2) return;
  const base = dimension === '2d' ? 'sampler2D' : dimension === '3d' ? 'sampler3D' : 'samplerCube';
  if (!isWebGL2) return base;
  return sampleType === 'sint' ? `i${base}` : sampleType === 'uint' ? `u${base}` : base;
}

export function getAttributeGlslType(itemSize: AttributeItemSize) {
  switch (itemSize) {
    case 1:
      return 'float';
    case 2:
      return 'vec2';
    case 3:
      return 'vec3';
    case 4:
      return 'vec4';
    case 16:
      return 'mat4';
  }
}

export type ShaderDeclarations = {
  /** generated declarations, injected at the start of the vertex/fragment shader */
  readonly vert: string;
  readonly frag: string;
  /** generated uniforms, hand-written declarations are removed, see `stripUniformDeclarations` */
  readonly uniformTypes: ReadonlyMap<string, string>;
  /** hand-written declarations with one of these types are removed (other types fail), see `stripSamplerDeclarations` */
  readonly samplerTypes: ReadonlyMap<string, readonly string[]>;
  /** generated attributes, hand-written declarations are removed, see `stripAttributeDeclarations` */
  readonly attributeTypes: ReadonlyMap<string, string>;
};

/**
 * GLSL declarations of the schema's uniforms (see `getShaderUniformDeclarations`), samplers, the data textures'
 * accessors (`getDataTextureAccessors`) and attributes (vertex shader only; with fixed locations, see
 * `getAttributeLayout`).
 * Define-dependent sampler dimensions are resolved with `defines`. Uniform blocks with `useBlocks` (WebGL2 only).
 */
export function getShaderDeclarations(
  schema: RenderableSchema,
  isWebGL2: boolean,
  defines: DefineLookup,
  useBlocks = isWebGL2,
): ShaderDeclarations {
  const samplers: [name: string, type: string, handWritten: string[]][] = [];
  let key = `${isWebGL2}|${useBlocks}`;
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'texture') continue;
    const dimensions = typeof spec.dimension === 'string' ? [spec.dimension] : TextureDimensions;
    const dimension =
      typeof spec.dimension === 'string'
        ? spec.dimension
        : TextureDimensions.find((d) => defines(`${(spec.dimension as { define: string }).define}_${d}`));
    if (typeof spec.dimension !== 'string') key += `|${k}:${dimension}`;
    const type = dimension && getSamplerGlslType(spec.sampleType, dimension, isWebGL2);
    if (!type) continue;
    // hand-written alternatives for other defines or WebGL1 are in inactive branches
    const handWritten = new Set<string>();
    for (const d of dimensions) {
      for (const webgl2 of [true, false]) {
        const t = getSamplerGlslType(spec.sampleType, d, webgl2);
        if (t) handWritten.add(t);
      }
    }
    samplers.push([k, type, Array.from(handWritten)]);
  }

  let cache = declarationsCache.get(schema);
  if (!cache) {
    cache = new Map();
    declarationsCache.set(schema, cache);
  }
  const cached = cache.get(key);
  if (cached) return cached;

  const uniforms = getShaderUniformDeclarations(schema, isWebGL2 && useBlocks);
  const shared = uniforms ? [uniforms.code] : [];
  const samplerTypes = new Map<string, readonly string[]>();
  for (const [name, type, handWritten] of samplers) {
    shared.push(`uniform highp ${type} ${name};`);
    samplerTypes.set(name, handWritten);
  }
  // the accessors of the data textures (`read_<name>(i)`), after their samplers and size uniforms
  const accessors = getDataTextureAccessors(schema, isWebGL2);
  if (accessors) shared.push(accessors);

  const attributes: string[] = [];
  const attributeTypes = new Map<string, string>();
  for (const { name, itemSize, location } of getAttributeLayout(schema).entries) {
    const type = getAttributeGlslType(itemSize);
    attributes.push(isWebGL2 ? `layout(location = ${location}) in ${type} ${name};` : `attribute ${type} ${name};`);
    attributeTypes.set(name, type);
  }

  const declarations: ShaderDeclarations = {
    vert: [...shared, ...attributes].join('\n'),
    frag: shared.join('\n'),
    uniformTypes: uniforms?.types ?? new Map(),
    samplerTypes,
    attributeTypes,
  };
  cache.set(key, declarations);
  return declarations;
}

const declarationsCache = new WeakMap<RenderableSchema, Map<string, ShaderDeclarations>>();
