/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import type { AttributeLayoutEntry } from '../attribute-layout.js';
import { type DefineLookup, type RenderableSchema, type TextureSampleType } from '../renderable/schema.js';
import type { DefineValues } from '../shader-code.js';
import { preprocessShader } from '../shader-preprocessor.js';
import { getAttributeLayout } from '../attribute-layout.js';
import {
  type BindingLayout,
  getBindingLayout,
  getUniformBlockLayout,
  type UniformBlockMember,
  UniformBlocks,
} from '../uniform-layout.js';
import { ShaderStage } from './constants.js';
import { DataLayerPagesName, getDataLayerAccessors, getReadDataPages } from './data-layers.js';
import {
  getWgslDrawBinding,
  getWgslDrawLayout,
  getWgslStorageBindings,
  type WgslStorageBinding,
} from './uniform-declarations.js';
import {
  getWgslDeclarations,
  getWgslGroup,
  getWgslSamplerBinding,
  getWgslVertexAttributes,
  WgslSamplers,
} from './wgsl-declarations.js';

const getNextProgramId = idFactory();

/** Preprocessor defines of the define values: string `dX: 'v'` → `dX_v`, numbers with their value, `true` booleans */
export function getWgslDefines(defineValues: DefineValues): Map<string, string> {
  const defines = new Map<string, string>();
  for (const name of Object.keys(defineValues)) {
    const v = defineValues[name].ref.value as unknown;
    if (typeof v === 'string') defines.set(`${name}_${v}`, '');
    else if (typeof v === 'number') defines.set(name, `${v}`);
    else if (v === true) defines.set(name, '');
  }
  return defines;
}

export function getWgslDefineLookup(defines: ReadonlyMap<string, string>): DefineLookup {
  return (name) => defines.has(name);
}

/**
 * Override constants declared by `source` (`override oX`), with the values of the number defines `dX`; with
 * `entryPoint`, only those the entry point reaches (WebKit rejects constants for overrides a stage does not use)
 */
export function getOverrideConstants(source: string, defineValues: DefineValues, entryPoint?: string) {
  const reached = entryPoint !== undefined ? getReachedDeclarations(source, entryPoint) : undefined;
  const constants: Record<string, number> = {};
  for (const m of source.matchAll(/\boverride\s+(o\w+)/g)) {
    if (reached && !reached.has(m[1])) continue;
    const value = defineValues[`d${m[1].slice(1)}`]?.ref.value as unknown;
    if (typeof value === 'number') constants[m[1]] = value;
  }
  return constants;
}

/** Module-scope declarations of `source` by name (functions, structs, constants, variables, overrides, aliases) */
function getModuleDeclarations(source: string) {
  const code = source.replace(/\/\/[^\n]*|\/\*[\s\S]*?\*\//g, ' ');
  const declarations = new Map<string, string>();
  const keyword = /(fn|struct|const|var|override|alias)\b(?:\s*<[^>]*>)?\s+(\w+)/y;
  let depth = 0;
  for (let i = 0; i < code.length; ++i) {
    const c = code[i];
    if (c === '{') ++depth;
    else if (c === '}') --depth;
    if (depth > 0 || (i > 0 && /\w/.test(code[i - 1]))) continue;
    keyword.lastIndex = i;
    const m = keyword.exec(code);
    if (!m) continue;
    const start = keyword.lastIndex;
    let end =
      m[1] === 'fn' || m[1] === 'struct' ? getBlockEnd(code, code.indexOf('{', start)) : code.indexOf(';', start);
    if (end < 0) end = code.length;
    declarations.set(m[2], code.slice(start, end + 1));
    i = end;
  }
  return declarations;
}

/** Index of the brace closing the block that opens at `open` */
function getBlockEnd(code: string, open: number) {
  if (open < 0) return -1;
  let depth = 0;
  for (let i = open; i < code.length; ++i) {
    if (code[i] === '{') ++depth;
    else if (code[i] === '}' && --depth === 0) return i;
  }
  return -1;
}

/** Names of the module-scope declarations that `entryPoint` uses, directly or through other declarations */
function getReachedDeclarations(source: string, entryPoint: string) {
  const declarations = getModuleDeclarations(source);
  const reached = new Set([entryPoint]);
  const queue = [entryPoint];
  while (queue.length > 0) {
    const text = declarations.get(queue.pop()!);
    if (!text) continue;
    for (const [name] of text.matchAll(/\b[A-Za-z_]\w*\b/g)) {
      if (declarations.has(name) && !reached.has(name)) {
        reached.add(name);
        queue.push(name);
      }
    }
  }
  return reached;
}

/**
 * Textures only read with `textureLoad` are bound as `unfilterable-float`, which also accepts unfilterable float and
 * depth formats (e.g. a depth texture copied like a color texture); filtered ones need `float`
 */
export function getTextureSampleType(sampleType: TextureSampleType, isFiltered: boolean): GPUTextureSampleType {
  return sampleType === 'float' && !isFiltered ? 'unfilterable-float' : sampleType;
}

/** Whether `source` samples `name` with `textureSample*` (and a filtering sampler) */
export function isFilteredTexture(source: string, name: string) {
  return new RegExp(`\\btextureSample\\w*\\(\\s*${name}\\b`).test(source);
}

const Visibility = ShaderStage.VERTEX | ShaderStage.FRAGMENT;

function getVertexFormat(itemSize: number): GPUVertexFormat {
  switch (itemSize) {
    case 1:
      return 'float32';
    case 2:
      return 'float32x2';
    case 3:
      return 'float32x3';
    case 4:
      return 'float32x4';
    default:
      throw new Error(`unsupported attribute item size '${itemSize}'`);
  }
}

/** Attributes read from one vertex buffer, interleaved if more than one */
export interface VertexBufferSlot {
  readonly key: string;
  readonly names: readonly string[];
  readonly itemSizes: readonly number[];
}

/**
 * One vertex buffer per attribute; with more attributes than `maxVertexBuffers`, the per-vertex attributes share one
 * interleaved buffer (cylinders have 9)
 */
export function getVertexBufferSlots(
  attributes: readonly AttributeLayoutEntry[],
  maxVertexBuffers: number,
): { slots: VertexBufferSlot[]; layouts: GPUVertexBufferLayout[] } {
  if (attributes.length <= maxVertexBuffers) {
    return {
      slots: attributes.map((a) => ({ key: a.name, names: [a.name], itemSizes: [a.itemSize] })),
      layouts: attributes.map(getVertexBufferLayout),
    };
  }
  const merged = attributes.filter((a) => !a.divisor && a.locationCount === 1);
  const separate = attributes.filter((a) => !merged.includes(a));
  if (separate.length + 1 > maxVertexBuffers) {
    throw new Error(`${attributes.length} attributes don't fit into ${maxVertexBuffers} vertex buffers`);
  }
  let offset = 0;
  const attrs: GPUVertexAttribute[] = [];
  for (const a of merged) {
    attrs.push({ shaderLocation: a.location, offset, format: getVertexFormat(a.itemSize) });
    offset += a.itemSize * 4;
  }
  return {
    slots: [
      ...separate.map((a) => ({ key: a.name, names: [a.name], itemSizes: [a.itemSize] })),
      {
        key: merged.map((a) => a.name).join('+'),
        names: merged.map((a) => a.name),
        itemSizes: merged.map((a) => a.itemSize),
      },
    ],
    layouts: [...separate.map(getVertexBufferLayout), { arrayStride: offset, stepMode: 'vertex', attributes: attrs }],
  };
}

export function getVertexBufferLayout(entry: AttributeLayoutEntry): GPUVertexBufferLayout {
  const stepMode: GPUVertexStepMode = entry.divisor ? 'instance' : 'vertex';
  if (entry.locationCount === 4) {
    return {
      arrayStride: 64,
      stepMode,
      attributes: [0, 1, 2, 3].map((i) => ({
        shaderLocation: entry.location + i,
        offset: i * 16,
        format: 'float32x4' as GPUVertexFormat,
      })),
    };
  }
  return {
    arrayStride: entry.itemSize * 4,
    stepMode,
    attributes: [{ shaderLocation: entry.location, offset: 0, format: getVertexFormat(entry.itemSize) }],
  };
}

export interface WebGPUProgram {
  readonly id: number;
  /** the name of the `ShaderCode`, for labels */
  readonly name: string;
  /** the render variant (`dRenderVariant`), e.g. `colorBlended` or `compute`, for labels */
  readonly variant: string;
  readonly module: GPUShaderModule;
  readonly code: string;
  readonly pipelineLayout: GPUPipelineLayout;
  /** one per group, empty where a group has no bindings */
  readonly bindGroupLayouts: readonly GPUBindGroupLayout[];
  readonly bindingLayout: BindingLayout;
  /** textures, samplers and storage buffers in the bind group layouts */
  readonly referenced: ReadonlySet<string>;
  /** the storage buffers the program reads, in the object's bind group */
  readonly storageBindings: readonly WgslStorageBinding[];
  readonly drawLayout: { members: UniformBlockMember[]; size: number } | undefined;
  /** binding of the draw uniforms in the object's bind group, see `getWgslDrawBinding` */
  readonly drawBinding: number;
  readonly attributes: readonly AttributeLayoutEntry[];
  /** the attributes per vertex buffer, see `getVertexBufferSlots` */
  readonly vertexBufferSlots: readonly VertexBufferSlot[];
  readonly vertexBuffers: readonly GPUVertexBufferLayout[];
  /** override constants of the vertex and fragment stage, each only those the stage uses */
  readonly vertexConstants: Record<string, number>;
  readonly fragmentConstants: Record<string, number>;
  readonly pipelines: Map<string, GPURenderPipeline>;
  /** the pipeline of the last lookup (`getPipelineFor`): draws of a program one after another mostly ask the same */
  lastPipeline?: { state: object; target: object; drawMode: string; pipeline: GPURenderPipeline };
}

export function createWebGPUProgram(
  device: GPUDevice,
  name: string,
  source: string,
  defineValues: DefineValues,
  schema: RenderableSchema,
  onError: (message: string) => void,
): WebGPUProgram {
  const defines = getWgslDefines(defineValues);
  const lookup = getWgslDefineLookup(defines);
  const body = preprocessShader(source, defines, { expandMacros: true });

  const identifiers = new Set(body.match(/\b[A-Za-z_]\w*\b/g) ?? []);
  // the data layers' pages are read by the generated accessors (`getDataLayerAccessors`), each layer from its page
  const dataPages = Number(defineValues[DataLayerPagesName]?.ref.value ?? 0);
  const readPages = getReadDataPages(schema, (n) => identifiers.has(n), dataPages);
  const isReferenced = (n: string) => identifiers.has(n) || readPages.has(n);

  const bindingLayout = getBindingLayout(schema, lookup);
  const declarations = getWgslDeclarations(schema, bindingLayout, lookup, isReferenced);
  if (declarations.unsupported.length > 0) {
    throw new Error(`uniforms without a block are not supported by WebGPU: ${declarations.unsupported.join(', ')}`);
  }
  const accessors = getDataLayerAccessors(schema, isReferenced, dataPages);
  const code = accessors ? `${declarations.code}\n\n${accessors}\n\n${body}` : `${declarations.code}\n\n${body}`;

  const module = device.createShaderModule({ code, label: name });
  module.getCompilationInfo().then((info) => {
    const lines = code.split('\n');
    for (const m of info.messages) {
      if (m.type !== 'error') continue;
      onError(`WGSL error in '${name}' line ${m.lineNum}: ${m.message}\n  ${lines[m.lineNum - 1]?.trim()}`);
    }
  });

  const referenced = new Set<string>();
  const groups: GPUBindGroupLayoutEntry[][] = [[], [], []];
  for (const block of UniformBlocks) {
    const tier = bindingLayout[block];
    const entries = groups[getWgslGroup(block)];
    if (tier.uniforms) {
      entries.push({
        binding: 0,
        visibility: Visibility,
        buffer: { type: 'uniform', hasDynamicOffset: true, minBindingSize: tier.uniforms.size },
      });
    }
    for (const t of tier.textures) {
      if (!isReferenced(t.name)) continue;
      referenced.add(t.name);
      entries.push({
        binding: t.binding,
        visibility: Visibility,
        texture: {
          sampleType: getTextureSampleType(t.sampleType, isFilteredTexture(body, t.name)),
          viewDimension: t.dimension,
        },
      });
    }
  }
  const samplerBinding = getWgslSamplerBinding(bindingLayout);
  for (let i = 0; i < WgslSamplers.length; ++i) {
    if (!isReferenced(WgslSamplers[i])) continue;
    referenced.add(WgslSamplers[i]);
    groups[0].push({ binding: samplerBinding + i, visibility: Visibility, sampler: { type: 'filtering' } });
  }
  const drawLayout = getWgslDrawLayout(schema);
  const drawBinding = getWgslDrawBinding(schema);
  if (drawLayout) {
    groups[getWgslGroup('object')].push({
      binding: drawBinding,
      visibility: Visibility,
      buffer: { type: 'uniform', hasDynamicOffset: true, minBindingSize: drawLayout.size },
    });
  }
  const storageBindings = getWgslStorageBindings(schema, lookup).filter((b) => isReferenced(b.name));
  for (const b of storageBindings) {
    referenced.add(b.name);
    groups[getWgslGroup('object')].push({
      binding: b.binding,
      visibility: Visibility,
      buffer: { type: 'read-only-storage' },
    });
  }

  const bindGroupLayouts = groups.map((entries, i) =>
    device.createBindGroupLayout({ entries, label: `${name} group ${i}` }),
  );
  const pipelineLayout = device.createPipelineLayout({ bindGroupLayouts, label: name });

  const attributes = getWgslVertexAttributes(schema, lookup);
  const vertexBuffers = getVertexBufferSlots(attributes, device.limits.maxVertexBuffers);

  return {
    id: getNextProgramId(),
    name,
    variant: String(defineValues.dRenderVariant?.ref.value ?? ''),
    module,
    code,
    pipelineLayout,
    bindGroupLayouts,
    bindingLayout,
    referenced,
    storageBindings,
    drawLayout,
    drawBinding,
    attributes,
    vertexBufferSlots: vertexBuffers.slots,
    vertexBuffers: vertexBuffers.layouts,
    vertexConstants: getOverrideConstants(body, defineValues, 'vertexMain'),
    fragmentConstants: getOverrideConstants(body, defineValues, 'fragmentMain'),
    pipelines: new Map(),
  };
}

const schemaKeys = new WeakMap<RenderableSchema, string>();

/** The parts of `schema` that generated declarations and layouts depend on, equal for schemas of the same shape */
function getSchemaKey(schema: RenderableSchema) {
  let key = schemaKeys.get(schema);
  if (key === undefined) {
    const parts = [getAttributeLayout(schema).hash];
    for (const block of UniformBlocks) parts.push(getUniformBlockLayout(block, schema)?.hash ?? '');
    const draw = getWgslDrawLayout(schema);
    parts.push(draw ? draw.members.map((m) => `${m.name}:${m.kind}`).join(',') : '');
    for (const k of Object.keys(schema)) {
      const spec = schema[k];
      if (spec.type === 'storage-buffer') {
        parts.push(`${k}:storage:${spec.element}:${spec.usedIf?.join('|') ?? ''}`);
        continue;
      }
      if (spec.type !== 'texture') continue;
      const dimension = typeof spec.dimension === 'string' ? spec.dimension : spec.dimension.define;
      parts.push(`${k}:${dimension}:${spec.sampleType}:${spec.variant ?? ''}:${spec.usedIf?.join('|') ?? ''}`);
    }
    key = parts.join(';');
    schemaKeys.set(schema, key);
  }
  return key;
}

/**
 * Programs (shader modules, layouts and their pipelines) shared by the render items of a context, by program
 * name, defines and schema shape
 */
export function createWebGPUProgramCache(
  device: GPUDevice,
  onError: (message: string) => void,
  onCreate?: (program: WebGPUProgram) => void,
) {
  const programs = new Map<string, WebGPUProgram>();
  return {
    get: (name: string, source: string, defineValues: DefineValues, schema: RenderableSchema) => {
      const defines = [...getWgslDefines(defineValues)].map(([k, v]) => (v ? `${k}=${v}` : k)).sort();
      const key = `${name}|${defines.join(',')}|${getSchemaKey(schema)}`;
      let program = programs.get(key);
      if (!program) {
        program = createWebGPUProgram(device, name, source, defineValues, schema, onError);
        programs.set(key, program);
        onCreate?.(program);
      }
      return program;
    },
    get size() {
      return programs.size;
    },
    clear: () => programs.clear(),
  };
}
export type WebGPUProgramCache = ReturnType<typeof createWebGPUProgramCache>;
