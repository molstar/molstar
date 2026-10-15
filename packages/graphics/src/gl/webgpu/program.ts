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
import { getWgslDrawLayout, WgslDrawGroup } from './uniform-declarations.js';
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

/** Override constants declared by `source` (`override oX`), with the values of the number defines `dX` */
export function getOverrideConstants(source: string, defineValues: DefineValues) {
  const constants: Record<string, number> = {};
  for (const m of source.matchAll(/\boverride\s+(o\w+)/g)) {
    const value = defineValues[`d${m[1].slice(1)}`]?.ref.value as unknown;
    if (typeof value === 'number') constants[m[1]] = value;
  }
  return constants;
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
  readonly module: GPUShaderModule;
  readonly code: string;
  readonly pipelineLayout: GPUPipelineLayout;
  /** one per group, empty where a group has no bindings */
  readonly bindGroupLayouts: readonly GPUBindGroupLayout[];
  readonly bindingLayout: BindingLayout;
  /** textures and samplers in the bind group layouts */
  readonly referenced: ReadonlySet<string>;
  readonly drawLayout: { members: UniformBlockMember[]; size: number } | undefined;
  readonly attributes: readonly AttributeLayoutEntry[];
  /** the attributes per vertex buffer, see `getVertexBufferSlots` */
  readonly vertexBufferSlots: readonly VertexBufferSlot[];
  readonly vertexBuffers: readonly GPUVertexBufferLayout[];
  readonly constants: Record<string, number>;
  readonly pipelines: Map<string, GPURenderPipeline>;
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
  const isReferenced = (n: string) => identifiers.has(n);

  const bindingLayout = getBindingLayout(schema, lookup);
  const declarations = getWgslDeclarations(schema, bindingLayout, lookup, isReferenced);
  if (declarations.unsupported.length > 0) {
    throw new Error(`uniforms without a block are not supported by WebGPU: ${declarations.unsupported.join(', ')}`);
  }
  const code = `${declarations.code}\n\n${body}`;

  const module = device.createShaderModule({ code, label: name });
  module.getCompilationInfo().then((info) => {
    const lines = code.split('\n');
    for (const m of info.messages) {
      if (m.type !== 'error') continue;
      onError(`WGSL error in '${name}' line ${m.lineNum}: ${m.message}\n  ${lines[m.lineNum - 1]?.trim()}`);
    }
  });

  const referenced = new Set<string>();
  const groups: GPUBindGroupLayoutEntry[][] = [[], [], [], []];
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
  if (drawLayout) {
    groups[WgslDrawGroup].push({
      binding: 0,
      visibility: Visibility,
      buffer: { type: 'uniform', hasDynamicOffset: true, minBindingSize: drawLayout.size },
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
    module,
    code,
    pipelineLayout,
    bindGroupLayouts,
    bindingLayout,
    referenced,
    drawLayout,
    attributes,
    vertexBufferSlots: vertexBuffers.slots,
    vertexBuffers: vertexBuffers.layouts,
    constants: getOverrideConstants(body, defineValues),
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
/** `ref.device` is read when a program is created, it changes when the context restores a lost device */
export function createWebGPUProgramCache(
  ref: { readonly device: GPUDevice },
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
        program = createWebGPUProgram(ref.device, name, source, defineValues, schema, onError);
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
