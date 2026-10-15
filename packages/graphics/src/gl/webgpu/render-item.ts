/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { idFactory } from '@molstar/core/util/id-factory';
import type {
  ProgramVariant,
  RenderItem,
  DrawMode,
  GlobalBindings,
  GraphicsProgram,
  MultiDrawBaseData,
  Transparency,
} from '../render-item.js';
import {
  type RenderableSchema,
  type RenderableValues,
  DefineSpec,
  getValueVersions,
  splitValues,
  type UniformSpec,
} from '../renderable/schema.js';
import type { TextureImage, TextureVolume } from '../renderable/util.js';
import type { DefineValues, ShaderCode } from '../shader-code.js';
import type { Texture, Textures, TextureValues } from '../texture.js';
import type { UniformKind, UniformsList } from '../uniform.js';
import {
  createUniformBlockData,
  createUniformBlockWriter,
  type UniformBlockLayout,
  getUniformBlockLayout,
} from '../uniform-layout.js';
import { BufferUsage } from './constants.js';
import type { WebGPUDrawContext } from './context.js';
import { asWebGPUGlobalBindings, getTextureEntries } from './global-bindings.js';
import type { VertexBufferSlot, WebGPUProgram } from './program.js';
import { getWgslProgram } from './shader/programs.js';
import { notSupported, WebGPUUnsupportedError } from './errors.js';
import { createWebGPUTexture, isLinearFiltered, isWebGPUTexture } from './texture.js';
import { createUniformSnapshot, type UniformSnapshot } from './uniform-ring.js';
import { WgslDrawGroup } from './uniform-declarations.js';
import { getWgslSamplerBinding, WgslSamplers } from './wgsl-declarations.js';

const getNextRenderItemId = idFactory();

function getRenderVariant(variant: ProgramVariant, transparency: Transparency | undefined): string {
  if (variant !== 'color' || !transparency) return variant;
  switch (transparency) {
    case 'blended':
      return 'colorBlended';
    case 'wboit':
      return 'colorWboit';
    case 'dpoit':
      return 'colorDpoit';
  }
}

type GPUData = { buffer: GPUBuffer; size: number };

/** Writes the bytes of `array`, a multiple of 4 */
function writeBytes(device: GPUDevice, buffer: GPUBuffer, array: ArrayBufferView) {
  const bytes = new Uint8Array(array.buffer as ArrayBuffer, array.byteOffset, Math.floor(array.byteLength / 4) * 4);
  device.queue.writeBuffer(buffer, 0, bytes);
}

function createGPUData(device: GPUDevice, array: ArrayBufferView, usage: number): GPUData {
  const size = Math.max(4, Math.ceil(array.byteLength / 4) * 4);
  const buffer = device.createBuffer({ size, usage: usage | BufferUsage.COPY_DST });
  if (array.byteLength > 0) writeBytes(device, buffer, array);
  return { buffer, size };
}

/** Rewrites `data` with `array` or recreates it when too small */
function updateGPUData(ctx: WebGPUDrawContext, data: GPUData, array: ArrayBufferView, usage: number): GPUData {
  const { device } = ctx;
  if (array.byteLength > data.size) {
    ctx.release(data.buffer);
    return createGPUData(device, array, usage);
  }
  if (array.byteLength > 0) writeBytes(device, data.buffer, array);
  return data;
}

/** Interleaves per-vertex `arrays` with `itemSizes` components each */
function interleave(arrays: readonly Float32Array[], itemSizes: readonly number[]) {
  let count = Infinity;
  let stride = 0;
  for (let i = 0; i < arrays.length; ++i) {
    count = Math.min(count, Math.floor(arrays[i].length / itemSizes[i]));
    stride += itemSizes[i];
  }
  const out = new Float32Array(count * stride);
  let offset = 0;
  for (let i = 0; i < arrays.length; ++i) {
    const array = arrays[i];
    const itemSize = itemSizes[i];
    for (let j = 0; j < count; ++j) {
      for (let k = 0; k < itemSize; ++k) out[j * stride + offset + k] = array[j * itemSize + k];
    }
    offset += itemSize;
  }
  return out;
}

function createTextures(ctx: WebGPUDrawContext, schema: RenderableSchema, values: TextureValues): Textures {
  const textures: Textures = [];
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'texture') continue;
    const value = values[k];
    if (!value) continue;
    if (spec.kind === 'texture') {
      textures.push([k, value.ref.value as Texture]);
    } else {
      const texture = createWebGPUTexture(ctx, spec.kind, spec.format, spec.dataType, spec.filter);
      texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
      textures.push([k, texture]);
    }
  }
  return textures;
}

type Block = { snapshot: UniformSnapshot; write: () => void; layout: { size: number } };

function createBlock(ctx: WebGPUDrawContext, layout: UniformBlockLayout, entries: UniformsList): Block {
  const snapshot = createUniformSnapshot(ctx.ring, layout.size);
  const data = createUniformBlockData(snapshot.bytes.buffer as ArrayBuffer, 0, layout.size);
  return { snapshot, write: createUniformBlockWriter(layout, data, entries, snapshot.markChanged), layout };
}

/**
 * Render item of graphics (`transparency` given) or compute renderables. The WGSL program is resolved at the
 * first render of a variant, so render items of unported programs can be created, rendering them fails.
 */
export function createWebGPURenderItem<T extends ProgramVariant>(
  ctx: WebGPUDrawContext,
  drawMode: DrawMode,
  shaderCode: ShaderCode,
  schema: RenderableSchema,
  values: RenderableValues,
  materialId: number,
  transparency: Transparency | undefined,
): RenderItem<T> {
  const id = getNextRenderItemId();
  const { stats } = ctx;

  if (schema.dRenderVariant === undefined) {
    Object.defineProperty(schema, 'dRenderVariant', { value: DefineSpec('string') });
  }

  const {
    attributeValues,
    defineValues,
    textureValues,
    materialTextureValues,
    uniformValues,
    materialUniformValues,
    bufferedUniformValues,
  } = splitValues(schema, values);

  const objectEntries: UniformsList = [];
  const materialEntries: UniformsList = [];
  const drawEntries: UniformsList = [];
  for (const e of [
    ...Object.entries(uniformValues),
    ...Object.entries(bufferedUniformValues),
    ...Object.entries(materialUniformValues),
  ]) {
    const { block } = schema[e[0]] as UniformSpec<UniformKind>;
    if (block === 'object') objectEntries.push(e);
    else if (block === 'material') materialEntries.push(e);
    else if (block === 'draw') drawEntries.push(e);
  }

  const objectLayout = getUniformBlockLayout('object', schema);
  const materialLayout = getUniformBlockLayout('material', schema);
  const objectBlock = objectLayout ? createBlock(ctx, objectLayout, objectEntries) : undefined;
  const materialBlock = materialLayout ? createBlock(ctx, materialLayout, materialEntries) : undefined;
  objectBlock?.write();
  materialBlock?.write();

  // 'draw' uniforms of the render item and of the global bindings, per global bindings instance; multi-draw
  // levels (`MultiDrawBaseData`) have their own block with their uniforms (e.g. `uLod`) replacing the item's
  type DrawBlocks = { base?: Block; levels: WeakMap<MultiDrawBaseData, Block> };
  const drawBlocks = new WeakMap<{ drawUniforms: UniformsList }, DrawBlocks>();
  const noGlobals = { drawUniforms: [] as UniformsList };
  function getDrawBlock(program: WebGPUProgram, globals: { drawUniforms: UniformsList }, level?: MultiDrawBaseData) {
    if (!program.drawLayout) return;
    let blocks = drawBlocks.get(globals);
    if (!blocks) drawBlocks.set(globals, (blocks = { levels: new WeakMap() }));
    const { members, size } = program.drawLayout;
    const layout: UniformBlockLayout = { block: 'object', members, size, hash: '' };
    if (!level) return (blocks.base ??= createBlock(ctx, layout, [...drawEntries, ...globals.drawUniforms]));
    let block = blocks.levels.get(level);
    if (!block) {
      const names = new Set(level.uniforms.map(([k]) => k));
      const entries = [...drawEntries.filter(([k]) => !names.has(k)), ...level.uniforms, ...globals.drawUniforms];
      block = createBlock(ctx, layout, entries);
      blocks.levels.set(level, block);
    }
    return block;
  }

  const versions = getValueVersions(values);

  /** per vertex buffer slot of the programs, created and updated on use */
  const vertexBuffers = new Map<string, GPUData & { versions: number[] }>();
  function getVertexBuffer(slot: VertexBufferSlot) {
    const cells = slot.names.map((name) => {
      const cell = attributeValues[name];
      if (!cell) throw new Error(`missing attribute '${name}' of '${shaderCode.name}'`);
      return cell;
    });
    let data = vertexBuffers.get(slot.key);
    if (data && cells.every((c, i) => c.ref.version === data!.versions[i])) return data.buffer;
    const arrays = cells.map((c) => c.ref.value as Float32Array);
    const array = arrays.length === 1 ? arrays[0] : interleave(arrays, slot.itemSizes);
    const gpuData = data
      ? updateGPUData(ctx, data, array, BufferUsage.VERTEX)
      : createGPUData(ctx.device, array, BufferUsage.VERTEX);
    data = { ...gpuData, versions: cells.map((c) => c.ref.version) };
    vertexBuffers.set(slot.key, data);
    return data.buffer;
  }
  let indexBuffer: GPUData | undefined;
  function createIndexBuffer() {
    if (values.elements?.ref.value) {
      indexBuffer = createGPUData(ctx.device, values.elements.ref.value, BufferUsage.INDEX);
    }
  }
  createIndexBuffer();

  const textures = createTextures(ctx, schema, textureValues);
  const materialTextures = createTextures(ctx, schema, materialTextureValues);

  let drawCount: number = values.drawCount.ref.value;
  let instanceCount: number = values.instanceCount.ref.value;
  stats.drawCount += drawCount;
  stats.instanceCount += instanceCount;
  stats.instancedDrawCount += instanceCount * drawCount;

  /**
   * Bound textures sampled with linear filtering, as the program key and as `dLinear_<texture>` defines (the filter is a property of the
   * texture as on WebGL, WGSL picks filtered or `textureLoad` reads by these defines), as a key and as define values
   */
  function getLinearTextures() {
    let key = '';
    for (const list of [textures, materialTextures]) {
      for (const [name, texture] of list) {
        if (isWebGPUTexture(texture) && isLinearFiltered(ctx.device, texture)) key += name + ',';
      }
    }
    return key;
  }

  function createProgram(variant: T, linearTextures: string) {
    const wgsl = getWgslProgram(shaderCode.name);
    if (!wgsl) notSupported(`the '${shaderCode.name}' program`);
    if (!wgsl.variants.includes(variant)) {
      notSupported(`render variant '${variant}' of '${shaderCode.name}'`);
    }
    const variantDefines: DefineValues = {
      ...(defineValues as DefineValues),
      ...Object.fromEntries(
        linearTextures
          .split(',')
          .filter((n) => n)
          .map((n) => [`dLinear_${n}`, ValueCell.create(true)]),
      ),
      dRenderVariant: ValueCell.create(getRenderVariant(variant, transparency)),
    };
    const unsupported = wgsl.getUnsupported?.(variantDefines);
    if (unsupported) notSupported(`${unsupported} of '${shaderCode.name}'`);
    return ctx.programs.get(shaderCode.name, wgsl.source, variantDefines, schema);
  }

  /** `null` if the program is not ported and the context is not strict, the render item then draws nothing */
  const programs = new Map<string, WebGPUProgram | null>();
  function getWebGPUProgram(variant: T) {
    const linearTextures = getLinearTextures();
    const programKey = `${variant}|${linearTextures}`;
    let program = programs.get(programKey);
    if (program === undefined) {
      try {
        program = createProgram(variant, linearTextures);
      } catch (e) {
        if (ctx.strict || !(e instanceof WebGPUUnsupportedError)) throw e;
        ctx.warnUnsupported(e.message);
        program = null;
      }
      programs.set(programKey, program);
    }
    return program;
  }

  // bind groups 1 (material) and 2 (object) per bind group layout
  type CachedGroup = { group: GPUBindGroup; buffer: GPUBuffer | undefined; key: string };
  const bindGroups = new Map<GPUBindGroupLayout, CachedGroup>();
  function getBindGroup(program: WebGPUProgram, index: 1 | 2, block: Block | undefined, tierTextures: Textures) {
    const tier = program.bindingLayout[index === 1 ? 'material' : 'object'];
    const slot = tier.uniforms && block ? block.snapshot.getSlot() : undefined;
    const tex = getTextureEntries(ctx.device, tierTextures, tier.textures, program.referenced);
    const layout = program.bindGroupLayouts[index];
    let cached = bindGroups.get(layout);
    if (!cached || cached.buffer !== slot?.buffer || cached.key !== tex.key) {
      const entries: GPUBindGroupEntry[] = [...tex.entries];
      if (slot) entries.unshift({ binding: 0, resource: { buffer: slot.buffer, offset: 0, size: block!.layout.size } });
      cached = { group: ctx.device.createBindGroup({ layout, entries }), buffer: slot?.buffer, key: tex.key };
      bindGroups.set(layout, cached);
    }
    return { group: cached.group, offsets: slot ? [slot.offset] : [] };
  }

  // group 0 without global bindings (compute renderables): only the shared samplers
  const samplerGroups = new Map<GPUBindGroupLayout, GPUBindGroup>();
  function getSamplerBindGroup(program: WebGPUProgram) {
    if (program.bindingLayout.global.uniforms) {
      throw new Error(`'${shaderCode.name}' uses the global block, render it with global bindings`);
    }
    const layout = program.bindGroupLayouts[0];
    let group = samplerGroups.get(layout);
    if (!group) {
      const samplerBinding = getWgslSamplerBinding(program.bindingLayout);
      const entries: GPUBindGroupEntry[] = [];
      for (let i = 0; i < WgslSamplers.length; ++i) {
        const name = WgslSamplers[i];
        if (program.referenced.has(name)) entries.push({ binding: samplerBinding + i, resource: ctx.samplers[name] });
      }
      group = ctx.device.createBindGroup({ layout, entries });
      samplerGroups.set(layout, group);
    }
    return group;
  }

  const emptyDrawGroups = new Map<GPUBindGroupLayout, GPUBindGroup>();
  function getDrawBindGroup(
    program: WebGPUProgram,
    globals: { drawUniforms: UniformsList },
    level?: MultiDrawBaseData,
  ) {
    const layout = program.bindGroupLayouts[WgslDrawGroup];
    const block = getDrawBlock(program, globals, level);
    if (!block) {
      let group = emptyDrawGroups.get(layout);
      if (!group) emptyDrawGroups.set(layout, (group = ctx.device.createBindGroup({ layout, entries: [] })));
      return { group, offsets: [] as number[] };
    }
    block.write();
    const slot = block.snapshot.getSlot();
    let cached = bindGroups.get(layout);
    if (!cached || cached.buffer !== slot.buffer) {
      const entries: GPUBindGroupEntry[] = [
        { binding: 0, resource: { buffer: slot.buffer, offset: 0, size: block.layout.size } },
      ];
      cached = { group: ctx.device.createBindGroup({ layout, entries }), buffer: slot.buffer, key: '' };
      bindGroups.set(layout, cached);
    }
    return { group: cached.group, offsets: [slot.offset] };
  }

  const graphicsPrograms = new Map<T, GraphicsProgram>();
  function getGraphicsProgram(variant: T): GraphicsProgram {
    let p = graphicsPrograms.get(variant);
    if (!p) {
      p = { id: id * 16 + graphicsPrograms.size, variant, isReady: () => true, link: () => {}, finalize: () => true };
      graphicsPrograms.set(variant, p);
    }
    return p;
  }

  let destroyed = false;

  /** after the context restored a lost device: GPU objects are created anew on use, data uploaded again */
  const restorable = {
    reset: () => {
      vertexBuffers.clear();
      createIndexBuffer();
      programs.clear();
      bindGroups.clear();
      samplerGroups.clear();
      emptyDrawGroups.clear();
    },
  };
  ctx.restorables.add(restorable);

  return {
    id,
    materialId,
    getByteCount: () => {
      let count = 0;
      for (const b of vertexBuffers.values()) count += b.size;
      if (indexBuffer) count += indexBuffer.size;
      for (const [, t] of textures) count += t.getByteCount();
      for (const [, t] of materialTextures) count += t.getByteCount();
      return count;
    },
    getProgram: getGraphicsProgram,
    setTransparency: (value: Transparency) => {
      if (value === transparency) return;
      transparency = value;
      programs.clear();
    },

    render: (variant: T, globals?: GlobalBindings, mdbDataList?: MultiDrawBaseData[]) => {
      if (drawCount === 0 || instanceCount === 0) return;
      const gpuGlobals = globals ? asWebGPUGlobalBindings(globals) : undefined;
      const program = getWebGPUProgram(variant);
      if (!program) return;

      const pipeline = ctx.getPipeline(program, drawMode);
      // compiling in the background, drawn once ready
      if (!pipeline) return;
      const pass = ctx.getPassEncoder();
      pass.setPipeline(pipeline);

      if (gpuGlobals) {
        const global = gpuGlobals.getBindGroup(program);
        pass.setBindGroup(0, global.group, [global.offset]);
      } else {
        pass.setBindGroup(0, getSamplerBindGroup(program));
      }
      materialBlock?.write();
      const material = getBindGroup(program, 1, materialBlock, materialTextures);
      pass.setBindGroup(1, material.group, material.offsets);
      objectBlock?.write();
      const object = getBindGroup(program, 2, objectBlock, textures);
      pass.setBindGroup(2, object.group, object.offsets);
      const slots = program.vertexBufferSlots;
      for (let i = 0, il = slots.length; i < il; ++i) {
        pass.setVertexBuffer(i, getVertexBuffer(slots[i]));
      }
      if (indexBuffer) pass.setIndexBuffer(indexBuffer.buffer, 'uint32');

      if (mdbDataList) {
        // culled and level-of-detail ranges, `firstInstance` offsets the instanced attributes
        for (const level of mdbDataList) {
          if (level.count === 0) continue;
          const draw = getDrawBindGroup(program, gpuGlobals ?? noGlobals, level);
          pass.setBindGroup(WgslDrawGroup, draw.group, draw.offsets);
          const { counts, instanceCounts, firsts, offsets, baseVertices, baseInstances } = level;
          for (let i = 0; i < level.count; ++i) {
            if (counts[i] === 0 || instanceCounts[i] === 0) continue;
            if (indexBuffer) {
              // element offsets are in bytes
              pass.drawIndexed(counts[i], instanceCounts[i], offsets[i] / 4, baseVertices[i] ?? 0, baseInstances[i]);
            } else {
              pass.draw(counts[i], instanceCounts[i], firsts[i], baseInstances[i]);
            }
            stats.calls.drawInstancedBase += 1;
            stats.calls.counts += instanceCounts[i];
          }
        }
      } else {
        const draw = getDrawBindGroup(program, gpuGlobals ?? noGlobals);
        pass.setBindGroup(WgslDrawGroup, draw.group, draw.offsets);
        if (indexBuffer) pass.drawIndexed(drawCount, instanceCount);
        else pass.draw(drawCount, instanceCount);
        stats.calls.drawInstanced += 1;
        stats.calls.counts += instanceCount;
      }
    },

    update: () => {
      let definesChanged = false;
      for (const k of Object.keys(defineValues)) {
        const v = defineValues[k];
        if (v.ref.version !== versions[k]) {
          definesChanged = true;
          versions[k] = v.ref.version;
        }
      }
      if (definesChanged) programs.clear();

      if (values.drawCount.ref.version !== versions.drawCount) {
        stats.drawCount += values.drawCount.ref.value - drawCount;
        drawCount = values.drawCount.ref.value;
        versions.drawCount = values.drawCount.ref.version;
      }
      if (values.instanceCount.ref.version !== versions.instanceCount) {
        stats.instanceCount += values.instanceCount.ref.value - instanceCount;
        instanceCount = values.instanceCount.ref.value;
        versions.instanceCount = values.instanceCount.ref.version;
      }

      if (indexBuffer && values.elements.ref.version !== versions.elements) {
        indexBuffer = updateGPUData(ctx, indexBuffer, values.elements.ref.value, BufferUsage.INDEX);
        versions.elements = values.elements.ref.version;
      }

      const updateTextures = (list: Textures, textureValues: TextureValues) => {
        for (let i = 0, il = list.length; i < il; ++i) {
          const [k, texture] = list[i];
          const value = textureValues[k];
          if (value.ref.version === versions[k]) continue;
          // textures of kind 'texture' are managed externally
          if (schema[k].kind !== 'texture') {
            texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
          } else {
            list[i][1] = value.ref.value as Texture;
          }
          versions[k] = value.ref.version;
        }
      };
      updateTextures(textures, textureValues);
      updateTextures(materialTextures, materialTextureValues);
    },

    destroy: () => {
      if (destroyed) return;
      for (const b of vertexBuffers.values()) ctx.release(b.buffer);
      if (indexBuffer) ctx.release(indexBuffer.buffer);
      for (const [k, t] of textures) if (schema[k].kind !== 'texture') t.destroy();
      for (const [k, t] of materialTextures) if (schema[k].kind !== 'texture') t.destroy();
      stats.drawCount -= drawCount;
      stats.instanceCount -= instanceCount;
      stats.instancedDrawCount -= instanceCount * drawCount;
      ctx.restorables.delete(restorable);
      destroyed = true;
    },
  };
}
