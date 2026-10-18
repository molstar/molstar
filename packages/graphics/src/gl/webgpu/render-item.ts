/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ReversedZeroToOne } from '../depth-range.js';
import { ValueCell } from '@molstar/core/util';
import { idFactory } from '@molstar/core/util/id-factory';
import type {
  ProgramVariant,
  RenderItem,
  DrawMode,
  GlobalBindings,
  GraphicsProgram,
  MultiDrawBaseData,
  ProgramItem,
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
  type UniformBlockWriter,
  getUniformBlockLayout,
} from '../uniform-layout.js';
import { BufferUsage } from './constants.js';
import { asWebGPUStorageBuffer } from './compute.js';
import type { DrawIndirect } from '../compute.js';
import type { WebGPUDrawContext, WebGPULinkable } from './context.js';
import { asWebGPUGlobalBindings, getTextureEntries } from './global-bindings.js';
import type { VertexBufferSlot, WebGPUProgram } from './program.js';
import { getWgslProgram } from './shader/programs.js';
import { notSupported, WebGPUUnsupportedError } from './errors.js';
import {
  asWebGPUTexture,
  createWebGPUTexture,
  isLinearFiltered,
  isWebGPUTexture,
  type WebGPUTexture,
} from './texture.js';
import { createUniformArena, type UniformArena, type UniformSlot, type UniformSnapshot } from './uniform-ring.js';
import { getWgslSamplerBinding, WgslSamplers } from './wgsl-declarations.js';

const getNextRenderItemId = idFactory();

/**
 * The program of `variant` of a render item of `shaderCode` and `schema` with `defineValues`, its `linearTextures` read
 * filtered (`dLinear_*`), drawn with `transparency`: shared by render items and program items (`createWebGPUProgramItem`)
 */
function createVariantProgram(
  ctx: WebGPUDrawContext,
  shaderCode: ShaderCode,
  schema: RenderableSchema,
  defineValues: DefineValues,
  variant: ProgramVariant,
  transparency: Transparency | undefined,
  linearTextures: readonly string[],
) {
  const wgsl = getWgslProgram(shaderCode.name);
  if (!wgsl) notSupported(`the '${shaderCode.name}' program`);
  if (!wgsl.variants.includes(variant)) {
    notSupported(`render variant '${variant}' of '${shaderCode.name}'`);
  }
  const variantDefines: DefineValues = {
    ...defineValues,
    ...Object.fromEntries(linearTextures.map((n) => [`dLinear_${n}`, ValueCell.create(true)])),
    dRenderVariant: ValueCell.create(getRenderVariant(variant, transparency)),
  };
  const unsupported = wgsl.getUnsupported?.(variantDefines);
  if (unsupported) notSupported(`${unsupported} of '${shaderCode.name}'`);
  return ctx.programs.get(shaderCode.name, wgsl.source, variantDefines, schema);
}

const getNextProgramItemId = idFactory();

/**
 * Programs of a render item without its data, see `GraphicsContext.createProgramItem`: `resources.preparePipelines`
 * creates them (from the same cache as render items) and their pipelines
 */
export function createWebGPUProgramItem(
  ctx: WebGPUDrawContext,
  drawMode: DrawMode,
  shaderCode: ShaderCode,
  schema: RenderableSchema,
  defineValues: DefineValues,
  linearTextures: readonly string[],
  transparency: Transparency,
): ProgramItem {
  const id = getNextProgramItemId();
  if (schema.dRenderVariant === undefined) {
    Object.defineProperty(schema, 'dRenderVariant', { value: DefineSpec('string') });
  }
  const programs = new Map<ProgramVariant, WebGPUProgram | null>();
  const getProgram = (variant: ProgramVariant) => {
    let program = programs.get(variant);
    if (program === undefined) {
      try {
        program = createVariantProgram(ctx, shaderCode, schema, defineValues, variant, transparency, linearTextures);
      } catch {
        // not ported or not supported: a render item of it reports that when drawn
        program = null;
      }
      programs.set(variant, program);
    }
    return program;
  };
  const item: ProgramItem = {
    id,
    getProgram: (variant) => ({ id: id * 16, variant, isReady: () => true, link: () => {}, finalize: () => true }),
    destroy: () => programs.clear(),
  };
  ctx.itemPreparers.set(item, (variant, state, target) => {
    const program = getProgram(variant);
    if (program) ctx.preparePipeline(program, drawMode, state, target);
  });
  return item;
}

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

type GPUData = { buffer: GPUBuffer; size: number; usage: number };

/** Writes the bytes of `array`, a multiple of 4 */
function writeBytes(device: GPUDevice, buffer: GPUBuffer, array: ArrayBufferView) {
  const bytes = new Uint8Array(array.buffer as ArrayBuffer, array.byteOffset, Math.floor(array.byteLength / 4) * 4);
  device.queue.writeBuffer(buffer, 0, bytes);
}

/** Vertex buffers count as attributes, index buffers as elements (`GraphicsResources.getByteCounts`) */
const byteKind = (usage: number) => (usage & BufferUsage.INDEX ? 'elements' : 'attribute');

function createGPUData(ctx: WebGPUDrawContext, array: ArrayBufferView, usage: number): GPUData {
  const { device } = ctx;
  const size = Math.max(4, Math.ceil(array.byteLength / 4) * 4);
  const buffer = device.createBuffer({ size, usage: usage | BufferUsage.COPY_DST });
  if (array.byteLength > 0) writeBytes(device, buffer, array);
  const kind = byteKind(usage);
  ctx.byteCounts[kind] += size;
  ctx.stats.resourceCounts[kind] += 1;
  return { buffer, size, usage };
}

function releaseGPUData(ctx: WebGPUDrawContext, data: GPUData) {
  ctx.release(data.buffer);
  const kind = byteKind(data.usage);
  ctx.byteCounts[kind] -= data.size;
  ctx.stats.resourceCounts[kind] -= 1;
}

/** Rewrites `data` with `array` or recreates it when too small */
function updateGPUData(ctx: WebGPUDrawContext, data: GPUData, array: ArrayBufferView, usage: number): GPUData {
  const { device } = ctx;
  if (array.byteLength > data.size) {
    releaseGPUData(ctx, data);
    return createGPUData(ctx, array, usage);
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

type Block = { snapshot: UniformSnapshot; write: UniformBlockWriter; layout: { size: number } };

function createBlock(snapshot: UniformSnapshot, layout: UniformBlockLayout, entries: UniformsList): Block {
  const data = createUniformBlockData(snapshot.bytes.buffer as ArrayBuffer, 0, layout.size);
  return { snapshot, write: createUniformBlockWriter(layout, data, entries, snapshot.markChanged), layout };
}

/** A bind group with what it was created for, see `isCurrent` */
type CachedGroup = {
  group: GPUBindGroup;
  buffer: GPUBuffer | undefined;
  drawBuffer: GPUBuffer | undefined;
  bindingsVersion: number;
  textures: (WebGPUTexture | undefined)[];
  views: (GPUTextureView | undefined)[];
  /** the bound storage buffers, see `WebGPUProgram.storageBindings` */
  storageBuffers: GPUBuffer[];
  offsets: number[];
};

/**
 * The material block of the render items of a material id, shared as the WebGL2 material buffers are: material
 * uniforms are the same for all items of a material, each item writes its values into the shared data; with the
 * bind groups of the block (materials have no textures)
 */
type SharedMaterial = {
  readonly key: string;
  readonly arena: UniformArena;
  readonly snapshot: UniformSnapshot;
  readonly groups: Map<GPUBindGroupLayout, CachedGroup>;
  refCount: number;
};
const SharedMaterials = new WeakMap<WebGPUDrawContext, Map<string, SharedMaterial>>();

function acquireMaterial(ctx: WebGPUDrawContext, materialId: number, layout: UniformBlockLayout): SharedMaterial {
  let materials = SharedMaterials.get(ctx);
  if (!materials) SharedMaterials.set(ctx, (materials = new Map()));
  // -1: a material of its own
  const key = materialId === -1 ? '' : `${materialId}|${layout.hash}`;
  let material = key ? materials.get(key) : undefined;
  if (!material) {
    const arena = createUniformArena(ctx.ring, ctx.release);
    material = { key, arena, snapshot: arena.createSnapshot(layout.size), groups: new Map(), refCount: 0 };
    if (key) materials.set(key, material);
  }
  material.refCount += 1;
  return material;
}

function releaseMaterial(ctx: WebGPUDrawContext, material: SharedMaterial) {
  material.refCount -= 1;
  if (material.refCount > 0) return;
  material.arena.destroy();
  if (material.key) SharedMaterials.get(ctx)?.delete(material.key);
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
    storageBufferValues,
    // the projections mapped to reversed `[0, 1]` depth; window rectangles as given (top-left origin)
  } = splitValues(schema, values, { depthRange: ReversedZeroToOne });

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

  /** the blocks of the item, their writers stop observing the cells when it is destroyed */
  const itemBlocks: Block[] = [];
  const ownBlock = (snapshot: UniformSnapshot, layout: UniformBlockLayout, entries: UniformsList) => {
    const block = createBlock(snapshot, layout, entries);
    itemBlocks.push(block);
    return block;
  };
  const objectLayout = getUniformBlockLayout('object', schema);
  const materialLayout = getUniformBlockLayout('material', schema);
  /** the object block and the draw blocks, in one buffer: one bind group (2) with dynamic offsets */
  const arena = createUniformArena(ctx.ring, ctx.release);
  const objectBlock = objectLayout
    ? ownBlock(arena.createSnapshot(objectLayout.size), objectLayout, objectEntries)
    : undefined;
  const material = materialLayout ? acquireMaterial(ctx, materialId, materialLayout) : undefined;
  const materialBlock = material ? ownBlock(material.snapshot, materialLayout!, materialEntries) : undefined;
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
    if (!level) {
      return (blocks.base ??= ownBlock(arena.createSnapshot(size), layout, [...drawEntries, ...globals.drawUniforms]));
    }
    let block = blocks.levels.get(level);
    if (!block) {
      // levels are kept by the renderable, replaced when their capacity grows (the slot of a replaced one stays unused)
      const names = new Set(level.uniforms.map(([k]) => k));
      const entries = [...drawEntries.filter(([k]) => !names.has(k)), ...level.uniforms, ...globals.drawUniforms];
      block = ownBlock(arena.createSnapshot(size), layout, entries);
      blocks.levels.set(level, block);
    }
    return block;
  }

  const versions = getValueVersions(values);

  /** per vertex buffer slot of the programs, created and updated on use; `data` is unset until uploaded */
  const vertexBuffers = new Map<string, { cells: ValueCell<any>[]; versions: number[]; data?: GPUData }>();
  const slotEntries = new Map<VertexBufferSlot, NonNullable<ReturnType<typeof vertexBuffers.get>>>();
  function getVertexBuffer(slot: VertexBufferSlot) {
    // by slot object first (stable per program), by key across programs
    let entry = slotEntries.get(slot) ?? vertexBuffers.get(slot.key);
    if (!entry) {
      const cells = slot.names.map((name) => {
        const cell = attributeValues[name];
        if (!cell) throw new Error(`missing attribute '${name}' of '${shaderCode.name}'`);
        return cell;
      });
      entry = { cells, versions: cells.map(() => -1) };
      vertexBuffers.set(slot.key, entry);
    }
    slotEntries.set(slot, entry);
    const { cells, versions } = entry;
    let changed = !entry.data;
    for (let i = 0, il = cells.length; i < il; ++i) {
      if (cells[i].ref.version !== versions[i]) changed = true;
    }
    if (!changed) return entry.data!.buffer;
    const arrays = cells.map((c) => c.ref.value as Float32Array);
    const array = arrays.length === 1 ? arrays[0] : interleave(arrays, slot.itemSizes);
    entry.data = entry.data
      ? updateGPUData(ctx, entry.data, array, BufferUsage.VERTEX)
      : createGPUData(ctx, array, BufferUsage.VERTEX);
    for (let i = 0, il = cells.length; i < il; ++i) versions[i] = cells[i].ref.version;
    return entry.data.buffer;
  }
  let indexBuffer: GPUData | undefined;
  function createIndexBuffer() {
    if (values.elements?.ref.value) {
      indexBuffer = createGPUData(ctx, values.elements.ref.value, BufferUsage.INDEX);
    }
  }
  createIndexBuffer();

  const textures = createTextures(ctx, schema, textureValues);
  const materialTextures = createTextures(ctx, schema, materialTextureValues);

  let drawCount: number = values.drawCount.ref.value;
  let instanceCount: number = values.instanceCount.ref.value;
  /** the draw arguments and instance count last written into them */
  let indirectWritten: { args: DrawIndirect; instanceCount: number } | undefined;
  function getDrawIndirect() {
    const args = (values.drawIndirect as ValueCell<DrawIndirect | undefined> | undefined)?.ref.value;
    if (!args) return undefined;
    if (indexBuffer) notSupported('indirect draws of indexed geometry');
    if (indirectWritten?.args !== args || indirectWritten.instanceCount !== instanceCount) {
      args.buffer.write(new Uint32Array([instanceCount]), args.byteOffset + 4);
      indirectWritten = { args, instanceCount };
    }
    return args;
  }
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
    const names = linearTextures.split(',').filter((n) => n);
    return createVariantProgram(ctx, shaderCode, schema, defineValues as DefineValues, variant, transparency, names);
  }

  /** `null` if the program is not ported and the context is not strict, the render item then draws nothing */
  const programs = new Map<string, WebGPUProgram | null>();
  /** the program of each variant for the current textures, see `invalidatePrograms` */
  const variantPrograms = new Map<T, WebGPUProgram | null>();
  function getWebGPUProgram(variant: T) {
    let program = variantPrograms.get(variant);
    if (program !== undefined) return program;
    const linearTextures = getLinearTextures();
    const programKey = `${variant}|${linearTextures}`;
    program = programs.get(programKey);
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
    variantPrograms.set(variant, program);
    return program;
  }

  /**
   * After a change of what selects the programs: the defines, the transparency, the bound textures (their filter selects
   * linear reads) or the device; `clearAll` also drops the programs of other texture configurations
   */
  function invalidatePrograms(clearAll: boolean) {
    variantPrograms.clear();
    if (clearAll) programs.clear();
    textureBindingsVersion += 1;
    ctx.linkables.add(linkable);
  }
  /** incremented when the bound textures or the programs change, bind groups resolve their textures again */
  let textureBindingsVersion = 0;

  // bind groups 1 (material) and 2 (object, with the draw uniforms) per bind group layout: the bound textures are
  // resolved when the textures or programs change, per draw only their views (recreated textures) and the uniform
  // buffers are compared
  const bindGroups = new Map<GPUBindGroupLayout, CachedGroup>();
  // materials have no textures: their bind groups are shared with the other items of the material
  const materialGroups = material && materialTextures.length === 0 ? material.groups : bindGroups;
  const noOffsets: number[] = [];
  const noStorageBuffers: GPUBuffer[] = [];
  const globalOffsets = [0];
  function isCurrent(
    cached: CachedGroup,
    buffer: GPUBuffer | undefined,
    drawBuffer: GPUBuffer | undefined,
    bindingsVersion: number,
  ) {
    if (cached.buffer !== buffer || cached.drawBuffer !== drawBuffer || cached.bindingsVersion !== bindingsVersion) {
      return false;
    }
    const { textures, views } = cached;
    for (let i = 0, il = textures.length; i < il; ++i) {
      if (textures[i]?.getView() !== views[i]) return false;
    }
    const { storageBuffers } = cached;
    for (let i = 0, il = storageBuffers.length; i < il; ++i) {
      if (storageBuffers[i] !== currentStorageBuffers[i]) return false;
    }
    return true;
  }
  /** the storage buffers of the program being bound, in the order of its `storageBindings` */
  const currentStorageBuffers: GPUBuffer[] = [];
  function getStorageBuffers(program: WebGPUProgram) {
    currentStorageBuffers.length = 0;
    for (const { name } of program.storageBindings) {
      const buffer = storageBufferValues[name]?.ref.value;
      if (!buffer) throw new Error(`missing storage buffer '${name}' of '${shaderCode.name}'`);
      currentStorageBuffers.push(asWebGPUStorageBuffer(buffer).getGPUBuffer());
    }
    return currentStorageBuffers;
  }
  function getBindGroup(
    program: WebGPUProgram,
    index: 1 | 2,
    block: Block | undefined,
    tierTextures: Textures,
    drawBlock?: Block,
  ) {
    const tier = program.bindingLayout[index === 1 ? 'material' : 'object'];
    const slot = tier.uniforms && block ? block.snapshot.getSlot() : undefined;
    let drawSlot: UniformSlot | undefined;
    if (drawBlock) {
      drawBlock.write();
      drawSlot = drawBlock.snapshot.getSlot();
    }
    const layout = program.bindGroupLayouts[index];
    const groups = index === 1 ? materialGroups : bindGroups;
    const bindingsVersion = groups === bindGroups ? textureBindingsVersion : 0;
    const storageBuffers = index === 2 ? getStorageBuffers(program) : noStorageBuffers;
    let cached = groups.get(layout);
    if (!cached || !isCurrent(cached, slot?.buffer, drawSlot?.buffer, bindingsVersion)) {
      const tex = getTextureEntries(ctx.device, tierTextures, tier.textures, program.referenced);
      const entries: GPUBindGroupEntry[] = [...tex.entries];
      if (slot) entries.unshift({ binding: 0, resource: { buffer: slot.buffer, offset: 0, size: block!.layout.size } });
      if (drawSlot) {
        entries.push({
          binding: program.drawBinding,
          resource: { buffer: drawSlot.buffer, offset: 0, size: drawBlock!.layout.size },
        });
      }
      if (index === 2) {
        for (let i = 0; i < program.storageBindings.length; ++i) {
          entries.push({ binding: program.storageBindings[i].binding, resource: { buffer: storageBuffers[i] } });
        }
      }
      const textures: (WebGPUTexture | undefined)[] = [];
      for (const { name } of tier.textures) {
        if (!program.referenced.has(name)) continue;
        const t = tierTextures.find((e) => e[0] === name)?.[1];
        textures.push(t ? asWebGPUTexture(t) : undefined);
      }
      cached = {
        group: ctx.device.createBindGroup({ layout, entries }),
        buffer: slot?.buffer,
        drawBuffer: drawSlot?.buffer,
        bindingsVersion,
        textures,
        views: textures.map((t) => t?.getView()),
        storageBuffers: storageBuffers.slice(),
        offsets: [],
      };
      groups.set(layout, cached);
    }
    // dynamic offsets in binding order: the uniform block (0), then the draw uniforms
    const { offsets } = cached;
    offsets.length = 0;
    if (slot) offsets.push(slot.offset);
    if (drawSlot) offsets.push(drawSlot.offset);
    return { group: cached.group, offsets: offsets.length ? offsets : noOffsets };
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

  // created ahead of the first draw by `linkPrograms`, which starts compiling the pipelines of known keys
  const linkable: WebGPULinkable = {
    link: (variant) => {
      if (!getWgslProgram(shaderCode.name)?.variants.includes(variant)) return;
      try {
        getWebGPUProgram(variant as T);
      } catch {
        // reported by the first draw of the variant
      }
    },
  };
  ctx.linkables.add(linkable);

  const item: RenderItem<T> = {
    id,
    materialId,
    getLinearTextures: () =>
      getLinearTextures()
        .split(',')
        .filter((n) => n),
    getByteCount: () => {
      let count = 0;
      for (const b of vertexBuffers.values()) count += b.data?.size ?? 0;
      if (indexBuffer) count += indexBuffer.size;
      for (const [, t] of textures) count += t.getByteCount();
      for (const [, t] of materialTextures) count += t.getByteCount();
      return count;
    },
    getProgram: getGraphicsProgram,
    setTransparency: (value: Transparency) => {
      if (value === transparency) return;
      transparency = value;
      invalidatePrograms(true);
    },

    render: (variant: T, globals?: GlobalBindings, mdbDataList?: MultiDrawBaseData[]) => {
      if (drawCount === 0 || instanceCount === 0) return;
      const gpuGlobals = globals ? asWebGPUGlobalBindings(globals) : undefined;
      const program = getWebGPUProgram(variant);
      if (!program) return;

      const pipeline = ctx.getPipeline(program, drawMode);
      // compiling in the background, drawn once ready
      if (!pipeline) return;
      const state = ctx.getPassState();
      const pass = state.encoder;
      state.setPipeline(pipeline);

      if (gpuGlobals) {
        const global = gpuGlobals.getBindGroup(program);
        globalOffsets[0] = global.offset;
        state.setBindGroup(0, global.group, globalOffsets);
      } else {
        state.setBindGroup(0, getSamplerBindGroup(program));
      }
      materialBlock?.write();
      const material = getBindGroup(program, 1, materialBlock, materialTextures);
      state.setBindGroup(1, material.group, material.offsets);
      objectBlock?.write();
      const slots = program.vertexBufferSlots;
      for (let i = 0, il = slots.length; i < il; ++i) {
        state.setVertexBuffer(i, getVertexBuffer(slots[i]));
      }
      if (indexBuffer) state.setIndexBuffer(indexBuffer.buffer, 'uint32');

      if (mdbDataList) {
        // culled and level-of-detail ranges, `firstInstance` offsets the instanced attributes
        for (const level of mdbDataList) {
          if (level.count === 0) continue;
          const drawBlock = getDrawBlock(program, gpuGlobals ?? noGlobals, level);
          const object = getBindGroup(program, 2, objectBlock, textures, drawBlock);
          state.setBindGroup(2, object.group, object.offsets);
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
        const drawBlock = getDrawBlock(program, gpuGlobals ?? noGlobals);
        const object = getBindGroup(program, 2, objectBlock, textures, drawBlock);
        state.setBindGroup(2, object.group, object.offsets);
        const indirect = getDrawIndirect();
        if (indirect) pass.drawIndirect(asWebGPUStorageBuffer(indirect.buffer).getGPUBuffer(), indirect.byteOffset);
        else if (indexBuffer) pass.drawIndexed(drawCount, instanceCount);
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
      if (definesChanged) invalidatePrograms(true);

      if (values.drawCount.ref.version !== versions.drawCount) {
        stats.drawCount += values.drawCount.ref.value - drawCount;
        stats.instancedDrawCount += instanceCount * (values.drawCount.ref.value - drawCount);
        drawCount = values.drawCount.ref.value;
        versions.drawCount = values.drawCount.ref.version;
      }
      if (values.instanceCount.ref.version !== versions.instanceCount) {
        stats.instanceCount += values.instanceCount.ref.value - instanceCount;
        stats.instancedDrawCount += (values.instanceCount.ref.value - instanceCount) * drawCount;
        instanceCount = values.instanceCount.ref.value;
        versions.instanceCount = values.instanceCount.ref.version;
      }

      if (indexBuffer && values.elements.ref.version !== versions.elements) {
        indexBuffer = updateGPUData(ctx, indexBuffer, values.elements.ref.value, BufferUsage.INDEX);
        versions.elements = values.elements.ref.version;
      }

      let texturesReplaced = false;
      const updateTextures = (list: Textures, textureValues: TextureValues) => {
        for (let i = 0, il = list.length; i < il; ++i) {
          const [k, texture] = list[i];
          const value = textureValues[k];
          if (value.ref.version === versions[k]) continue;
          // textures of kind 'texture' are managed externally
          if (schema[k].kind !== 'texture') {
            texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
          } else if (list[i][1] !== value.ref.value) {
            list[i][1] = value.ref.value as Texture;
            texturesReplaced = true;
          }
          versions[k] = value.ref.version;
        }
      };
      updateTextures(textures, textureValues);
      updateTextures(materialTextures, materialTextureValues);
      // another texture may be filtered differently (a program of other linear reads)
      if (texturesReplaced && !definesChanged) invalidatePrograms(false);
    },

    destroy: () => {
      if (destroyed) return;
      for (const b of vertexBuffers.values()) if (b.data) releaseGPUData(ctx, b.data);
      if (indexBuffer) releaseGPUData(ctx, indexBuffer);
      for (const [k, t] of textures) if (schema[k].kind !== 'texture') t.destroy();
      for (const [k, t] of materialTextures) if (schema[k].kind !== 'texture') t.destroy();
      stats.drawCount -= drawCount;
      stats.instanceCount -= instanceCount;
      stats.instancedDrawCount -= instanceCount * drawCount;
      for (const b of itemBlocks) b.write.dispose();
      arena.destroy();
      if (material) releaseMaterial(ctx, material);
      ctx.linkables.delete(linkable);
      destroyed = true;
    },
  };
  // for `resources.preparePipelines`
  ctx.itemPreparers.set(item as RenderItem<any>, (variant, state, target) => {
    let program: WebGPUProgram | null;
    try {
      program = getWebGPUProgram(variant as T);
    } catch {
      return; // reported by the first draw of the variant
    }
    if (program) ctx.preparePipeline(program, drawMode, state, target);
  });
  return item;
}
