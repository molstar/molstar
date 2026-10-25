/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isQuadIndexing, QuadIndicesPerItem } from '../quad-indices.js';
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
import { IndirectRecordBytes } from '../render-item.js';
import {
  type RenderableSchema,
  type RenderableValues,
  DefineSpec,
  splitValues,
  type TextureSpec,
  type UniformSpec,
} from '../renderable/schema.js';
import type { DataArray, TextureImage, TextureVolume } from '../renderable/util.js';
import type { DefineType, DefineValues, ShaderCode } from '../shader-code.js';
import type { Texture, TextureKind, Textures, TextureValues } from '../texture.js';
import { type SharedTexture, type SharedTextureImage, updateSharedTexture } from '../shared-texture.js';
import type { UniformKind, UniformsList } from '../uniform.js';
import {
  createUniformBlockData,
  createUniformBlockWriter,
  type UniformBlockLayout,
  type UniformBlockWriter,
  getUniformBlockLayout,
} from '../uniform-layout.js';
import { BufferUsage } from './constants.js';
import type { BufferRange } from './buffer-pages.js';
import { SegmentedArray } from '../segmented-array.js';
import {
  type DataLayer,
  DataLayerNames,
  DataLayerElementBytes,
  DataLayers,
  DataLayerPagesName,
  DataPageNames,
  getDataLayerPages,
  getDataLayerSchema,
  getDataOffsetName,
  MaxDataPages,
} from './data-layers.js';
import { TextureChanges } from '../renderable/texture-changes.js';
import { asWebGPUStorageBuffer, isWebGPUStorageBuffer } from './compute.js';
import type { DrawIndirect } from '../compute.js';
import type { WebGPUDrawContext, WebGPULinkable } from './context.js';
import { asWebGPUGlobalBindings, getTextureEntries } from './global-bindings.js';
import type { VertexBufferSlot, WebGPUProgram } from './program.js';
import { getWgslProgram, resolveWgslSource } from './shader/programs.js';
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
  return ctx.programs.get(shaderCode.name, resolveWgslSource(wgsl.source, variantDefines), variantDefines, schema);
}

const getNextProgramItemId = idFactory();

/**
 * The schema the programs of a graphics item are made from: the theme layers (`DataSpec`) read from the data pages
 * (`getDataLayerSchema`); with the render variant define
 */
function getGraphicsSchema(schema: RenderableSchema) {
  const graphicsSchema = getDataLayerSchema(schema);
  for (const s of [schema, graphicsSchema]) {
    if (s.dRenderVariant === undefined) {
      Object.defineProperty(s, 'dRenderVariant', { value: DefineSpec('string') });
    }
  }
  return graphicsSchema;
}

/**
 * Programs of a render item without its data, see `GraphicsContext.createProgramItem`: `resources.preparePipelines`
 * creates them (from the same cache as render items) and their pipelines
 */
export function createWebGPUProgramItem(
  ctx: WebGPUDrawContext,
  drawMode: DrawMode,
  shaderCode: ShaderCode,
  renderableSchema: RenderableSchema,
  defineValues: DefineValues,
  linearTextures: readonly string[],
  transparency: Transparency,
): ProgramItem {
  const id = getNextProgramItemId();
  const schema = getGraphicsSchema(renderableSchema);
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

/**
 * Vertex or index data of a render item: a range of the context's geometry pages (`geometryPages`); of a segmented
 * array, the parts written
 */
type GPUData = { range: BufferRange; usage: number; segmented?: SegmentedArray.Uploaded };

type GPUDataSource = ArrayBufferView | SegmentedArray;

/** Writes `array` into `data`'s range: of a segmented array the parts changed since `data` was written */
function writeGPUData(ctx: WebGPUDrawContext, data: GPUData, array: GPUDataSource) {
  if (SegmentedArray.is(array)) {
    const bpe = array.BYTES_PER_ELEMENT;
    data.segmented = SegmentedArray.writeChanged(array, data.segmented, (part, offset) =>
      ctx.geometryPages.write(data.range, part, offset * bpe),
    );
  } else {
    data.segmented = undefined;
    ctx.geometryPages.write(data.range, array);
  }
}

/** Vertex buffers count as attributes, index buffers as elements (`GraphicsResources.getByteCounts`) */
const byteKind = (usage: number) => (usage & BufferUsage.INDEX ? 'elements' : 'attribute');

function createGPUData(ctx: WebGPUDrawContext, array: GPUDataSource, usage: number): GPUData {
  const range = ctx.geometryPages.allocate(array.byteLength);
  const data: GPUData = { range, usage };
  writeGPUData(ctx, data, array);
  const kind = byteKind(usage);
  ctx.byteCounts[kind] += range.size;
  ctx.stats.resourceCounts[kind] += 1;
  return data;
}

function releaseGPUData(ctx: WebGPUDrawContext, data: GPUData) {
  ctx.geometryPages.free(data.range);
  const kind = byteKind(data.usage);
  ctx.byteCounts[kind] -= data.range.size;
  ctx.stats.resourceCounts[kind] -= 1;
}

/** Rewrites `data` with `array` or moves it to a larger range */
function updateGPUData(ctx: WebGPUDrawContext, data: GPUData, array: GPUDataSource, usage: number): GPUData {
  if (array.byteLength > data.range.size) {
    releaseGPUData(ctx, data);
    return createGPUData(ctx, array, usage);
  }
  writeGPUData(ctx, data, array);
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

/** `shared` receives the textures of specs with `shared`, by name (`.v6/plans/gl/shared-textures.md` §4) */
function createTextures(
  ctx: WebGPUDrawContext,
  schema: RenderableSchema,
  values: TextureValues,
  shared: Map<string, SharedTexture>,
): Textures {
  const textures: Textures = [];
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'texture') continue;
    const value = values[k];
    if (!value) continue;
    if (spec.kind === 'texture') {
      textures.push([k, value.ref.value as Texture]);
    } else if (spec.shared) {
      const held = ctx.sharedTextures.acquire(spec, value.ref.value as SharedTextureImage);
      shared.set(k, held);
      textures.push([k, held.texture]);
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
/** Stands for an absent buffer in the keys of `sharedGroups` */
const NoBuffer = {};
type GroupTrie = { group?: GPUBindGroup; next: WeakMap<object, GroupTrie> };
/**
 * Bind groups binding buffers only, by their layout and buffers (a trie of weak maps: dropped with a destroyed layout
 * or buffer): one group for the items binding the same uniform arenas and storage pages, so that their draws bind one
 * group at their dynamic offsets (`shared-draw-resources.md`, shared object groups)
 */
const sharedGroups = new WeakMap<GPUBindGroupLayout, GroupTrie>();
function getSharedGroup(
  device: GPUDevice,
  layout: GPUBindGroupLayout,
  entries: GPUBindGroupEntry[],
  buffer: GPUBuffer | undefined,
  drawBuffer: GPUBuffer | undefined,
  storageBuffers: readonly GPUBuffer[],
) {
  let node = sharedGroups.get(layout);
  if (!node) sharedGroups.set(layout, (node = { next: new WeakMap() }));
  const step = (key: object) => {
    let child = node!.next.get(key);
    if (!child) node!.next.set(key, (child = { next: new WeakMap() }));
    node = child;
  };
  step(buffer ?? NoBuffer);
  step(drawBuffer ?? NoBuffer);
  for (let i = 0, il = storageBuffers.length; i < il; ++i) step(storageBuffers[i]);
  return (node.group ??= device.createBindGroup({ layout, entries }));
}

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
  /** what `getBindGroup` returns for the group, one object */
  result: { group: GPUBindGroup; offsets: number[] };
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
  renderableSchema: RenderableSchema,
  values: RenderableValues,
  materialId: number,
  transparency: Transparency | undefined,
): RenderItem<T> {
  const id = getNextRenderItemId();
  const { stats } = ctx;

  // compute programs read their textures as they are
  if (renderableSchema.dRenderVariant === undefined) {
    Object.defineProperty(renderableSchema, 'dRenderVariant', { value: DefineSpec('string') });
  }
  const schema = transparency !== undefined ? getGraphicsSchema(renderableSchema) : renderableSchema;
  // the page of each of the item's layers, a define of its programs only where not 0 (all in one): programs of items
  // within the binding limit are those of before, as precompile lists name them (`updateLayers`)
  const layerPages = ValueCell.create(0);
  /** the layout changed since `update` last looked: programs reading the layers from other pages */
  let layerPagesChanged = false;

  const {
    attributeValues,
    defineValues,
    textureValues,
    materialTextureValues,
    dataValues,
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
  /**
   * The theme layers in a range of the data pages (`ctx.dataPages`), one after another (16-byte aligned), their word
   * offsets in the object block; moved to a larger range when one outgrows it. Layers together beyond the binding limit
   * (`ctx.dataBindingSize`) in up to `MaxDataPages` ranges, each within it, bound as pages of their own
   * (`DataPageNames`, the layout in `layerPages`)
   */
  const layers: DataLayer[] =
    schema === renderableSchema ? [] : DataLayerNames.filter((n) => dataValues[n] !== undefined);
  const layerOffsets = layers.map((l) => ValueCell.create(0));
  for (let i = 0; i < layers.length; ++i) objectEntries.push([getDataOffsetName(layers[i]), layerOffsets[i]]);
  const layerVersions = layers.map(() => -1);
  /** of segmented layers (merged objects: the members' layers as parts), the parts written */
  const layerParts: (SegmentedArray.Uploaded | undefined)[] = layers.map(() => undefined);
  /** per page its range, per layer its page and offset in it as written */
  const layerRanges: (BufferRange | undefined)[] = new Array(MaxDataPages).fill(undefined);
  const layerPageOf = layers.map(() => -1);
  const layerPageOffsets = layers.map(() => -1);
  // scratch of `updateLayers`, allocated once: per layer its size, page and offset, per page its bytes, cursor and whether moved
  const nextSizes = layers.map(() => 0);
  const nextPageOf = layers.map(() => 0);
  const nextOffsets = layers.map(() => 0);
  const pageBytes: number[] = new Array(MaxDataPages).fill(0);
  const pageCursors: number[] = new Array(MaxDataPages).fill(0);
  const movedPages: boolean[] = new Array(MaxDataPages).fill(false);
  /** false where a layer alone exceeds the limit, or the layers more than the pages: not drawn (warned) */
  let layersFit = true;
  /**
   * Writes the elements `ranges` (`[start, end)` pairs, all where undefined) of segmented `array` at `offset`: whole
   * words, composed from the parts (a word may hold elements of two)
   */
  function writeSegmentedLayer(
    range: BufferRange,
    array: SegmentedArray,
    ranges: readonly number[] | undefined,
    offset: number,
  ) {
    const bpe = array.BYTES_PER_ELEMENT;
    const write = (start: number, end: number) => {
      const s = Math.floor((start * bpe) / 4) * 4;
      const e = Math.min(array.byteLength, Math.ceil((end * bpe) / 4) * 4);
      if (e > s) ctx.dataPages.write(range, SegmentedArray.read(array, s / bpe, Math.ceil(e / bpe)), offset + s);
    };
    if (ranges) {
      for (let k = 0; k < ranges.length; k += 2) write(ranges[k], ranges[k + 1]);
    } else {
      for (const p of array.parts) if (p.count > 0) write(p.offset, p.offset + p.count);
    }
  }
  const layerImage = (layer: DataLayer) => dataValues[layer].ref.value as DataArray<Uint8Array>;
  /** the bytes of the elements of `layer` (its array may be longer) */
  const layerBytes = (layer: DataLayer) => layerImage(layer).count * DataLayerElementBytes[DataLayers[layer]];
  const alignLayer = (bytes: number) => Math.ceil(bytes / 16) * 16;
  /**
   * The pages of the layers of `sizes` (bytes) into `pageOf` and their bytes into `pageBytes`: all in the first within
   * the limit, else packed largest first; false where they don't fit
   */
  function getLayerPages(sizes: readonly number[], pageOf: number[]) {
    const limit = ctx.dataBindingSize;
    let total = 0;
    for (let i = 0, il = sizes.length; i < il; ++i) total += sizes[i];
    pageBytes.fill(0);
    if (total <= limit) {
      pageOf.fill(0);
      pageBytes[0] = total;
      return true;
    }
    const order = sizes.map((_, i) => i).sort((a, b) => sizes[b] - sizes[a]);
    for (const i of order) {
      const p = pageBytes.findIndex((s) => s + sizes[i] <= limit);
      if (p < 0) return false;
      pageOf[i] = p;
      pageBytes[p] += sizes[i];
    }
    return true;
  }
  function updateLayers() {
    if (layers.length === 0) return;
    // nothing to do where no layer changed (their sizes with them); not fitting is checked again on every update
    if (layersFit) {
      let changed = false;
      for (let i = 0, il = layers.length; i < il; ++i) {
        if (dataValues[layers[i]].ref.version !== layerVersions[i]) {
          changed = true;
          break;
        }
      }
      if (!changed) return;
    }
    const sizes = nextSizes;
    const pageOf = nextPageOf;
    for (let i = 0, il = layers.length; i < il; ++i) sizes[i] = alignLayer(layerBytes(layers[i]));
    if (!getLayerPages(sizes, pageOf)) {
      if (layersFit) {
        ctx.warnUnsupported(
          `the layers of a '${shaderCode.name}' item (${sizes.join(', ')} bytes) exceed ${MaxDataPages} storage ` +
            `bindings of ${ctx.dataBindingSize} bytes`,
        );
      }
      layersFit = false;
      return;
    }
    layersFit = true;
    const pages = getDataLayerPages(layers, pageOf);
    if (pages !== layerPages.ref.value) {
      ValueCell.update(layerPages, pages);
      if (pages === 0) delete defineValues[DataLayerPagesName];
      else defineValues[DataLayerPagesName] = layerPages as unknown as ValueCell<DefineType>; // a number define
      layerPagesChanged = true;
    }
    // per page its layers one after another, a range per page (moved to a larger one when outgrown)
    const pageOffsets = nextOffsets;
    pageCursors.fill(0);
    for (let i = 0, il = layers.length; i < il; ++i) {
      pageOffsets[i] = pageCursors[pageOf[i]];
      pageCursors[pageOf[i]] += sizes[i];
    }
    for (let p = 0; p < MaxDataPages; ++p) {
      const range = layerRanges[p];
      movedPages[p] = false;
      if (range && (pageBytes[p] === 0 || pageBytes[p] > range.size)) {
        ctx.dataPages.free(range);
        ctx.byteCounts.texture -= range.size;
        layerRanges[p] = undefined;
      }
      if (pageBytes[p] > 0 && !layerRanges[p]) {
        const allocated = ctx.dataPages.allocate(pageBytes[p]);
        ctx.byteCounts.texture += allocated.size;
        layerRanges[p] = allocated;
        movedPages[p] = true;
      }
    }
    for (let i = 0; i < layers.length; ++i) {
      const range = layerRanges[pageOf[i]]!;
      const offset = pageOffsets[i];
      // written again where its range or its place in it changed (another layer before it grew)
      const moved = movedPages[pageOf[i]] || layerPageOf[i] !== pageOf[i] || layerPageOffsets[i] !== offset;
      layerPageOf[i] = pageOf[i];
      layerPageOffsets[i] = offset;
      const cell = dataValues[layers[i]];
      const image = layerImage(layers[i]);
      const { array } = image;
      if (SegmentedArray.is(array)) {
        if (moved || cell.ref.version !== layerVersions[i]) {
          const ranges = moved ? undefined : SegmentedArray.changedRanges(array, layerParts[i]);
          if (!ranges) layerParts[i] = SegmentedArray.uploaded(array);
          writeSegmentedLayer(range, array, ranges, offset);
          layerVersions[i] = cell.ref.version;
        }
      } else if (moved || cell.ref.version !== layerVersions[i]) {
        layerParts[i] = undefined;
        // the changed elements only where logged (`TextureChanges`), whole words of the page
        const elementSize = DataLayerElementBytes[DataLayers[layers[i]]];
        // whole words of the elements (the array may be longer, reused; its length is whole words)
        const byteCount = Math.min(array.byteLength, Math.ceil(layerBytes(layers[i]) / 4) * 4);
        const ranges = moved ? undefined : TextureChanges.get(image, layerVersions[i], cell.ref.version);
        if (ranges) {
          for (let k = 0; k < ranges.length; k += 2) {
            const start = Math.floor((ranges[k] * elementSize) / 4) * 4;
            const end = Math.min(byteCount, Math.ceil((ranges[k + 1] * elementSize) / 4) * 4);
            const bytes = new Uint8Array(array.buffer, array.byteOffset + start, end - start);
            ctx.dataPages.write(range, bytes, offset + start);
          }
        } else {
          ctx.dataPages.write(range, new Uint8Array(array.buffer, array.byteOffset, byteCount), offset);
        }
        layerVersions[i] = cell.ref.version;
      }
      ValueCell.updateIfChanged(layerOffsets[i], (range.offset + offset) / 4);
    }
  }
  updateLayers();
  // the programs are created after, with the layout
  layerPagesChanged = false;

  const objectLayout = getUniformBlockLayout('object', schema);
  const materialLayout = getUniformBlockLayout('material', schema);
  /** the object block and the draw blocks, in one buffer: one bind group (2) with dynamic offsets */
  const arena = createUniformArena(ctx.ring, ctx.release, ctx.uniformPages);
  const objectBlock = objectLayout
    ? ownBlock(arena.createSnapshot(objectLayout.size), objectLayout, objectEntries)
    : undefined;
  const material = materialLayout ? acquireMaterial(ctx, materialId, materialLayout) : undefined;
  const materialBlock = material ? ownBlock(material.snapshot, materialLayout!, materialEntries) : undefined;
  objectBlock?.write();
  materialBlock?.write();

  // 'draw' uniforms of the render item and of the global bindings, per global bindings instance; multi-draw
  // levels (`MultiDrawBaseData`) have their own block with their uniforms (e.g. `uLod`) replacing the item's; the
  // draws of solid interior caps a block per cap (`GlobalBindings.setSolidInteriorCap`), its values constant, so that
  // the item's other draws keep their slot and a recorded pass replays
  type CapBlocks = Map<number, Block>;
  type DrawBlocks = {
    base?: Block;
    levels: WeakMap<MultiDrawBaseData, Block>;
    caps: CapBlocks;
    levelCaps: WeakMap<MultiDrawBaseData, CapBlocks>;
  };
  type DrawGlobals = { drawUniforms: UniformsList; solidInteriorCap?: { pass: number; clip: number } };
  const drawBlocks = new WeakMap<DrawGlobals, DrawBlocks>();
  const noGlobals: DrawGlobals = { drawUniforms: [] };
  function getDrawBlock(program: WebGPUProgram, globals: DrawGlobals, level?: MultiDrawBaseData) {
    if (!program.drawLayout) return;
    let blocks = drawBlocks.get(globals);
    if (!blocks) {
      blocks = { levels: new WeakMap(), caps: new Map(), levelCaps: new WeakMap() };
      drawBlocks.set(globals, blocks);
    }
    const { members, size } = program.drawLayout;
    // the layout and entries for a block created: not per draw
    const getLayout = (): UniformBlockLayout => ({ block: 'object', members, size, hash: '' });
    /** the item's or the level's entries */
    const getEntries = (): UniformsList => {
      if (!level) return [...drawEntries, ...globals.drawUniforms];
      const names = new Set(level.uniforms.map(([k]) => k));
      return [...drawEntries.filter(([k]) => !names.has(k)), ...level.uniforms, ...globals.drawUniforms];
    };
    const cap = globals.solidInteriorCap;
    if (cap && cap.pass !== 0) {
      let caps = level ? blocks.levelCaps.get(level) : blocks.caps;
      if (!caps) blocks.levelCaps.set(level!, (caps = new Map()));
      // the clip index is -1 (near plane) or a clip object's
      const key = cap.pass * 256 + cap.clip + 1;
      let block = caps.get(key);
      if (!block) {
        const entries: UniformsList = getEntries().map(([k, v]) =>
          k === 'uSolidInteriorPass'
            ? [k, ValueCell.create(cap.pass)]
            : k === 'uSolidInteriorClip'
              ? [k, ValueCell.create(cap.clip)]
              : [k, v],
        );
        block = ownBlock(arena.createSnapshot(size), getLayout(), entries);
        caps.set(key, block);
      }
      return block;
    }
    if (!level) return (blocks.base ??= ownBlock(arena.createSnapshot(size), getLayout(), getEntries()));
    let block = blocks.levels.get(level);
    if (!block) {
      // levels are kept by the renderable, replaced when their capacity grows (the slot of a replaced one stays unused)
      block = ownBlock(arena.createSnapshot(size), getLayout(), getEntries());
      blocks.levels.set(level, block);
    }
    return block;
  }

  /** per vertex buffer slot of the programs, created and updated on use; `data` is unset until uploaded */
  const vertexBuffers = new Map<string, { cells: ValueCell<any>[]; versions: number[]; data?: GPUData }>();
  const slotEntries = new Map<VertexBufferSlot, NonNullable<ReturnType<typeof vertexBuffers.get>>>();
  function getVertexBuffer(slot: VertexBufferSlot) {
    // by slot object first (stable per program), by key across programs
    let entry = slotEntries.get(slot);
    if (!entry) {
      entry = vertexBuffers.get(slot.key);
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
    }
    const { cells, versions } = entry;
    // written on the GPU (instanced attributes of GPU-culled merged objects): bound as they are
    const value = cells.length === 1 ? cells[0].ref.value : undefined;
    if (isWebGPUStorageBuffer(value)) {
      const buffer = value.getGPUBuffer();
      return { buffer, offset: 0, size: buffer.size };
    }
    let changed = !entry.data;
    for (let i = 0, il = cells.length; i < il; ++i) {
      if (cells[i].ref.version !== versions[i]) changed = true;
    }
    if (!changed) return entry.data!.range;
    const arrays = cells.map((c) => c.ref.value as Float32Array | SegmentedArray<Float32Array>);
    const array =
      arrays.length === 1
        ? arrays[0]
        : interleave(
            arrays.map((a) => (SegmentedArray.is(a) ? SegmentedArray.flatten(a) : a)),
            slot.itemSizes,
          );
    entry.data = entry.data
      ? updateGPUData(ctx, entry.data, array, BufferUsage.VERTEX)
      : createGPUData(ctx, array, BufferUsage.VERTEX);
    for (let i = 0, il = cells.length; i < il; ++i) versions[i] = cells[i].ref.version;
    return entry.data.range;
  }
  let indexBuffer: GPUData | undefined;
  /** quads without elements of their own draw indexed by the shared quad indices (`quad-indices.ts`) */
  const quads = !!shaderCode.quads && !values.elements?.ref.value && isQuadIndexing();
  function createIndexBuffer() {
    if (values.elements?.ref.value) {
      indexBuffer = createGPUData(ctx, values.elements.ref.value, BufferUsage.INDEX);
    }
  }
  createIndexBuffer();

  // the textures of specs with `shared`, by name, released instead of destroyed (`.v6/plans/gl/shared-textures.md` §4)
  const sharedTextures = new Map<string, SharedTexture>();
  const textures = createTextures(ctx, schema, textureValues, sharedTextures);
  const materialTextures = createTextures(ctx, schema, materialTextureValues, sharedTextures);
  // the versions of the values as last seen by `update`, per list in its order: compared for every item on every update,
  // not looked up by name in an object of all the values (a dictionary of a hundred keys); the layers' define is not
  // among the defines (`layerPagesChanged`)
  const defineEntries = Object.entries(defineValues);
  const defineVersions = defineEntries.map(([, v]) => v.ref.version);
  const textureVersions = textures.map(([k]) => textureValues[k].ref.version);
  const materialTextureVersions = materialTextures.map(([k]) => materialTextureValues[k].ref.version);
  let drawCountVersion = values.drawCount.ref.version;
  let instanceCountVersion = values.instanceCount.ref.version;
  let elementsVersion = values.elements?.ref.version ?? -1;
  /** loads the changed textures of `list`; whether one was replaced by another */
  function updateTextures(list: Textures, textureValues: TextureValues, versions: number[]) {
    let replaced = false;
    for (let i = 0, il = list.length; i < il; ++i) {
      const [k, texture] = list[i];
      const value = textureValues[k];
      if (value.ref.version === versions[i]) continue;
      const held = sharedTextures.get(k);
      if (held) {
        // the same image uploads once for all its items, another image is another texture (§4)
        const next = updateSharedTexture(
          ctx.sharedTextures,
          held,
          schema[k] as TextureSpec<TextureKind>,
          value.ref.value as SharedTextureImage,
        );
        if (next !== held) {
          sharedTextures.set(k, next);
          list[i][1] = next.texture;
          replaced = true;
        }
      } else if (schema[k].kind !== 'texture') {
        // textures of kind 'texture' are managed externally
        texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
      } else if (list[i][1] !== value.ref.value) {
        list[i][1] = value.ref.value as Texture;
        replaced = true;
      }
      versions[i] = value.ref.version;
    }
    return replaced;
  }

  let drawCount: number = values.drawCount.ref.value;
  let instanceCount: number = values.instanceCount.ref.value;
  /** the draw arguments and instance count last written into them */
  let indirectWritten: { args: DrawIndirect; instanceCount: number } | undefined;
  function getDrawIndirect() {
    const args = (values.drawIndirect as ValueCell<DrawIndirect | undefined> | undefined)?.ref.value;
    if (!args) return undefined;
    if (indexBuffer) notSupported('indirect draws of indexed geometry');
    if (args.count !== undefined) return args;
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
  /**
   * per program the cells of its storage bindings (a number: the data page of that index, `DataPageNames`), resolved
   * once, not by name per draw
   */
  const storageCells = new WeakMap<WebGPUProgram, (ValueCell<any> | number)[]>();
  function getStorageBuffers(program: WebGPUProgram) {
    let cells = storageCells.get(program);
    if (!cells) {
      cells = program.storageBindings.map(({ name }) => {
        const page = (DataPageNames as readonly string[]).indexOf(name);
        if (page >= 0) return page;
        const cell = storageBufferValues[name];
        if (!cell) throw new Error(`missing storage buffer '${name}' of '${shaderCode.name}'`);
        return cell;
      });
      storageCells.set(program, cells);
    }
    currentStorageBuffers.length = cells.length;
    for (let i = 0, il = cells.length; i < il; ++i) {
      const cell = cells[i];
      if (typeof cell === 'number') {
        currentStorageBuffers[i] = layerRanges[cell]!.buffer;
        continue;
      }
      const buffer = cell.ref.value;
      if (!buffer)
        throw new Error(`missing storage buffer '${program.storageBindings[i].name}' of '${shaderCode.name}'`);
      currentStorageBuffers[i] = asWebGPUStorageBuffer(buffer).getGPUBuffer();
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
      // an object group binding buffers only is shared by the items binding the same (their slots dynamic offsets)
      const shared = index === 2 && textures.length === 0;
      cached = {
        group: shared
          ? getSharedGroup(ctx.device, layout, entries, slot?.buffer, drawSlot?.buffer, storageBuffers)
          : ctx.device.createBindGroup({ layout, entries }),
        buffer: slot?.buffer,
        drawBuffer: drawSlot?.buffer,
        bindingsVersion,
        textures,
        views: textures.map((t) => t?.getView()),
        storageBuffers: storageBuffers.slice(),
        offsets: [],
        result: { group: undefined as unknown as GPUBindGroup, offsets: noOffsets },
      };
      cached.result.group = cached.group;
      groups.set(layout, cached);
    }
    // dynamic offsets in binding order: the uniform block (0), then the draw uniforms
    const { offsets, result } = cached;
    offsets.length = 0;
    if (slot) offsets.push(slot.offset);
    if (drawSlot) offsets.push(drawSlot.offset);
    // one result per cached group (read before the next call)
    result.offsets = offsets.length ? offsets : noOffsets;
    return result;
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
      for (const b of vertexBuffers.values()) count += b.data?.range.size ?? 0;
      if (indexBuffer) count += indexBuffer.range.size;
      // a shared texture by the item's share (§5)
      for (const [k, t] of textures) count += sharedTextures.get(k)?.getByteCount() ?? t.getByteCount();
      for (const [k, t] of materialTextures) count += sharedTextures.get(k)?.getByteCount() ?? t.getByteCount();
      for (const range of layerRanges) if (range) count += range.size;
      return count;
    },
    getProgram: getGraphicsProgram,
    setTransparency: (value: Transparency) => {
      if (value === transparency) return;
      transparency = value;
      invalidatePrograms(true);
    },

    render: (variant: T, globals?: GlobalBindings, mdbDataList?: MultiDrawBaseData[]) => {
      if (drawCount === 0 || instanceCount === 0 || !layersFit) return;
      const gpuGlobals = globals ? asWebGPUGlobalBindings(globals) : undefined;
      const program = getWebGPUProgram(variant);
      if (!program) return;

      const pipeline = ctx.getPipeline(program, drawMode);
      // compiling in the background, drawn once ready
      if (!pipeline) {
        stats.skippedDraws += 1;
        return;
      }
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
        const range = getVertexBuffer(slots[i]);
        state.setVertexBuffer(i, range.buffer, range.offset, range.size);
      }
      if (indexBuffer) {
        const { range } = indexBuffer;
        state.setIndexBuffer(range.buffer, 'uint32', range.offset, range.size);
      }
      // quads without elements of their own draw indexed by the shared quad indices (`quad-indices.ts`), unindexed
      // past their limit
      const quad = quads ? ctx.quadIndices.get(drawCount / QuadIndicesPerItem) : undefined;
      if (quad) state.setIndexBuffer(quad.buffer, 'uint32', 0, quad.buffer.size);
      const indexed = !!indexBuffer || !!quad;

      if (mdbDataList) {
        // culled and level-of-detail ranges, `firstInstance` offsets the instanced attributes
        for (const level of mdbDataList) {
          if (level.count === 0) continue;
          const drawBlock = getDrawBlock(program, gpuGlobals ?? noGlobals, level);
          const object = getBindGroup(program, 2, objectBlock, textures, drawBlock);
          state.setBindGroup(2, object.group, object.offsets);
          const { counts, instanceCounts, firsts, offsets, baseVertices, baseInstances, indirect } = level;
          if (indirect) {
            // arguments written on the GPU (GPU culling): a draw per record, drawn or not
            const buffer = asWebGPUStorageBuffer(indirect.buffer).getGPUBuffer();
            // indexed records of quads need the shared indices bound (within their limit, see the culler)
            const indexedRecords = indirect.indexed ?? !!indexBuffer;
            if (indexedRecords && !indexed) continue;
            const stride = indexedRecords ? IndirectRecordBytes.indexed : IndirectRecordBytes.vertex;
            for (let i = 0; i < level.count; ++i) {
              if (indexedRecords) pass.drawIndexedIndirect(buffer, indirect.byteOffset + i * stride);
              else pass.drawIndirect(buffer, indirect.byteOffset + i * stride);
            }
            stats.calls.drawInstancedBase += level.count;
            continue;
          }
          for (let i = 0; i < level.count; ++i) {
            if (counts[i] === 0 || instanceCounts[i] === 0) continue;
            if (indexed) {
              // element offsets are in bytes (those of quads are their first vertex times 4)
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
        if (indirect?.count !== undefined) pass.draw(indirect.count, instanceCount);
        else if (indirect)
          pass.drawIndirect(asWebGPUStorageBuffer(indirect.buffer).getGPUBuffer(), indirect.byteOffset);
        else if (indexed) pass.drawIndexed(drawCount, instanceCount);
        else pass.draw(drawCount, instanceCount);
        stats.calls.drawInstanced += 1;
        stats.calls.counts += instanceCount;
      }
    },

    update: () => {
      let definesChanged = false;
      for (let i = 0, il = defineEntries.length; i < il; ++i) {
        const v = defineEntries[i][1];
        if (v.ref.version !== defineVersions[i]) {
          definesChanged = true;
          defineVersions[i] = v.ref.version;
        }
      }
      if (definesChanged) invalidatePrograms(true);

      if (values.drawCount.ref.version !== drawCountVersion) {
        stats.drawCount += values.drawCount.ref.value - drawCount;
        stats.instancedDrawCount += instanceCount * (values.drawCount.ref.value - drawCount);
        drawCount = values.drawCount.ref.value;
        drawCountVersion = values.drawCount.ref.version;
      }
      if (values.instanceCount.ref.version !== instanceCountVersion) {
        stats.instanceCount += values.instanceCount.ref.value - instanceCount;
        stats.instancedDrawCount += (values.instanceCount.ref.value - instanceCount) * drawCount;
        instanceCount = values.instanceCount.ref.value;
        instanceCountVersion = values.instanceCount.ref.version;
      }

      if (indexBuffer && values.elements.ref.version !== elementsVersion) {
        indexBuffer = updateGPUData(ctx, indexBuffer, values.elements.ref.value, BufferUsage.INDEX);
        elementsVersion = values.elements.ref.version;
      }

      let texturesReplaced = updateTextures(textures, textureValues, textureVersions);
      if (updateTextures(materialTextures, materialTextureValues, materialTextureVersions)) texturesReplaced = true;
      updateLayers();
      // the layers' page layout changed (`layerPages`): programs reading them from other pages
      if (layerPagesChanged) {
        layerPagesChanged = false;
        if (!definesChanged) invalidatePrograms(true);
      }
      // another texture may be filtered differently (a program of other linear reads)
      if (texturesReplaced && !definesChanged) invalidatePrograms(false);
    },

    destroy: () => {
      if (destroyed) return;
      for (const b of vertexBuffers.values()) if (b.data) releaseGPUData(ctx, b.data);
      if (indexBuffer) releaseGPUData(ctx, indexBuffer);
      for (const [k, t] of [...textures, ...materialTextures]) {
        const held = sharedTextures.get(k);
        if (held) held.release();
        else if (schema[k].kind !== 'texture') t.destroy();
      }
      sharedTextures.clear();
      for (let p = 0; p < layerRanges.length; ++p) {
        const range = layerRanges[p];
        if (!range) continue;
        ctx.dataPages.free(range);
        ctx.byteCounts.texture -= range.size;
        layerRanges[p] = undefined;
      }
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
