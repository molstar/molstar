/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isQuadIndexing, QuadIndicesPerItem } from '../quad-indices.js';
import type { UniformBlockWriter } from '../uniform-layout.js';
import { createAttributeBuffers, type ElementsBuffer, type AttributeBuffers } from './buffer.js';
import type { AttributeKind } from '../buffer.js';
import { createTextures, type GLTexture } from './texture.js';
import type { Texture, TextureKind } from '../texture.js';
import { type SharedTexture, type SharedTextureImage, updateSharedTexture } from '../shared-texture.js';
import {
  ComputeRenderVariants,
  type ComputeRenderVariant,
  type DrawMode,
  GraphicsRenderVariants,
  type GraphicsRenderVariant,
  type GlobalBindings,
  type MultiDrawBaseData,
  type ProgramItem,
  type RenderItem,
  type Transparency,
} from '../render-item.js';
import type { WebGLContext } from './context.js';
import type { ShaderCode, DefineValues } from '../shader-code.js';
import { getContextDefineLookup } from './shader-code.js';
import { getRenderPassTargetHeight } from './render-pass-target.js';
import type { Program, Programs } from './program.js';
import {
  type RenderableSchema,
  type RenderableValues,
  type AttributeSpec,
  isUsed,
  splitValues,
  DefineSpec,
  TextureSpec,
  UniformSpec,
} from '../renderable/schema.js';
import { idFactory } from '@molstar/core/util/id-factory';
import { ValueCell } from '@molstar/core/util';
import { createTextureImage, type DataArray, type TextureImage, type TextureVolume } from '../renderable/util.js';
import {
  getDataTextureDimName,
  getDataTextureItemSize,
  getDataTextureNames,
  getDataTextureOffsetName,
  getDataTextureSchema,
} from './data-textures.js';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { SegmentedArray } from '../segmented-array.js';
import type { DataTexturePages, DataTextureRange } from './data-pages.js';
import { InstanceDataVec4s, packInstanceData } from '../renderable/merged.js';
import { TextureChanges } from '../renderable/texture-changes.js';
import { isDebugMode, isTimingMode } from '@molstar/core/util/debug';
import type { VertexArray } from './vertex-array.js';
import { fillSerial } from '@molstar/core/util/array';
import { deepClone } from '@molstar/core/util/object';
import { cloneUniformValues, type UniformBuffer, type UniformKind, type UniformsList } from '../uniform.js';
import { getUniformBlockLayout } from '../uniform-layout.js';
import { getAttributeLayout } from '../attribute-layout.js';
import { checkRenderPass } from './render-pass.js';
import { asWebGLGlobalBindings } from './global-bindings.js';

// Handle Firefox's preference [webgl.max-vert-ids-per-draw] which defaults to 30_000_000
// since FF119, see https://bugzilla.mozilla.org/show_bug.cgi?id=1849433
const MaxDrawCount = 30_000_000;

const getNextRenderItemId = idFactory();

export function getDrawMode(ctx: WebGLContext, drawMode: DrawMode) {
  const { gl } = ctx;
  switch (drawMode) {
    case 'points':
      return gl.POINTS;
    case 'lines':
      return gl.LINES;
    case 'line-strip':
      return gl.LINE_STRIP;
    case 'line-loop':
      return gl.LINE_LOOP;
    case 'triangles':
      return gl.TRIANGLES;
    case 'triangle-strip':
      return gl.TRIANGLE_STRIP;
    case 'triangle-fan':
      return gl.TRIANGLE_FAN;
  }
}

export interface WebGLRenderItem<T extends string> extends RenderItem<T> {
  getProgram: (variant: T) => Program;
}

//

function createProgramVariant(
  ctx: WebGLContext,
  variant: string,
  defineValues: DefineValues,
  shaderCode: ShaderCode,
  schema: RenderableSchema,
) {
  defineValues = { ...defineValues, dRenderVariant: ValueCell.create(variant) };
  if (schema.dRenderVariant === undefined) {
    Object.defineProperty(schema, 'dRenderVariant', { value: DefineSpec('string') });
  }
  return ctx.resources.program(defineValues, shaderCode, schema);
}

//

function createValueChanges() {
  return {
    attributes: false,
    defines: false,
    elements: false,
    textures: false,
  };
}
type ValueChanges = ReturnType<typeof createValueChanges>;

function resetValueChanges(valueChanges: ValueChanges) {
  valueChanges.attributes = false;
  valueChanges.defines = false;
  valueChanges.elements = false;
  valueChanges.textures = false;
}

//

function getRenderVariant(variant: string, transparency: Transparency | undefined): string {
  if (variant === 'color') {
    switch (transparency) {
      case 'blended':
        return 'colorBlended';
      case 'wboit':
        return 'colorWboit';
      case 'dpoit':
        return 'colorDpoit';
    }
  }
  return variant;
}

const getNextProgramItemId = idFactory();

/**
 * Programs of a render item without its data, see `GraphicsContext.createProgramItem`: created as render items create
 * theirs (the same cache), linked by `resources.preparePipelines`; WebGL programs don't read linear filtering
 */
export function createWebGLProgramItem(
  ctx: WebGLContext,
  shaderCode: ShaderCode,
  schema: RenderableSchema,
  defineValues: DefineValues,
  transparency: Transparency,
): ProgramItem {
  const programs = new Map<GraphicsRenderVariant, Program>();
  // the schema and defines render items make their programs with: the same programs
  const adapted = getProgramItemSchema(ctx, schema, defineValues);
  const getProgram = (variant: GraphicsRenderVariant) => {
    let program = programs.get(variant);
    if (!program) {
      const renderVariant = getRenderVariant(variant, transparency);
      program = createProgramVariant(ctx, renderVariant, adapted.defineValues, shaderCode, adapted.schema);
      programs.set(variant, program);
    }
    return program;
  };
  return {
    id: getNextProgramItemId(),
    getProgram,
    destroy: () => {
      for (const program of programs.values()) program.destroy();
      programs.clear();
    },
  };
}

export function createGraphicsRenderItem(
  ctx: WebGLContext,
  drawMode: DrawMode,
  shaderCode: ShaderCode,
  schema: RenderableSchema,
  values: RenderableValues,
  materialId: number,
  transparency: Transparency,
): WebGLRenderItem<GraphicsRenderVariant> {
  return createRenderItem(ctx, drawMode, shaderCode, schema, values, materialId, GraphicsRenderVariants, transparency);
}

export function createComputeRenderItem(
  ctx: WebGLContext,
  drawMode: DrawMode,
  shaderCode: ShaderCode,
  schema: RenderableSchema,
  values: RenderableValues,
  materialId = -1,
): WebGLRenderItem<ComputeRenderVariant> {
  return createRenderItem(ctx, drawMode, shaderCode, schema, values, materialId, ComputeRenderVariants, undefined);
}

/**
 * Vertices whose emulated vertex id (`aVertex`, a float attribute where `gl_VertexID` is missing: WebGL1) is exact:
 * 2^24, e.g. about 2.8M points (6 vertices each); a known limit of WebGL1
 */
const MaxEmulatedVertexCount = 2 ** 24;

/** A data array of a render item in its range of the data pages, with its uniforms */
type DataTextureState = {
  readonly pages: DataTexturePages;
  /** components per element */
  readonly itemSize: number;
  range: DataTextureRange;
  readonly dim: ValueCell<Vec2>;
  readonly offset: ValueCell<number>;
  /** of a segmented array, its parts written */
  parts: SegmentedArray.Uploaded | undefined;
};

/**
 * The schema and define values the programs of a render item are made with (`createRenderItem`), from those a program
 * item is given (the renderable's): the loop bounds as number defines on WebGL1, their values from the defines that
 * mirror them (`dClipObjectCount` for `uClipObjectCount`; throws for one without, as the direct volume's step count,
 * which only the data gives); the instance data of merged objects on WebGL2 without base instances; the data arrays as
 * data textures. Define values in the adapted schema's order, as render items have them (`splitValues`)
 */
function getProgramItemSchema(ctx: WebGLContext, schema: RenderableSchema, defineValues: DefineValues) {
  let adapted = schema;
  let values = defineValues;
  if (!ctx.isWebGL2) {
    adapted = getWebGL1Schema(schema);
    if (adapted !== schema) {
      values = {};
      for (const k of Object.keys(adapted)) {
        if (adapted[k].type !== 'define' || !(k in defineValues || schema[k].type === 'uniform')) continue;
        if (schema[k].type === 'uniform') {
          const mirror = defineValues[`d${k.slice(1)}`];
          if (!mirror) throw new Error(`the value of the loop bound '${k}' is not among the defines (WebGL1)`);
          values[k] = ValueCell.create(mirror.ref.value);
        } else {
          values[k] = defineValues[k];
        }
      }
    }
  }
  const { multiDrawInstancedBaseVertexBaseInstance, drawInstancedBaseVertexBaseInstance } = ctx.extensions;
  const instanceData =
    ctx.isWebGL2 &&
    !multiDrawInstancedBaseVertexBaseInstance &&
    !drawInstancedBaseVertexBaseInstance &&
    !!values.dSegmented?.ref.value &&
    adapted.aSegment !== undefined;
  if (instanceData) {
    adapted = getInstanceDataSchema(adapted);
    values = { ...values, dInstanceData: ValueCell.create(true) };
  }
  return { schema: getDataTextureSchema(adapted), defineValues: values };
}

/** The WebGL1 schemas of schemas, one per schema: the program cache keys declarations by the schema object */
const webGL1Schemas = new WeakMap<RenderableSchema, RenderableSchema>();

/** Loop bounds (`UniformSpec.loopBound`) as number defines: WebGL1 loops need constant bounds */
function getWebGL1Schema(schema: RenderableSchema): RenderableSchema {
  const cached = webGL1Schemas.get(schema);
  if (cached) return cached;
  let adapted: { [k: string]: RenderableSchema[string] } | undefined;
  for (const k of Object.keys(schema)) {
    const spec = schema[k];
    if (spec.type !== 'uniform' || !(spec as UniformSpec<UniformKind>).loopBound) continue;
    adapted ??= { ...schema };
    adapted[k] = DefineSpec('number');
  }
  const result = adapted ?? schema;
  webGL1Schemas.set(schema, result);
  return result;
}

/** The instance data schemas of schemas, one per schema: the program cache keys declarations by the schema object */
const instanceDataSchemas = new WeakMap<RenderableSchema, RenderableSchema>();

/**
 * Merged render items (`dSegmented`) on WebGL2 without base instances read their instance data (transforms, instance
 * ids, segment data) from a texture at the draw's first instance (`uInstanceBase`) plus `gl_InstanceID`
 * (`dInstanceData`, see `common-vert-params`): re-pointing the instanced attributes per draw of a segment cost about
 * 40% of a frame of 1000 merged units
 */
export function getInstanceDataSchema(schema: RenderableSchema): RenderableSchema {
  let result = instanceDataSchemas.get(schema);
  if (!result) {
    result = {
      ...schema,
      tInstanceData: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
      uInstanceBase: UniformSpec('i', undefined, 'draw'),
      dInstanceData: DefineSpec('boolean'),
    };
    instanceDataSchemas.set(schema, result);
  }
  return result;
}

/** Writes the instance data of merged `values` into the image of `tInstanceData`, grown where too small */
function updateInstanceData(values: RenderableValues) {
  const count: number = values.instanceCount.ref.value;
  const image = createTextureImage(
    Math.max(1, count * InstanceDataVec4s),
    4,
    Float32Array,
    (values.tInstanceData.ref.value as TextureImage<Float32Array>).array,
  );
  packInstanceData(values, image.array);
  ValueCell.update(values.tInstanceData, image);
}

/**
 * Creates a render item
 *
 * - assumes that `values.drawCount` and `values.instanceCount` exist
 */
export function createRenderItem<T extends string>(
  ctx: WebGLContext,
  drawMode: DrawMode,
  shaderCode: ShaderCode,
  schema: RenderableSchema,
  values: RenderableValues,
  materialId: number,
  renderVariants: T[],
  transparency: Transparency | undefined,
): WebGLRenderItem<T> {
  const id = getNextRenderItemId();
  const { stats, state, resources } = ctx;
  if (!ctx.isWebGL2) schema = getWebGL1Schema(schema);
  const {
    instancedArrays,
    vertexArrayObject,
    multiDrawInstancedBaseVertexBaseInstance,
    drawInstancedBaseVertexBaseInstance,
  } = ctx.extensions;
  const instanceData =
    ctx.isWebGL2 &&
    !multiDrawInstancedBaseVertexBaseInstance &&
    !drawInstancedBaseVertexBaseInstance &&
    !!values.dSegmented?.ref.value &&
    schema.aSegment !== undefined;
  /** the versions of the instance values `tInstanceData` was written from */
  const instanceDataVersions: number[] = [];
  const instanceDataKeys = ['instanceCount', 'aTransform', 'aInstance', 'aSegment', 'aSegmentSphere', 'aSegmentLod'];
  const syncInstanceData = () => {
    let changed = false;
    for (let i = 0; i < instanceDataKeys.length; ++i) {
      const v = values[instanceDataKeys[i]].ref.version;
      if (v !== instanceDataVersions[i]) {
        instanceDataVersions[i] = v;
        changed = true;
      }
    }
    if (changed) updateInstanceData(values);
  };
  if (instanceData) {
    schema = getInstanceDataSchema(schema);
    values = {
      ...values,
      tInstanceData: ValueCell.create(createTextureImage(1, 4, Float32Array)),
      uInstanceBase: ValueCell.create(0),
      dInstanceData: ValueCell.create(true),
    };
    syncInstanceData();
  }

  // the data arrays (`DataSpec`) in ranges of the context's data pages, the texture's size and the offset object uniforms
  // of this backend (`data-pages.ts`, `data-textures.ts`)
  const dataTextures = new Map<string, DataTextureState>();
  const dataSchema = getDataTextureSchema(schema);
  if (dataSchema !== schema) {
    schema = dataSchema;
    const uniforms: { [k: string]: ValueCell<any> } = {};
    for (const name of getDataTextureNames(schema)) {
      if (!values[name]) continue;
      const format = (schema[name] as TextureSpec<TextureKind>).data!;
      const pages = resources.dataPages(format);
      const range = pages.allocate((values[name].ref.value as DataArray<any>).count);
      const dim = ValueCell.create(Vec2.create(range.width, range.height));
      const offset = ValueCell.create(range.offset);
      uniforms[getDataTextureDimName(name)] = dim;
      uniforms[getDataTextureOffsetName(name)] = offset;
      dataTextures.set(name, { pages, itemSize: getDataTextureItemSize(format), range, dim, offset, parts: undefined });
    }
    values = { ...values, ...uniforms };
  }
  /**
   * Writes the element `ranges` of `data` (all where undefined) into the range of the data array `name`, moved to a
   * larger one when outgrown; of a segmented array (merged objects), the parts changed since its last write
   */
  const writeData = (name: string, data: DataArray<any>, ranges?: readonly number[]) => {
    const d = dataTextures.get(name)!;
    let all = ranges === undefined;
    if (data.count > d.range.size) {
      d.pages.free(d.range);
      d.range = d.pages.allocate(data.count);
      ValueCell.updateIfChanged(d.offset, d.range.offset);
      if (d.dim.ref.value[0] !== d.range.width || d.dim.ref.value[1] !== d.range.height) {
        ValueCell.update(d.dim, Vec2.create(d.range.width, d.range.height));
      }
      all = true;
    }
    const array = data.array as unknown;
    if (SegmentedArray.is(array)) {
      // the parts' ranges are by component, the texels by element
      const parts = all ? undefined : SegmentedArray.changedRanges(array, d.parts);
      if (parts) {
        const elements: number[] = [];
        for (let k = 0; k < parts.length; k += 2) {
          elements.push(Math.floor(parts[k] / d.itemSize), Math.ceil(parts[k + 1] / d.itemSize));
        }
        d.pages.write(d.range, data, elements);
      } else {
        d.parts = SegmentedArray.uploaded(array);
        d.pages.write(d.range, data);
      }
    } else {
      d.parts = undefined;
      d.pages.write(d.range, data, all ? undefined : ranges);
    }
  };

  // filter out unsupported variants
  renderVariants = renderVariants.filter((v) => {
    if (v === 'tracing') return !!ctx.extensions.drawBuffers;
    return true;
  });

  // emulate gl_VertexID when needed
  const contextDefines = getContextDefineLookup(ctx.extensions);
  const emulateVertexId = schema.aVertex?.type === 'attribute' && isUsed(schema.aVertex, contextDefines);
  // without elements, every drawn vertex needs an id (e.g. 6 per point quad with one vertex of data per point)
  const getEmulatedVertexCount = () =>
    schema.elements
      ? values.uVertexCount.ref.value
      : Math.max(values.uVertexCount.ref.value, values.drawCount.ref.value);
  /** The serial vertex ids of `aVertex` for at least `vertexCount` vertices, see `MaxEmulatedVertexCount` */
  const ensureVertexIds = (vertexCount: number) => {
    if (values.aVertex.ref.value.length >= vertexCount) return;
    if (vertexCount > MaxEmulatedVertexCount) {
      console.warn(
        `${vertexCount} vertices with an emulated vertex id (WebGL1): ids above ${MaxEmulatedVertexCount} are not exact`,
      );
    }
    ValueCell.update(values.aVertex, fillSerial(new Float32Array(vertexCount)));
  };
  if (emulateVertexId) ensureVertexIds(getEmulatedVertexCount());

  const {
    attributeValues,
    defineValues,
    textureValues,
    materialTextureValues,
    uniformValues,
    materialUniformValues,
    bufferedUniformValues,
    // the projections mapped to the context's clip depth, window rectangles flipped to GL's bottom-left origin
  } = splitValues(schema, values, {
    depthRange: ctx.depthRangeConversion,
    targetHeight: () => getRenderPassTargetHeight(ctx),
  });

  // on WebGL2, uniforms assigned to a block are backed by uniform buffers
  const useUniformBuffers = ctx.uniformBlocks;
  const getBlock = (k: string) => {
    const { block } = schema[k] as UniformSpec<UniformKind>;
    return !useUniformBuffers || block === 'draw' ? undefined : block;
  };

  const uniformValueEntries: UniformsList = [];
  const objectBlockValueEntries: UniformsList = [];
  for (const e of Object.entries(uniformValues)) {
    const block = getBlock(e[0]);
    if (block === undefined) {
      uniformValueEntries.push(e);
    } else if (block === 'object') {
      objectBlockValueEntries.push(e);
    } else {
      throw new Error(`uniform '${e[0]}' of block '${block}' cannot be set by a render item`);
    }
  }
  // 'material' variant uniforms are either classic or in the 'material' block, see `UniformSpec`
  const materialUniformValueEntries: UniformsList = [];
  const materialBlockValueEntries: UniformsList = [];
  for (const e of Object.entries(materialUniformValues)) {
    if (getBlock(e[0]) === 'material') {
      materialBlockValueEntries.push(e);
    } else {
      materialUniformValueEntries.push(e);
    }
  }
  // 'buffered' variant uniforms are either classic or in the 'object' block, see `UniformSpec`
  const backBufferUniformValueEntries = Object.entries(bufferedUniformValues);
  const frontBufferUniformValueEntries = Object.entries(cloneUniformValues(bufferedUniformValues));
  const frontBufferClassicValueEntries: UniformsList = [];
  for (const e of frontBufferUniformValueEntries) {
    if (getBlock(e[0]) === 'object') {
      objectBlockValueEntries.push(e);
    } else {
      frontBufferClassicValueEntries.push(e);
    }
  }
  const defineValueEntries = Object.entries(defineValues);

  // created whenever the program declares the block, members without values stay zero
  let objectUniforms: { buffer: UniformBuffer; write: UniformBlockWriter } | undefined;
  const objectLayout = useUniformBuffers ? getUniformBlockLayout('object', schema) : undefined;
  if (objectLayout) {
    const buffer = resources.uniformBuffer(objectLayout);
    objectUniforms = { buffer, write: buffer.createWriter(objectBlockValueEntries) };
    objectUniforms.write();
  }
  let materialUniforms: { buffer: UniformBuffer; write: UniformBlockWriter } | undefined;
  const materialLayout = useUniformBuffers ? getUniformBlockLayout('material', schema) : undefined;
  if (materialLayout) {
    const buffer = resources.materialUniformBuffer(materialId, materialLayout);
    materialUniforms = { buffer, write: buffer.createWriter(materialBlockValueEntries) };
    materialUniforms.write();
  }

  const glDrawMode = getDrawMode(ctx, drawMode);

  const programs: Programs = {};
  for (const rv of renderVariants) {
    programs[rv] = createProgramVariant(ctx, getRenderVariant(rv, transparency), defineValues, shaderCode, schema);
  }

  // the textures of specs with `shared`, by name, released instead of destroyed (`.v6/plans/gl/shared-textures.md` §4)
  const sharedTextures = new Map<string, SharedTexture>();
  const textures = createTextures(ctx, schema, textureValues, sharedTextures);
  for (const [name, d] of dataTextures) {
    writeData(name, textureValues[name].ref.value as unknown as DataArray<any>);
    textures.push([name, d.range.texture]);
  }
  const materialTextures = createTextures(ctx, schema, materialTextureValues, sharedTextures);
  const getTextureByteCount = (k: string, texture: Texture) =>
    sharedTextures.get(k)?.getByteCount() ?? texture.getByteCount();
  const attributeBuffers = createAttributeBuffers(ctx, schema, attributeValues, contextDefines);
  const instanceBuffers: AttributeBuffers = [];
  // the versions of the values as last seen by `update`, per list in its order: compared for every item on every update,
  // not looked up by name in an object of all the values (a dictionary of a hundred keys)
  const versionsOf = (list: readonly (readonly [string, unknown])[], values: Record<string, ValueCell<any>>) =>
    list.map(([k]) => values[k].ref.version);
  const defineVersions = versionsOf(defineValueEntries, defineValues);
  const attributeVersions = versionsOf(attributeBuffers, attributeValues);
  const textureVersions = versionsOf(textures, textureValues);
  const materialTextureVersions = versionsOf(materialTextures, materialTextureValues);
  const backBufferUniformVersions = versionsOf(backBufferUniformValueEntries, bufferedUniformValues);
  let drawCountVersion = values.drawCount.ref.version;
  let instanceCountVersion = values.instanceCount.ref.version;
  let elementsVersion = values.elements?.ref.version ?? -1;
  for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
    const ab = attributeBuffers[i];
    if (ab[1].divisor === 1) instanceBuffers.push(ab);
  }
  /**
   * per program whether draws without base instances offset its instanced attributes: not where it reads none (the
   * instance data from a texture, `dInstanceData`); the names stay, a buffer replaced by an update keeps its entry
   */
  const offsetsInstances = new WeakMap<Program, boolean>();
  const getOffsetsInstances = (program: Program) => {
    let value = offsetsInstances.get(program);
    if (value === undefined) {
      value = program.readsAnyAttribute(instanceBuffers);
      offsetsInstances.set(program, value);
    }
    return value;
  };

  let elementsBuffer: ElementsBuffer | undefined;
  const elements = values.elements;
  if (elements && elements.ref.value) {
    elementsBuffer = resources.elements(elements.ref.value);
  }
  /** quads without elements of their own draw indexed by the shared quad indices (`quad-indices.ts`) */
  const quads = shaderCode.quads && !elementsBuffer && isQuadIndexing();

  // attribute locations are fixed per schema, so all variants share one vertex array
  const attributeLayout = getAttributeLayout(schema);
  let vertexArray: VertexArray | null = vertexArrayObject
    ? resources.vertexArray(attributeLayout, attributeBuffers, elementsBuffer)
    : null;

  let drawCount: number = values.drawCount.ref.value;
  let instanceCount: number = values.instanceCount.ref.value;

  stats.drawCount += drawCount;
  stats.instanceCount += instanceCount;
  stats.instancedDrawCount += instanceCount * drawCount;

  const valueChanges = createValueChanges();

  let destroyed = false;
  let currentProgramId = -1;

  return {
    id,
    materialId,
    getByteCount() {
      let count = 0;
      for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
        count += attributeBuffers[i][1].getByteCount();
      }
      if (elementsBuffer) count += elementsBuffer.getByteCount();
      if (objectUniforms) count += objectUniforms.buffer.getByteCount();
      // a shared texture by the item's share (§5), a data array by its range of the page
      for (let i = 0, il = textures.length; i < il; ++i) {
        const d = dataTextures.get(textures[i][0]);
        count += d ? d.range.size * d.pages.bytesPerElement : getTextureByteCount(textures[i][0], textures[i][1]);
      }
      for (let i = 0, il = materialTextures.length; i < il; ++i) {
        count += getTextureByteCount(materialTextures[i][0], materialTextures[i][1]);
      }
      return count;
    },
    getProgram: (variant: T) => programs[variant],
    getLinearTextures: () => {
      const names: string[] = [];
      for (const list of [textures, materialTextures]) {
        for (const [name, texture] of list) if ((texture as GLTexture).filter === ctx.gl.LINEAR) names.push(name);
      }
      return names;
    },
    setTransparency: (value: Transparency) => {
      if (value === transparency) return;

      transparency = value;
      for (const rv of renderVariants) {
        programs[rv].destroy();
        programs[rv] = createProgramVariant(ctx, getRenderVariant(rv, transparency), defineValues, shaderCode, schema);
      }
    },

    render: (variant: T, globals?: GlobalBindings, mdbDataList?: MultiDrawBaseData[]) => {
      if (drawCount === 0 || instanceCount === 0) return;
      checkRenderPass(state, 'draw');

      const program = programs[variant];
      if (globals) asWebGLGlobalBindings(globals).apply(program);
      if (program.id === currentProgramId && state.currentRenderItemId === id) {
        program.setUniforms(uniformValueEntries);
        if (objectUniforms) {
          objectUniforms.write();
          objectUniforms.buffer.sync();
          objectUniforms.buffer.bind();
        }
        program.bindTextures(textures);
      } else {
        if (
          program.id !== state.currentProgramId ||
          program.id !== currentProgramId ||
          materialId === -1 ||
          materialId !== state.currentMaterialId
        ) {
          // console.log('program.id changed or materialId changed/-1', materialId)
          if (program.id !== state.currentProgramId) program.use();
          program.setUniforms(materialUniformValueEntries);
          if (materialUniforms) {
            materialUniforms.write();
            materialUniforms.buffer.sync();
            materialUniforms.buffer.bind();
          }
          program.bindTextures(materialTextures);
          state.currentMaterialId = materialId;
          currentProgramId = program.id;
        }
        program.setUniforms(uniformValueEntries);
        program.setUniforms(frontBufferClassicValueEntries);
        if (objectUniforms) {
          objectUniforms.write();
          objectUniforms.buffer.sync();
          objectUniforms.buffer.bind();
        }
        program.bindTextures(textures);
        if (vertexArray) {
          vertexArray.bind();
          // need to bind elements buffer explicitly since it is not always recorded in the VAO
          if (elementsBuffer) elementsBuffer.bind();
        } else {
          if (elementsBuffer) elementsBuffer.bind();
          program.bindAttributes(attributeBuffers);
        }
        state.currentRenderItemId = id;
      }
      // the shared quad indices, bound on every draw (shared, and replaced when grown); unindexed past their limit
      const eb = quads ? resources.quadIndices.get(drawCount / QuadIndicesPerItem)?.buffer : elementsBuffer;
      if (quads && eb) eb.bind();
      if (isDebugMode) {
        ctx.checkFramebufferStatus(`Framebuffer error rendering item id ${id}`);
      }
      if (mdbDataList) {
        // the instanced attributes' offset: 0 where a render begins (each ends with it, a vertex array starts with it),
        // changed only where a draw's first instance differs (one instance per object: never)
        let instanceOffset = 0;
        for (const mdbData of mdbDataList) {
          if (mdbData.count === 0) continue;
          // `indirectInstanceRanges` is false on WebGL: lists with GPU-written arguments are not made for it
          if (mdbData.indirect) throw new Error('indirect draw lists are not supported on WebGL');

          program.setUniforms(mdbData.uniforms);
          // console.log(mdbData.uniforms)
          if (multiDrawInstancedBaseVertexBaseInstance) {
            if (eb) {
              multiDrawInstancedBaseVertexBaseInstance.multiDrawElementsInstancedBaseVertexBaseInstance(
                glDrawMode,
                mdbData.counts,
                0,
                eb._dataType,
                mdbData.offsets,
                0,
                mdbData.instanceCounts,
                0,
                mdbData.baseVertices,
                0,
                mdbData.baseInstances,
                0,
                mdbData.count,
              );
            } else {
              multiDrawInstancedBaseVertexBaseInstance.multiDrawArraysInstancedBaseInstance(
                glDrawMode,
                mdbData.firsts,
                0,
                mdbData.counts,
                0,
                mdbData.instanceCounts,
                0,
                mdbData.baseInstances,
                0,
                mdbData.count,
              );
            }
          } else if (drawInstancedBaseVertexBaseInstance) {
            if (eb) {
              for (let i = 0; i < mdbData.count; ++i) {
                if (mdbData.counts[i] > 0) {
                  program.uniform('uDrawId', i);
                  drawInstancedBaseVertexBaseInstance.drawElementsInstancedBaseVertexBaseInstance(
                    glDrawMode,
                    mdbData.counts[i],
                    eb._dataType,
                    mdbData.offsets[i],
                    mdbData.instanceCounts[i],
                    mdbData.baseVertices[i],
                    mdbData.baseInstances[i],
                  );
                }
              }
            } else {
              for (let i = 0; i < mdbData.count; ++i) {
                if (mdbData.counts[i] > 0) {
                  program.uniform('uDrawId', i);
                  drawInstancedBaseVertexBaseInstance.drawArraysInstancedBaseInstance(
                    glDrawMode,
                    mdbData.firsts[i],
                    mdbData.counts[i],
                    mdbData.instanceCounts[i],
                    mdbData.baseInstances[i],
                  );
                }
              }
            }
          } else {
            const offsets = getOffsetsInstances(program);
            if (eb) {
              for (let i = 0; i < mdbData.count; ++i) {
                if (mdbData.counts[i] > 0) {
                  program.uniform('uDrawId', i);
                  if (instanceData) program.uniform('uInstanceBase', mdbData.baseInstances[i]);
                  if (offsets && mdbData.baseInstances[i] !== instanceOffset) {
                    instanceOffset = mdbData.baseInstances[i];
                    program.offsetAttributes(instanceBuffers, instanceOffset);
                  }
                  instancedArrays.drawElementsInstanced(
                    glDrawMode,
                    mdbData.counts[i],
                    eb._dataType,
                    mdbData.offsets[i],
                    mdbData.instanceCounts[i],
                  );
                }
              }
            } else {
              for (let i = 0; i < mdbData.count; ++i) {
                if (mdbData.counts[i] > 0) {
                  program.uniform('uDrawId', i);
                  if (instanceData) program.uniform('uInstanceBase', mdbData.baseInstances[i]);
                  if (offsets && mdbData.baseInstances[i] !== instanceOffset) {
                    instanceOffset = mdbData.baseInstances[i];
                    program.offsetAttributes(instanceBuffers, instanceOffset);
                  }
                  instancedArrays.drawArraysInstanced(
                    glDrawMode,
                    mdbData.firsts[i],
                    mdbData.counts[i],
                    mdbData.instanceCounts[i],
                  );
                }
              }
            }
          }
          if (isTimingMode) {
            if (multiDrawInstancedBaseVertexBaseInstance) {
              stats.calls.multiDrawInstancedBase += 1;
            } else if (drawInstancedBaseVertexBaseInstance) {
              stats.calls.drawInstancedBase += mdbData.count;
            } else {
              stats.calls.drawInstanced += mdbData.count;
            }
            for (let i = 0; i < mdbData.count; ++i) {
              stats.calls.counts += mdbData.instanceCounts[i];
            }
          }
        }
        if (!multiDrawInstancedBaseVertexBaseInstance && !drawInstancedBaseVertexBaseInstance) {
          // offsets are recorded in the vertex array, later draws without base instances expect 0
          if (instanceOffset !== 0) program.offsetAttributes(instanceBuffers, 0);
        }
      } else {
        let offset = 0;
        while (true) {
          const count = Math.min(drawCount - offset, MaxDrawCount);
          if (eb) {
            instancedArrays.drawElementsInstanced(glDrawMode, count, eb._dataType, offset * eb._bpe, instanceCount);
          } else {
            instancedArrays.drawArraysInstanced(glDrawMode, offset, count, instanceCount);
          }
          offset += count;
          if (offset >= drawCount) break;
        }
        if (isTimingMode) {
          stats.calls.drawInstanced += 1;
          stats.calls.counts += instanceCount;
        }
      }
      if (isDebugMode) {
        ctx.checkError(`Draw error rendering item id ${id}`);
      }
    },
    update: () => {
      resetValueChanges(valueChanges);
      // values may have changed that the fast path for repeated draws doesn't rebind
      // (attributes, elements, textures, front-buffer and material uniforms)
      currentProgramId = -1;

      if (emulateVertexId) ensureVertexIds(getEmulatedVertexCount());
      if (instanceData) syncInstanceData();

      for (let i = 0, il = defineValueEntries.length; i < il; ++i) {
        const value = defineValueEntries[i][1];
        if (value.ref.version !== defineVersions[i]) {
          valueChanges.defines = true;
          defineVersions[i] = value.ref.version;
        }
      }

      if (valueChanges.defines) {
        // console.log('some defines changed, need to rebuild programs');
        for (const rv of renderVariants) {
          programs[rv].destroy();
          programs[rv] = createProgramVariant(
            ctx,
            getRenderVariant(rv, transparency),
            defineValues,
            shaderCode,
            schema,
          );
        }
      }

      if (values.drawCount.ref.version !== drawCountVersion) {
        // console.log('drawCount version changed');
        stats.drawCount += values.drawCount.ref.value - drawCount;
        stats.instancedDrawCount += instanceCount * values.drawCount.ref.value - instanceCount * drawCount;
        drawCount = values.drawCount.ref.value;
        drawCountVersion = values.drawCount.ref.version;
      }
      if (values.instanceCount.ref.version !== instanceCountVersion) {
        // console.log('instanceCount version changed');
        stats.instanceCount += values.instanceCount.ref.value - instanceCount;
        stats.instancedDrawCount += values.instanceCount.ref.value * drawCount - instanceCount * drawCount;
        instanceCount = values.instanceCount.ref.value;
        instanceCountVersion = values.instanceCount.ref.version;
      }

      for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
        const [k, buffer] = attributeBuffers[i];
        const value = attributeValues[k];
        if (value.ref.version !== attributeVersions[i]) {
          if (buffer.length >= value.ref.value.length) {
            // console.log('attribute array large enough to update', buffer.id, k, value.ref.id, value.ref.version);
            buffer.updateSubData(value.ref.value, 0, buffer.length);
          } else {
            // console.log('attribute array too small, need to create new attribute', buffer.id, k, value.ref.id, value.ref.version);
            buffer.destroy();
            const { itemSize, divisor } = schema[k] as AttributeSpec<AttributeKind>;
            attributeBuffers[i][1] = resources.attribute(value.ref.value, itemSize, divisor);
            valueChanges.attributes = true;
          }
          attributeVersions[i] = value.ref.version;
        }
      }

      if (elementsBuffer && values.elements.ref.version !== elementsVersion) {
        if (elementsBuffer.length >= values.elements.ref.value.length) {
          // console.log('elements array large enough to update', values.elements.ref.id, values.elements.ref.version);
          elementsBuffer.updateSubData(values.elements.ref.value, 0, elementsBuffer.length);
        } else {
          // console.log('elements array to small, need to create new elements', values.elements.ref.id, values.elements.ref.version);
          elementsBuffer.destroy();
          elementsBuffer = resources.elements(values.elements.ref.value);
          valueChanges.elements = true;
        }
        elementsVersion = values.elements.ref.version;
      }

      if (valueChanges.attributes || valueChanges.elements) {
        // console.log('buffers changed, update vao');
        if (vertexArray) vertexArray.destroy();
        vertexArray = vertexArrayObject
          ? resources.vertexArray(attributeLayout, attributeBuffers, elementsBuffer)
          : null;
      }

      for (let i = 0, il = textures.length; i < il; ++i) {
        const [k, texture] = textures[i];
        const value = textureValues[k];
        if (value.ref.version !== textureVersions[i]) {
          const held = sharedTextures.get(k);
          if (held) {
            // the same image uploads once for all its items, another image is another texture (§4)
            const next = updateSharedTexture(
              resources.sharedTextures,
              held,
              schema[k] as TextureSpec<TextureKind>,
              value.ref.value as SharedTextureImage,
            );
            if (next !== held) {
              sharedTextures.set(k, next);
              textures[i][1] = next.texture;
            }
            valueChanges.textures = true;
          } else if ((schema[k] as TextureSpec<TextureKind>).data !== undefined) {
            // a data array (`DataSpec`): the changed elements only where logged (`TextureChanges`)
            const data = value.ref.value as unknown as DataArray<any>;
            writeData(k, data, TextureChanges.get(data, textureVersions[i], value.ref.version));
            textures[i][1] = dataTextures.get(k)!.range.texture;
            valueChanges.textures = true;
          } else if (schema[k].kind !== 'texture') {
            // update of textures with kind 'texture' is done externally
            texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
            valueChanges.textures = true;
          } else {
            textures[i][1] = value.ref.value as Texture;
          }
          textureVersions[i] = value.ref.version;
        }
      }

      for (let i = 0, il = materialTextures.length; i < il; ++i) {
        const [k, texture] = materialTextures[i];
        const value = materialTextureValues[k];
        if (value.ref.version !== materialTextureVersions[i]) {
          const held = sharedTextures.get(k);
          if (held) {
            // the same image uploads once for all its items, another image is another texture (§4)
            const next = updateSharedTexture(
              resources.sharedTextures,
              held,
              schema[k] as TextureSpec<TextureKind>,
              value.ref.value as SharedTextureImage,
            );
            if (next !== held) {
              sharedTextures.set(k, next);
              materialTextures[i][1] = next.texture;
            }
            valueChanges.textures = true;
          } else if (schema[k].kind !== 'texture') {
            // update of textures with kind 'texture' is done externally
            // console.log('materialTexture version changed, uploading image', k);
            texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
            valueChanges.textures = true;
          } else {
            materialTextures[i][1] = value.ref.value as Texture;
          }
          materialTextureVersions[i] = value.ref.version;
        }
      }

      for (let i = 0, il = backBufferUniformValueEntries.length; i < il; ++i) {
        const uniform = backBufferUniformValueEntries[i][1];
        if (uniform.ref.version !== backBufferUniformVersions[i]) {
          ValueCell.update(frontBufferUniformValueEntries[i][1], deepClone(uniform.ref.value));
          backBufferUniformVersions[i] = uniform.ref.version;
        }
      }
    },
    destroy: () => {
      if (!destroyed) {
        for (const k of renderVariants) {
          programs[k].destroy();
        }
        if (vertexArray) vertexArray.destroy();
        for (const list of [textures, materialTextures]) {
          for (const [k, texture] of list) {
            const held = sharedTextures.get(k);
            const d = dataTextures.get(k);
            if (held) held.release();
            else if (d) d.pages.free(d.range);
            // lifetime of textures with kind 'texture' is defined externally
            else if (schema[k].kind !== 'texture') texture.destroy();
          }
        }
        sharedTextures.clear();
        attributeBuffers.forEach(([_, buffer]) => buffer.destroy());
        if (elementsBuffer) elementsBuffer.destroy();
        if (objectUniforms) {
          objectUniforms.write.dispose();
          objectUniforms.buffer.destroy();
        }
        if (materialUniforms) {
          materialUniforms.write.dispose();
          materialUniforms.buffer.destroy();
        }

        stats.drawCount -= drawCount;
        stats.instanceCount -= instanceCount;
        stats.instancedDrawCount -= instanceCount * drawCount;

        destroyed = true;
      }
    },
  };
}
