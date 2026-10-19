/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { UniformBlockWriter } from '../uniform-layout.js';
import { createAttributeBuffers, type ElementsBuffer, type AttributeBuffers } from './buffer.js';
import type { AttributeKind } from '../buffer.js';
import { createTextures, type GLTexture } from './texture.js';
import type { Texture } from '../texture.js';
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
  getValueVersions,
  isUsed,
  splitValues,
  DefineSpec,
  type UniformSpec,
} from '../renderable/schema.js';
import { idFactory } from '@molstar/core/util/id-factory';
import { ValueCell } from '@molstar/core/util';
import type { TextureImage, TextureVolume } from '../renderable/util.js';
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
  const getProgram = (variant: GraphicsRenderVariant) => {
    let program = programs.get(variant);
    if (!program) {
      program = createProgramVariant(ctx, getRenderVariant(variant, transparency), defineValues, shaderCode, schema);
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

  const versions = getValueVersions(values);

  const glDrawMode = getDrawMode(ctx, drawMode);

  const programs: Programs = {};
  for (const rv of renderVariants) {
    programs[rv] = createProgramVariant(ctx, getRenderVariant(rv, transparency), defineValues, shaderCode, schema);
  }

  const textures = createTextures(ctx, schema, textureValues);
  const materialTextures = createTextures(ctx, schema, materialTextureValues);
  const attributeBuffers = createAttributeBuffers(ctx, schema, attributeValues, contextDefines);
  const instanceBuffers: AttributeBuffers = [];
  for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
    const ab = attributeBuffers[i];
    if (ab[1].divisor === 1) instanceBuffers.push(ab);
  }

  let elementsBuffer: ElementsBuffer | undefined;
  const elements = values.elements;
  if (elements && elements.ref.value) {
    elementsBuffer = resources.elements(elements.ref.value);
  }

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
      for (let i = 0, il = textures.length; i < il; ++i) {
        count += textures[i][1].getByteCount();
      }
      for (let i = 0, il = materialTextures.length; i < il; ++i) {
        count += materialTextures[i][1].getByteCount();
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
      if (isDebugMode) {
        ctx.checkFramebufferStatus(`Framebuffer error rendering item id ${id}`);
      }
      if (mdbDataList) {
        for (const mdbData of mdbDataList) {
          if (mdbData.count === 0) continue;

          program.setUniforms(mdbData.uniforms);
          // console.log(mdbData.uniforms)
          if (multiDrawInstancedBaseVertexBaseInstance) {
            if (elementsBuffer) {
              multiDrawInstancedBaseVertexBaseInstance.multiDrawElementsInstancedBaseVertexBaseInstance(
                glDrawMode,
                mdbData.counts,
                0,
                elementsBuffer._dataType,
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
            if (elementsBuffer) {
              for (let i = 0; i < mdbData.count; ++i) {
                if (mdbData.counts[i] > 0) {
                  program.uniform('uDrawId', i);
                  drawInstancedBaseVertexBaseInstance.drawElementsInstancedBaseVertexBaseInstance(
                    glDrawMode,
                    mdbData.counts[i],
                    elementsBuffer._dataType,
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
            if (elementsBuffer) {
              for (let i = 0; i < mdbData.count; ++i) {
                if (mdbData.counts[i] > 0) {
                  program.uniform('uDrawId', i);
                  program.offsetAttributes(instanceBuffers, mdbData.baseInstances[i]);
                  instancedArrays.drawElementsInstanced(
                    glDrawMode,
                    mdbData.counts[i],
                    elementsBuffer._dataType,
                    mdbData.offsets[i],
                    mdbData.instanceCounts[i],
                  );
                }
              }
            } else {
              for (let i = 0; i < mdbData.count; ++i) {
                if (mdbData.counts[i] > 0) {
                  program.uniform('uDrawId', i);
                  program.offsetAttributes(instanceBuffers, mdbData.baseInstances[i]);
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
          program.offsetAttributes(instanceBuffers, 0);
        }
      } else {
        let offset = 0;
        while (true) {
          const count = Math.min(drawCount - offset, MaxDrawCount);
          if (elementsBuffer) {
            instancedArrays.drawElementsInstanced(
              glDrawMode,
              count,
              elementsBuffer._dataType,
              offset * elementsBuffer._bpe,
              instanceCount,
            );
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

      for (let i = 0, il = defineValueEntries.length; i < il; ++i) {
        const [k, value] = defineValueEntries[i];
        if (value.ref.version !== versions[k]) {
          // console.log('define version changed', k);
          valueChanges.defines = true;
          versions[k] = value.ref.version;
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

      if (values.drawCount.ref.version !== versions.drawCount) {
        // console.log('drawCount version changed');
        stats.drawCount += values.drawCount.ref.value - drawCount;
        stats.instancedDrawCount += instanceCount * values.drawCount.ref.value - instanceCount * drawCount;
        drawCount = values.drawCount.ref.value;
        versions.drawCount = values.drawCount.ref.version;
      }
      if (values.instanceCount.ref.version !== versions.instanceCount) {
        // console.log('instanceCount version changed');
        stats.instanceCount += values.instanceCount.ref.value - instanceCount;
        stats.instancedDrawCount += values.instanceCount.ref.value * drawCount - instanceCount * drawCount;
        instanceCount = values.instanceCount.ref.value;
        versions.instanceCount = values.instanceCount.ref.version;
      }

      for (let i = 0, il = attributeBuffers.length; i < il; ++i) {
        const [k, buffer] = attributeBuffers[i];
        const value = attributeValues[k];
        if (value.ref.version !== versions[k]) {
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
          versions[k] = value.ref.version;
        }
      }

      if (elementsBuffer && values.elements.ref.version !== versions.elements) {
        if (elementsBuffer.length >= values.elements.ref.value.length) {
          // console.log('elements array large enough to update', values.elements.ref.id, values.elements.ref.version);
          elementsBuffer.updateSubData(values.elements.ref.value, 0, elementsBuffer.length);
        } else {
          // console.log('elements array to small, need to create new elements', values.elements.ref.id, values.elements.ref.version);
          elementsBuffer.destroy();
          elementsBuffer = resources.elements(values.elements.ref.value);
          valueChanges.elements = true;
        }
        versions.elements = values.elements.ref.version;
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
        if (value.ref.version !== versions[k]) {
          // update of textures with kind 'texture' is done externally
          if (schema[k].kind !== 'texture') {
            // console.log('texture version changed, uploading image', k);
            texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
            valueChanges.textures = true;
          } else {
            textures[i][1] = value.ref.value as Texture;
          }
          versions[k] = value.ref.version;
        }
      }

      for (let i = 0, il = materialTextures.length; i < il; ++i) {
        const [k, texture] = materialTextures[i];
        const value = materialTextureValues[k];
        if (value.ref.version !== versions[k]) {
          // update of textures with kind 'texture' is done externally
          if (schema[k].kind !== 'texture') {
            // console.log('materialTexture version changed, uploading image', k);
            texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
            valueChanges.textures = true;
          } else {
            materialTextures[i][1] = value.ref.value as Texture;
          }
          versions[k] = value.ref.version;
        }
      }

      for (let i = 0, il = backBufferUniformValueEntries.length; i < il; ++i) {
        const [k, uniform] = backBufferUniformValueEntries[i];
        if (uniform.ref.version !== versions[k]) {
          // console.log('back-buffer uniform version changed, updating front-buffer', k);
          ValueCell.update(frontBufferUniformValueEntries[i][1], deepClone(uniform.ref.value));
          versions[k] = uniform.ref.version;
        }
      }
    },
    destroy: () => {
      if (!destroyed) {
        for (const k of renderVariants) {
          programs[k].destroy();
        }
        if (vertexArray) vertexArray.destroy();
        textures.forEach(([k, texture]) => {
          // lifetime of textures with kind 'texture' is defined externally
          if (schema[k].kind !== 'texture') texture.destroy();
        });
        materialTextures.forEach(([k, texture]) => {
          // lifetime of textures with kind 'texture' is defined externally
          if (schema[k].kind !== 'texture') texture.destroy();
        });
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
