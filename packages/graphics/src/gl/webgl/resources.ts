/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject } from 'rxjs';
import type { SegmentedArray } from '../segmented-array.js';
import { type ProgramProps, createProgram, type Program } from './program.js';
import { type ShaderType, createShader, type Shader, type ShaderProps } from './shader.js';
import { type GLRenderingContext, isWebGL2 } from './compat.js';
import { createFramebuffer } from './framebuffer.js';
import type { WebGLExtensions } from './extensions.js';
import type { WebGLState } from './state.js';
import {
  type AttributeBuffer,
  createAttributeBuffer,
  type ElementsBuffer,
  createElementsBuffer,
  type AttributeBuffers,
  type PixelPackBuffer,
  createPixelPackBuffer,
} from './buffer.js';
import type { ArrayType, AttributeItemSize, ElementsType, UsageHint } from '../buffer.js';
import { createReferenceCache, type ReferenceItem } from '@molstar/core/util/reference-cache';
import type { WebGLParameters } from './context.js';
import type { GraphicsStats } from '../stats.js';
import { hashString, hashFnv32a } from '@molstar/core/data/util';
import type { DefineValues, ShaderCode } from '../shader-code.js';
import type { RenderableSchema } from '../renderable/schema.js';
import { createRenderbuffer } from './renderbuffer.js';
import type { Renderbuffer, RenderbufferAttachment, RenderbufferFormat } from '../renderbuffer.js';
import { createTexture, createCubeTexture, type GLTexture, registerGLTexture } from './texture.js';
import {
  type CubeFaces,
  decodeImage,
  type Texture,
  type TextureFilter,
  type TextureFormat,
  type TextureKind,
  type TextureType,
} from '../texture.js';
import { ValueCell } from '@molstar/core/util';
import { createPendingWork } from '../pending-work.js';
import { type VertexArray, createVertexArray } from './vertex-array.js';
import {
  createUniformArena,
  createUniformBuffer,
  createUniformRingBuffer,
  type UniformArena,
} from './uniform-buffer.js';
import type { UniformBuffer } from '../uniform.js';
import type { UniformBlockLayout } from '../uniform-layout.js';
import type { AttributeLayout } from '../attribute-layout.js';
import { getShaderDeclarations } from './shader-declarations.js';
import { createSharedTextures, type SharedTextures } from '../shared-texture.js';
import { createQuadIndices, getMaxQuadItemCount, type QuadIndices } from '../quad-indices.js';
import { getDefineLookup, getDefineValueNames } from '../renderable/schema.js';
import { now } from '@molstar/core/util/now';
import type { ProgramVariant, RenderItem } from '../render-item.js';
import { isTimingMode } from '@molstar/core/util/debug';
import type { DataFormat } from '../renderable/schema.js';
import { createDataTexturePages, type DataTexturePages } from './data-pages.js';
import type { GraphicsLimits, GraphicsResources, PipelineRequest } from '../context.js';

function defineValueHash(v: boolean | number | string): number {
  return typeof v === 'boolean' ? (v ? 1 : 0) : typeof v === 'number' ? v * 10000 : hashString(v);
}

/** The largest shared quad index buffer (`quad-indices.ts`): WebGL buffers are sized by the driver */
const QuadIndicesMaxBytes = 2 ** 31;

function wrapCached<T extends Resource>(resourceItem: ReferenceItem<T>) {
  const wrapped = {
    ...resourceItem.value,
    destroy: () => {
      resourceItem.free();
    },
  };

  return wrapped;
}

//

interface Resource {
  destroy: () => void;
}

type ResourceName = keyof GraphicsStats['resourceCounts'];

export interface WebGLResources extends GraphicsResources {
  attribute: (
    array: ArrayType | SegmentedArray,
    itemSize: AttributeItemSize,
    divisor: number,
    usageHint?: UsageHint,
  ) => AttributeBuffer;
  elements: (array: ElementsType | SegmentedArray, usageHint?: UsageHint) => ElementsBuffer;
  pixelPack: (format: TextureFormat, type: TextureType) => PixelPackBuffer;
  program: (defineValues: DefineValues, shaderCode: ShaderCode, schema: RenderableSchema) => Program;
  shader: (type: ShaderType, source: string) => Shader;
  texture: (kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter) => GLTexture;
  cubeTexture: (faces: CubeFaces, mipmaps: boolean, onload?: (errored?: boolean) => void) => GLTexture;
  /** Textures shared by render items, from `texture` (`.v6/plans/gl/shared-textures.md` §2) */
  sharedTextures: SharedTextures;
  /** The index buffer of quad items (spheres, points), shared by all (`quad-indices.ts`) */
  quadIndices: QuadIndices<ElementsBuffer>;
  /** The pages of the data arrays (`DataSpec`) of a format, shared by render items (`data-pages.ts`) */
  dataPages: (format: DataFormat) => DataTexturePages;
  vertexArray: (
    layout: AttributeLayout,
    attributeBuffers: AttributeBuffers,
    elementsBuffer?: ElementsBuffer,
  ) => VertexArray;

  destroy: () => void;
}

export function createResources(
  gl: GLRenderingContext,
  state: WebGLState,
  stats: GraphicsStats,
  extensions: WebGLExtensions,
  parameters: WebGLParameters,
  uniformBlocks: boolean,
  /** whether the context is lost: its objects are gone, also once the canvas' context is restored (on the same `gl`) */
  isLost: () => boolean = () => false,
): WebGLResources {
  const sets: { [k in ResourceName]: Set<Resource> } = {
    attribute: new Set<Resource>(),
    elements: new Set<Resource>(),
    pixelPack: new Set<Resource>(),
    uniform: new Set<Resource>(),
    framebuffer: new Set<Resource>(),
    program: new Set<Resource>(),
    renderbuffer: new Set<Resource>(),
    shader: new Set<Resource>(),
    texture: new Set<Resource>(),
    cubeTexture: new Set<Resource>(),
    vertexArray: new Set<Resource>(),
    // compute resources, WebGPU only
    storageBuffer: new Set<Resource>(),
    kernel: new Set<Resource>(),
  };

  function wrap<T extends Resource>(name: ResourceName, resource: T) {
    sets[name].add(resource);
    stats.resourceCounts[name] += 1;
    return {
      ...resource,
      destroy: () => {
        if (!sets[name].has(resource)) return;
        // deleting objects of a lost context through the restored `gl` is an error (INVALID_OPERATION)
        if (!isLost()) resource.destroy();
        sets[name].delete(resource);
        stats.resourceCounts[name] -= 1;
      },
    };
  }

  const shaderCache = createReferenceCache(
    (props: ShaderProps) => JSON.stringify(props),
    (props: ShaderProps) => wrap('shader', createShader(gl, props)),
    (shader: Shader) => {
      shader.destroy();
    },
  );

  function getShader(type: ShaderType, source: string) {
    return wrapCached(shaderCache.get({ type, source }));
  }

  const ownTextures = new WeakSet<Texture>();
  const pendingPrograms = new Set<Program>();
  const pending = createPendingWork();
  const outOfMemory = new Subject<string>();

  /** programs linking in parallel (`KHR_parallel_shader_compile`), pending work until they finished */
  const compilingPrograms = new Set<Program>();
  let compiled: (() => void) | undefined;
  function pollCompiling() {
    for (const p of compilingPrograms) if (p.isCompiled()) compilingPrograms.delete(p);
    if (compilingPrograms.size > 0) {
      setTimeout(pollCompiling, 4);
    } else {
      compiled?.();
      compiled = undefined;
    }
  }
  function trackCompiling(program: Program) {
    if (!extensions.parallelShaderCompile || program.isCompiled()) return;
    compilingPrograms.add(program);
    if (compiled) return;
    pending.track(new Promise<void>((resolve) => (compiled = resolve)), 'programs');
    setTimeout(pollCompiling, 0);
  }
  function clearCompiling() {
    compilingPrograms.clear();
    compiled?.();
    compiled = undefined;
  }
  /**
   * Hashes of the generated declarations per schema and defined defines: generating them is most of a program lookup,
   * which render items do per variant (cache hits too)
   */
  const declarationHashes = new WeakMap<RenderableSchema, Map<string, readonly [number, number]>>();
  function getDeclarationHashes(schema: RenderableSchema, defineValues: DefineValues) {
    let bySchema = declarationHashes.get(schema);
    if (!bySchema) declarationHashes.set(schema, (bySchema = new Map()));
    const names = getDefineValueNames(defineValues);
    const key = names.join(' ');
    let hashes = bySchema.get(key);
    if (!hashes) {
      const declarations = getShaderDeclarations(schema, isWebGL2(gl), getDefineLookup(names), uniformBlocks);
      hashes = [hashString(declarations.vert), hashString(declarations.frag)];
      bySchema.set(key, hashes);
    }
    return hashes;
  }

  const programCache = createReferenceCache(
    (props: ProgramProps) => {
      const array = [props.shaderCode.id];
      const variant = (props.defineValues.dRenderVariant?.ref.value || '') as string;
      Object.keys(props.defineValues).forEach((k) => {
        if (!props.shaderCode.ignoreDefine?.(k, variant, props.defineValues)) {
          array.push(hashString(k), defineValueHash(props.defineValues[k].ref.value));
        }
      });
      // generated declarations (incl. attribute locations) depend on the schema and the defines
      const [vert, frag] = getDeclarationHashes(props.schema, props.defineValues);
      array.push(vert, frag);
      return hashFnv32a(array).toString();
    },
    (props: ProgramProps) => {
      const program = createProgram(gl, state, extensions, parameters, getShader, props, uniformBlocks);
      if (program.variant !== 'compute') {
        pendingPrograms.add(program);
      }
      return wrap('program', program);
    },
    (program: Program) => {
      pendingPrograms.delete(program);
      program.destroy();
    },
  );

  let uniformArena: UniformArena | undefined;
  function getUniformArena() {
    if (!isWebGL2(gl)) {
      throw new Error('WebGL2 is required for uniform buffers');
    }
    if (!uniformArena) uniformArena = createUniformArena(gl, state);
    return uniformArena;
  }

  function uniformBuffer(layout: UniformBlockLayout) {
    return wrap('uniform', createUniformBuffer(getUniformArena(), layout));
  }

  const materialUniformBufferCache = createReferenceCache(
    (props: { materialId: number; layout: UniformBlockLayout }) => `${props.materialId}|${props.layout.hash}`,
    (props: { materialId: number; layout: UniformBlockLayout }) => uniformBuffer(props.layout),
    (buffer: UniformBuffer) => {
      buffer.destroy();
    },
  );

  const limits: GraphicsLimits = {
    maxTextureSize: parameters.maxTextureSize,
    max3dTextureSize: parameters.max3dTextureSize,
    maxRenderbufferSize: parameters.maxRenderbufferSize,
    maxDrawBuffers: parameters.maxDrawBuffers,
    maxTextureImageUnits: parameters.maxTextureImageUnits,
    maxQuadItemCount: getMaxQuadItemCount(QuadIndicesMaxBytes),
  };

  function texture(kind: TextureKind, format: TextureFormat, type: TextureType, filter: TextureFilter) {
    const texture = registerGLTexture(
      wrap('texture', createTexture(gl, extensions, kind, format, type, filter, state)),
    );
    ownTextures.add(texture);
    return texture;
  }

  // gone with the resources: a restored context gets new ones (§6)
  const sharedTextures = createSharedTextures((spec) => texture(spec.kind, spec.format, spec.dataType, spec.filter));
  // counted with the elements buffers
  const quadIndices = createQuadIndices(
    (indices) => wrap('elements', createElementsBuffer(gl, indices)),
    (buffer) => buffer.destroy(),
    QuadIndicesMaxBytes,
  );

  // gone with the resources, their textures with the others
  const dataPages: { [F in DataFormat]?: DataTexturePages } = {};
  const createDataPageTexture = (format: DataFormat) =>
    format === 'rgba32f'
      ? texture('image-float32', 'rgba', 'float', 'nearest')
      : texture('image-uint8', format === 'r8' ? 'r' : 'rgba', 'ubyte', 'nearest');

  return {
    limits,
    dataPages: (format: DataFormat) =>
      (dataPages[format] ??= createDataTexturePages(format, createDataPageTexture, parameters.maxTextureSize)),
    attribute: (
      array: ArrayType | SegmentedArray,
      itemSize: AttributeItemSize,
      divisor: number,
      usageHint?: UsageHint,
    ) => {
      return wrap('attribute', createAttributeBuffer(gl, state, extensions, array, itemSize, divisor, usageHint));
    },
    elements: (array: ElementsType | SegmentedArray, usageHint?: UsageHint) => {
      return wrap('elements', createElementsBuffer(gl, array, usageHint));
    },
    pixelPack: (format: TextureFormat, type: TextureType) => {
      if (!isWebGL2(gl)) {
        throw new Error('WebGL2 is required for pixel-pack buffers');
      }
      return wrap('pixelPack', createPixelPackBuffer(gl, extensions, format, type));
    },
    uniformBuffer,
    uniformRingBuffer: (layout: UniformBlockLayout) => {
      return wrap('uniform', createUniformRingBuffer(getUniformArena(), layout));
    },
    materialUniformBuffer: (materialId: number, layout: UniformBlockLayout) => {
      // -1 means unique material, no sharing possible
      if (materialId === -1) return uniformBuffer(layout);
      return wrapCached(materialUniformBufferCache.get({ materialId, layout }));
    },
    framebuffer: () => {
      return wrap('framebuffer', createFramebuffer(gl, extensions, state));
    },
    program: (defineValues: DefineValues, shaderCode: ShaderCode, schema: RenderableSchema) => {
      return wrapCached(programCache.get({ defineValues, shaderCode, schema }));
    },
    renderbuffer: (format: RenderbufferFormat, attachment: RenderbufferAttachment, width: number, height: number) => {
      return wrap('renderbuffer', createRenderbuffer(gl, format, attachment, width, height, state));
    },
    shader: getShader,
    texture,
    sharedTextures,
    quadIndices,
    ownsTexture: (texture: Texture) => ownTextures.has(texture),
    cubeTexture: (faces: CubeFaces, mipmaps: boolean, onload?: (errored?: boolean) => void) => {
      // settles once: a face that fails to load reports an error of its own
      let settle: (() => void) | undefined;
      pending.track(new Promise<void>((resolve) => (settle = resolve)), 'cube texture');
      // `onload` before settling, so that work waiting for the faces sees what `onload` did
      const texture = createCubeTexture(gl, faces, mipmaps, (errored) => {
        onload?.(errored);
        settle?.();
        settle = undefined;
      });
      return registerGLTexture(wrap('cubeTexture', texture));
    },
    vertexArray: (layout: AttributeLayout, attributeBuffers: AttributeBuffers, elementsBuffer?: ElementsBuffer) => {
      return wrap('vertexArray', createVertexArray(gl, state, extensions, layout, attributeBuffers, elementsBuffer));
    },

    getByteCounts: () => {
      let texture = 0;
      sets.texture.forEach((r) => {
        texture += (r as Texture).getByteCount();
      });

      let cubeTexture = 0;
      sets.cubeTexture.forEach((r) => {
        cubeTexture += (r as Texture).getByteCount();
      });

      let attribute = 0;
      sets.attribute.forEach((r) => {
        attribute += (r as AttributeBuffer).getByteCount();
      });

      let elements = 0;
      sets.elements.forEach((r) => {
        elements += (r as ElementsBuffer).getByteCount();
      });

      let pixelPack = 0;
      sets.pixelPack.forEach((r) => {
        pixelPack += (r as PixelPackBuffer).getByteCount();
      });

      const uniform = uniformArena?.getByteCount() ?? 0;

      let renderbuffer = 0;
      sets.renderbuffer.forEach((r) => {
        renderbuffer += (r as Renderbuffer).getByteCount();
      });

      // no storage buffers on WebGL
      return { texture, cubeTexture, attribute, elements, pixelPack, uniform, renderbuffer, storageBuffer: 0 };
    },

    pending,
    outOfMemory,
    loadImage: (texture: Texture, source: string, cell?: ValueCell<Texture>, label = 'image') => {
      const load = decodeImage(source).then((image) => {
        // uploaded with `UNPACK_COLORSPACE_CONVERSION_WEBGL` `NONE`
        texture.load(image);
        if (cell) ValueCell.update(cell, texture);
      });
      return pending.track(load, label);
    },

    linkPrograms: (variants?: ProgramVariant[]) => {
      for (const p of pendingPrograms) {
        if (variants && !variants.includes(p.variant)) continue;
        p.link();
        trackCompiling(p);
      }
    },
    preparePipelines: (requests: Iterable<PipelineRequest>) => {
      for (const { item, variant } of requests) {
        const program = (item as RenderItem<ProgramVariant>).getProgram(variant) as Program;
        program.link();
        trackCompiling(program);
      }
      return pending.whenReady();
    },
    finalizePrograms: (variants?: ProgramVariant[], isSynchronous?: boolean) => {
      let isReady = true;
      let pendingCount = 0;
      for (const p of pendingPrograms) {
        if (p.isReady()) pendingPrograms.delete(p);
        if (!variants || variants.includes(p.variant)) {
          isReady = false;
          pendingCount += 1;
        }
      }
      if (isReady) return true;

      let linkStatus = true;
      let finalizedCount = 0;
      const t = now();
      for (const p of pendingPrograms) {
        if (variants && !variants.includes(p.variant)) continue;
        if (!p.finalize(isSynchronous)) {
          linkStatus = false;
        } else {
          pendingPrograms.delete(p);
          finalizedCount += 1;
        }
        if (!isSynchronous && now() - t > 16) {
          linkStatus = false;
          break;
        }
      }

      if (isTimingMode) {
        console.log(
          `Finalized ${finalizedCount} of ${pendingCount} programs (${variants ? variants.join(', ') : 'all'}) in ${(now() - t).toFixed(2)} ms`,
        );
      }

      return linkStatus;
    },

    destroy: () => {
      clearCompiling();
      pending.clear();
      sets.attribute.forEach((r) => r.destroy());
      sets.elements.forEach((r) => r.destroy());
      sets.pixelPack.forEach((r) => r.destroy());
      sets.uniform.forEach((r) => r.destroy());
      sets.framebuffer.forEach((r) => r.destroy());
      sets.renderbuffer.forEach((r) => r.destroy());
      sets.shader.forEach((r) => r.destroy());
      sets.program.forEach((r) => r.destroy());
      sets.vertexArray.forEach((r) => r.destroy());
      sets.texture.forEach((r) => r.destroy());
      sets.cubeTexture.forEach((r) => r.destroy());

      shaderCache.clear();
      programCache.clear();
      materialUniformBufferCache.clear();
      uniformArena?.destroy();
      pendingPrograms.clear();
    },
  };
}
