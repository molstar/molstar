/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ValueCell } from '@molstar/core/util';
import { idFactory } from '@molstar/core/util/id-factory';
import type { DefineSpec, KindValue, TextureSpec, UniformSpec, UsedIf } from './renderable/schema.js';
import type { DefineKind } from './shader-code.js';
import type { Texture, TextureFormat, TextureKind, TextureType } from './texture.js';
import type { UniformKind } from './uniform.js';

/**
 * Compute shaders: kernels dispatched over a grid of workgroups, reading and writing storage buffers and storage
 * textures. Available where `GraphicsContext.compute` is defined (WebGPU); algorithms keep a fallback without it.
 */

const getNextKernelCodeId = idFactory();

/** A kernel by name, like `ShaderCode` for programs: the source is registered per backend (WGSL in `gl/webgpu/`) */
export interface KernelCode {
  readonly id: number;
  readonly name: string;
  /** invocations per workgroup, declared in the kernel as `WorkgroupSizeX`, `WorkgroupSizeY`, `WorkgroupSizeZ` */
  readonly workgroupSize: readonly [number, number, number];
  /** sources of kernels generated at runtime (e.g. grid3d), otherwise backends look the name up in their registry */
  readonly sources?: { readonly wgsl?: string };
}

export function KernelCode(
  name: string,
  workgroupSize: [number, number?, number?],
  sources?: KernelCode['sources'],
): KernelCode {
  const [x, y = 1, z = 1] = workgroupSize;
  return { id: getNextKernelCodeId(), name, workgroupSize: [x, y, z], sources };
}

export type StorageBufferAccess = 'read' | 'read-write';

/**
 * A storage buffer binding, declared as a runtime-sized array of `element`: a WGSL type such as `u32`, `f32`,
 * `vec4<f32>`, `atomic<u32>`, or a struct declared by the kernel. In render schemas (WebGPU only, WebGL ignores them):
 * read-only, bound with the object, for the render items that reference it (`usedIf`)
 */
export type StorageBufferSpec = {
  type: 'storage-buffer';
  kind: 'storage-buffer';
  element: string;
  access: StorageBufferAccess;
  usedIf?: UsedIf;
};
export function StorageBufferSpec(
  element: string,
  access: StorageBufferAccess = 'read',
  usedIf?: UsedIf,
): StorageBufferSpec {
  return { type: 'storage-buffer', kind: 'storage-buffer', element, access, usedIf };
}

export type StorageTextureAccess = 'write' | 'read' | 'read-write';

/** A storage texture binding of a texture with `format` and `dataType` (formats with storage support, e.g. `rgba` `ubyte`) */
export type StorageTextureSpec = {
  type: 'storage-texture';
  kind: 'storage-texture';
  format: TextureFormat;
  dataType: TextureType;
  dimension: '2d' | '3d';
  access: StorageTextureAccess;
};
export function StorageTextureSpec(
  format: TextureFormat,
  dataType: TextureType,
  dimension: '2d' | '3d',
  access: StorageTextureAccess = 'write',
): StorageTextureSpec {
  return { type: 'storage-texture', kind: 'storage-texture', format, dataType, dimension, access };
}

/** Uniforms (one block), sampled textures, defines, storage buffers and storage textures of a kernel */
export type KernelSchema = {
  readonly [k: string]:
    | UniformSpec<UniformKind>
    | TextureSpec<TextureKind>
    | DefineSpec<DefineKind>
    | StorageBufferSpec
    | StorageTextureSpec;
};

export type KernelKindValue = KindValue & { 'storage-buffer': StorageBuffer; 'storage-texture': Texture };
export type KernelValues<S extends KernelSchema> = {
  readonly [k in keyof S]: ValueCell<KernelKindValue[S[k]['kind']]>;
};

/** GPU memory for kernels; a multiple of 4 bytes */
export interface StorageBuffer {
  readonly id: number;
  readonly byteLength: number;
  /** Writes `data` at `byteOffset`, before the commands of the current submission run (a multiple of 4 bytes) */
  write: (data: ArrayBufferView, byteOffset?: number) => void;
  /** Zeroes the range, in command order */
  clear: (byteOffset?: number, byteLength?: number) => void;
  /** Copies a range into `target`, in command order */
  copyTo: (target: StorageBuffer, srcByteOffset?: number, dstByteOffset?: number, byteLength?: number) => void;
  /** The contents after the commands recorded so far; tracked as pending work */
  readAsync: (byteOffset?: number, byteLength?: number) => Promise<ArrayBuffer>;
  getByteCount: () => number;
  destroy: () => void;
}

/**
 * Draw arguments in a storage buffer at `byteOffset` (a multiple of 4): four u32, the vertex count, instance count,
 * first vertex and first instance (zero). The producer writes the vertex count on the GPU, the render item drawing with
 * them the instance count. For non-indexed draws.
 */
export interface DrawIndirect {
  readonly buffer: StorageBuffer;
  readonly byteOffset: number;
  /**
   * The vertex count in the buffer once read back: drawn directly with it then, an indirect draw costs far more GPU
   * time than a direct one in Chrome (`TextureMesh.setDrawnVertexCount`)
   */
  readonly count?: number;
}

export interface ComputeKernel {
  /** Applies changed define values (a different variant of the kernel) */
  update: () => void;
  /** Dispatches `x` × `y` × `z` workgroups in the open compute pass */
  dispatch: (x: number, y?: number, z?: number) => void;
  /** Dispatches the workgroup counts of three `u32` in `buffer` at `byteOffset` */
  dispatchIndirect: (buffer: StorageBuffer, byteOffset?: number) => void;
  getByteCount: () => number;
  destroy: () => void;
}

/**
 * A kernel to prepare (`ComputeContext.prepareKernels`): the pipeline `createKernel` uses for the same code, schema and
 * define values
 */
export interface KernelRequest {
  readonly code: KernelCode;
  readonly schema: KernelSchema;
  /** the values of the schema's defines */
  readonly defines?: { readonly [k: string]: string | number | boolean };
}

export interface ComputeContext {
  /** Largest storage buffer a kernel can bind, in bytes */
  readonly maxStorageBufferSize: number;
  /** Bytes of the storage buffers created and not destroyed yet */
  getStorageBufferByteCount: () => number;
  /**
   * Whether the backend has the kernel; with `schema`, also whether its storage buffers fit in a stage (WebGPU's
   * `maxStorageBuffersPerShaderStage`: 8 under the core limits)
   */
  isKernelSupported: (code: KernelCode, schema?: KernelSchema) => boolean;
  createKernel: <S extends KernelSchema>(code: KernelCode, schema: S, values: KernelValues<S>) => ComputeKernel;
  /**
   * Compiles the pipelines of `requests` in the background (tracked as pending work) and resolves when they are ready;
   * kernels created afterwards use them, a kernel created before compiles its pipeline synchronously, blocking the
   * queue until compiled
   */
  prepareKernels: (requests: Iterable<KernelRequest>) => Promise<void>;
  /**
   * Whether the pipelines of `requests` are ready (prepared or created before): awaiting `prepareKernels` yields even
   * then, which lets other work run in between (e.g. work sharing a texture)
   */
  hasKernels: (requests: Iterable<KernelRequest>) => boolean;
  createBuffer: (byteLength: number, label?: string) => StorageBuffer;
  /**
   * As `createBuffer`, resolving once the device created it, or to `undefined` when it ran out of memory (WebGPU: an
   * error scope): a buffer that failed would invalidate the commands using it, those recorded with them too
   */
  tryCreateBuffer: (byteLength: number, label?: string) => Promise<StorageBuffer | undefined>;
  /** Whether `buffer` was created by this context (not by one it replaced, e.g. before a restore or backend switch) */
  ownsBuffer: (buffer: StorageBuffer) => boolean;
  /**
   * Ends an open render or compute pass and begins a compute pass with `label` (timed like render passes); dispatches
   * go into it until `endComputePass`, a render pass, a copy or a submit ends it
   */
  beginComputePass: (options: { label: string }) => void;
  endComputePass: () => void;
}
