/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { TriTable } from '@molstar/graphics/geo/util/marching-cubes/tables';
import {
  type ComputeContext,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
} from '../../compute.js';
import type { GraphicsContext } from '../../context.js';
import { DefineSpec, TextureSpec, UniformSpec } from '../../renderable/schema.js';
import type { Texture } from '../../texture.js';

/**
 * Marching cubes with compute kernels, into chunks of a vertex arena (`computeChunkedIsosurface`): the inputs, the
 * kernels and the tables they share. Results follow the WebGL isosurface: positions and normals in the space of the
 * transform; indexed, a group per vertex (of the nearer edge end, as the CPU marching cubes), a triangle's of its first.
 */

export const WorkgroupSize: [number, number, number] = [4, 4, 4];

export const McChunkCountKernel = KernelCode('mc-chunk-count', WorkgroupSize);
export const McChunkEmitKernel = KernelCode('mc-chunk-emit', WorkgroupSize);
export const McChunkAllocateKernel = KernelCode('mc-chunk-allocate', [1]);
export const McChunkFreeKernel = KernelCode('mc-chunk-free', [1]);
export const McChunkPushKernel = KernelCode('mc-chunk-push', [1]);

const InputUniforms = {
  uGridDim: UniformSpec('iv3'),
  uIsoValue: UniformSpec('f'),
};

/** Volume data in a storage buffer of f32, cell `(x, y, z)` at `dot((x, y, z), uStrides)`; groups are the data index */
export const McBufferInput = {
  ...InputUniforms,
  uStrides: UniformSpec('iv3'),
  /** points of the data per axis: grid positions past them wrap (periodic volumes, the grid one point larger) */
  uDataDim: UniformSpec('iv3'),
  volume: StorageBufferSpec('f32', 'read'),
  dInput: DefineSpec('string', ['buffer', 'texture']),
};

/** Volume data in a 3d texture, the value in `dChannel`; with `dPackedGroup`, the group packed in rgb */
export const McTextureInput = {
  ...InputUniforms,
  tVolume: TextureSpec('texture', 'rgba', 'float', 'nearest', undefined, { dimension: '3d' }),
  dInput: DefineSpec('string', ['buffer', 'texture']),
  dChannel: DefineSpec('string', ['red', 'alpha']),
  dPackedGroup: DefineSpec('boolean'),
};

export const McChunkCountOutput = {
  triCounts: StorageBufferSpec('u32', 'read'),
  table: StorageBufferSpec('atomic<u32>', 'read-write'),
  bases: StorageBufferSpec('u32', 'read-write'),
  /** the first vertex of each grid point */
  pointBases: StorageBufferSpec('u32', 'read-write'),
};

export const McChunkAllocateSchema = {
  /** the box of the surface's quantized positions: its min and its extent per step */
  uBoxMin: UniformSpec('v3'),
  uBoxStep: UniformSpec('v3'),
  /** the vertex pages of the table, the triangle pages after them */
  uVertexCapacity: UniformSpec('i'),
  table: StorageBufferSpec('u32', 'read-write'),
  free: StorageBufferSpec('u32', 'read-write'),
  drawArgs: StorageBufferSpec('u32', 'read-write'),
};

export const McChunkEmitOutput = {
  uInvert: UniformSpec('b'),
  uTransform: UniformSpec('m4'),
  uTransformAdjoint: UniformSpec('m3'),
  /** the box of the quantized positions: its min and the steps per unit */
  uBoxMin: UniformSpec('v3'),
  uBoxScale: UniformSpec('v3'),
  triCounts: StorageBufferSpec('u32', 'read'),
  triTable: StorageBufferSpec('u32', 'read'),
  bases: StorageBufferSpec('u32', 'read'),
  pointBases: StorageBufferSpec('u32', 'read'),
  table: StorageBufferSpec('u32', 'read'),
  vertices: StorageBufferSpec('u32', 'read-write'),
};

export const McChunkFreeSchema = {
  table: StorageBufferSpec('u32', 'read-write'),
  free: StorageBufferSpec('u32', 'read-write'),
};

export const McChunkPushSchema = {
  uFirst: UniformSpec('i'),
  uCount: UniformSpec('i'),
  free: StorageBufferSpec('u32', 'read-write'),
};

export type IsosurfaceInput =
  | {
      readonly kind: 'buffer';
      /** f32 values */
      readonly volume: StorageBuffer;
      /** data offset per axis step */
      readonly strides: Vec3;
      /** the data's points per axis where the grid wraps around them (one point more per axis); the grid's by default */
      readonly dataDim?: Vec3;
    }
  | {
      readonly kind: 'texture';
      /** a 3d texture */
      readonly texture: Texture;
      readonly channel: 'red' | 'alpha';
      readonly packedGroup: boolean;
    };

export function isComputeIsosurfaceSupported(
  graphics: GraphicsContext,
): graphics is GraphicsContext & { compute: ComputeContext } {
  const { compute } = graphics;
  if (!compute) return false;
  return [McChunkCountKernel, McChunkAllocateKernel, McChunkEmitKernel, McChunkFreeKernel, McChunkPushKernel].every(
    (k) => compute.isKernelSupported(k),
  );
}

type Shared = { triCounts: StorageBuffer; triTable: StorageBuffer };
const shared = new WeakMap<ComputeContext, Shared>();

/** The marching cubes tables per context */
export function getShared(compute: ComputeContext): Shared {
  let s = shared.get(compute);
  if (!s) {
    const triCounts = new Uint32Array(256);
    const triTable = new Uint32Array(256 * 16).fill(255);
    for (let i = 0; i < 256; ++i) {
      triCounts[i] = TriTable[i].length / 3;
      for (let j = 0; j < TriTable[i].length; ++j) triTable[i * 16 + j] = TriTable[i][j];
    }
    s = {
      triCounts: compute.createBuffer(triCounts.byteLength, 'mc tri counts'),
      triTable: compute.createBuffer(triTable.byteLength, 'mc tri table'),
    };
    s.triCounts.write(triCounts);
    s.triTable.write(triTable);
    shared.set(compute, s);
  }
  return s;
}

/** What selects the kernels of an input: its kind and how its values are read */
export type IsosurfaceInputKind =
  | { readonly kind: 'buffer' }
  | { readonly kind: 'texture'; readonly channel: 'red' | 'alpha'; readonly packedGroup: boolean };

export function getInputKernelSchema(input: IsosurfaceInputKind): {
  schema: typeof McBufferInput | typeof McTextureInput;
  defines: { [k: string]: string | boolean };
} {
  if (input.kind === 'buffer') {
    return { schema: McBufferInput, defines: { dInput: 'buffer' } };
  }
  return {
    schema: McTextureInput,
    defines: { dInput: 'texture', dChannel: input.channel, dPackedGroup: input.packedGroup },
  };
}

export function getInputValues(input: IsosurfaceInput, gridDim: Vec3, isoValue: number) {
  const { schema, defines } = getInputKernelSchema(input);
  const common = { uGridDim: ValueCell.create(Vec3.clone(gridDim)), uIsoValue: ValueCell.create(isoValue) };
  const defineValues = Object.fromEntries(Object.entries(defines).map(([k, v]) => [k, ValueCell.create(v)]));
  if (input.kind === 'buffer') {
    return {
      schema,
      values: {
        ...common,
        uStrides: ValueCell.create(Vec3.clone(input.strides)),
        uDataDim: ValueCell.create(Vec3.clone(input.dataDim ?? gridDim)),
        volume: ValueCell.create(input.volume),
        ...defineValues,
      },
    };
  }
  return { schema, values: { ...common, tVolume: ValueCell.create(input.texture), ...defineValues } };
}

/** The kernels of `computeChunkedIsosurface` for `input`, see `ComputeContext.prepareKernels` */
export function getIsosurfaceKernels(input: IsosurfaceInputKind): KernelRequest[] {
  const { schema, defines } = getInputKernelSchema(input);
  return [
    { code: McChunkCountKernel, schema: { ...schema, ...McChunkCountOutput }, defines },
    { code: McChunkAllocateKernel, schema: McChunkAllocateSchema },
    { code: McChunkEmitKernel, schema: { ...schema, ...McChunkEmitOutput }, defines },
    { code: McChunkFreeKernel, schema: McChunkFreeSchema },
    { code: McChunkPushKernel, schema: McChunkPushSchema },
  ];
}
