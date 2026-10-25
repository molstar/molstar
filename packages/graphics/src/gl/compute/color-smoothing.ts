/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec2, type Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import {
  type ComputeContext,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
  StorageTextureSpec,
} from '../compute.js';
import type { GraphicsContext } from '../context.js';
import { DefineSpec, TextureSpec, UniformSpec } from '../renderable/schema.js';
import type { Texture } from '../texture.js';
import { getDispatchSize } from './prefix-sum.js';

/**
 * Color smoothing of texture meshes in the vertex arena (see `computeChunkedIsosurface`) with compute kernels: the
 * colors of sampled vertices, read through the page table, are accumulated into the cells
 * of a grid, weighted by distance (`color-smoothing-accumulate`, fixed-point atomic sums), then normalized into the 2d
 * atlas of the grid (`color-smoothing-normalize`), as the point sprite accumulation on WebGL. Grids whose sums exceed
 * the storage buffer limit are processed in slabs of slices.
 */

export const ColorAccumulateKernel = KernelCode('color-smoothing-accumulate', [64]);
export const ColorNormalizeKernel = KernelCode('color-smoothing-normalize', [8, 8]);

export const ColorAccumulateSchema = {
  /** an upper bound of the vertex count, the drawn count is in the page table */
  uVertexCount: UniformSpec('i'),
  uSampleCount: UniformSpec('i'),
  uInstanceCount: UniformSpec('i'),
  uStride: UniformSpec('i'),
  uGroupCount: UniformSpec('i'),
  uGridDim: UniformSpec('iv3'),
  /** first slice and slice count of the slab */
  uSlab: UniformSpec('iv2'),
  uBboxMin: UniformSpec('v3'),
  uBboxSize: UniformSpec('v3'),
  uResolution: UniformSpec('f'),
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  dColorType: DefineSpec('string', ['group', 'groupInstance']),
  isoTable: StorageBufferSpec('u32', 'read'),
  isoVertices: StorageBufferSpec('u32', 'read'),
  transforms: StorageBufferSpec('mat4x4<f32>', 'read'),
  instances: StorageBufferSpec('f32', 'read'),
  sums: StorageBufferSpec('atomic<u32>', 'read-write'),
};

/** The kernels of `computeColorSmoothing` for `colorType`, see `ComputeContext.prepareKernels` */
export function getColorSmoothingKernels(colorType: 'group' | 'groupInstance'): KernelRequest[] {
  return [
    { code: ColorAccumulateKernel, schema: ColorAccumulateSchema, defines: { dColorType: colorType } },
    { code: ColorNormalizeKernel, schema: ColorNormalizeSchema },
  ];
}

export const ColorNormalizeSchema = {
  uGridDim: UniformSpec('iv3'),
  uSlab: UniformSpec('iv2'),
  uTexCols: UniformSpec('i'),
  sums: StorageBufferSpec('u32', 'read'),
  outColor: StorageTextureSpec('rgba', 'ubyte', '2d', 'write'),
};

/** Per grid cell: weighted rgba and the weight */
const SumsPerCell = 5;

/** Per-dimension dispatch limit of the WebGPU core */
const MaxWorkgroups = 65535;

export function isComputeColorSmoothingSupported(
  graphics: GraphicsContext,
): graphics is GraphicsContext & { compute: ComputeContext } {
  const { compute } = graphics;
  return (
    !!compute && compute.isKernelSupported(ColorAccumulateKernel) && compute.isKernelSupported(ColorNormalizeKernel)
  );
}

export interface ColorSmoothingSamples {
  /** an upper bound, the drawn count is in the page table */
  vertexCount: number;
  /** every `stride`-th vertex is sampled */
  stride: number;
  instanceCount: number;
  groupCount: number;
  transformBuffer: Float32Array;
  instanceBuffer: Float32Array;
  /** the page table and positions of the surface in the vertex arena */
  chunks: { table: StorageBuffer; vertices: StorageBuffer };
  colorData: Texture;
  colorType: 'group' | 'groupInstance';
}

export interface ColorSmoothingGrid {
  gridDim: Vec3;
  /** slices per row of the atlas */
  texCols: number;
  /** size of the atlas, `texCols` × `gridDim[0]` wide */
  width: number;
  height: number;
  bboxMin: Vec3;
  bboxSize: Vec3;
  resolution: number;
}

/** Writes the smoothed colors of `samples` on `grid` into the rgba `texture` (defined to the atlas size) */
export function computeColorSmoothing(
  compute: ComputeContext,
  samples: ColorSmoothingSamples,
  grid: ColorSmoothingGrid,
  texture: Texture,
) {
  const { gridDim, texCols, width, height } = grid;
  texture.define(width, height);

  const sliceSize = gridDim[0] * gridDim[1] * SumsPerCell * 4;
  const slabDepth = Math.min(gridDim[2], Math.floor(compute.maxStorageBufferSize / sliceSize));
  if (slabDepth < 1) throw new Error(`a color smoothing grid slice exceeds the storage buffer limit`);
  const sums = compute.createBuffer(slabDepth * sliceSize, 'color smoothing sums');
  const transforms = compute.createBuffer(samples.transformBuffer.byteLength, 'color smoothing transforms');
  transforms.write(samples.transformBuffer);
  const instances = compute.createBuffer(samples.instanceBuffer.byteLength, 'color smoothing instances');
  instances.write(samples.instanceBuffer);

  const sampleCount = Math.round(samples.vertexCount / samples.stride);
  const gridDimValue = ValueCell.create(gridDim);
  const slab = ValueCell.create(Vec2.create(0, slabDepth));
  const accumulateKernel = compute.createKernel(ColorAccumulateKernel, ColorAccumulateSchema, {
    isoTable: ValueCell.create(samples.chunks.table),
    isoVertices: ValueCell.create(samples.chunks.vertices),
    uVertexCount: ValueCell.create(samples.vertexCount),
    uSampleCount: ValueCell.create(sampleCount),
    uInstanceCount: ValueCell.create(samples.instanceCount),
    uStride: ValueCell.create(samples.stride),
    uGroupCount: ValueCell.create(samples.groupCount),
    uGridDim: gridDimValue,
    uSlab: slab,
    uBboxMin: ValueCell.create(grid.bboxMin),
    uBboxSize: ValueCell.create(grid.bboxSize),
    uResolution: ValueCell.create(grid.resolution),
    tColor: ValueCell.create(samples.colorData),
    dColorType: ValueCell.create(samples.colorType),
    transforms: ValueCell.create(transforms),
    instances: ValueCell.create(instances),
    sums: ValueCell.create(sums),
  });
  const normalizeKernel = compute.createKernel(ColorNormalizeKernel, ColorNormalizeSchema, {
    uGridDim: gridDimValue,
    uSlab: slab,
    uTexCols: ValueCell.create(texCols),
    sums: ValueCell.create(sums),
    outColor: ValueCell.create(texture),
  });

  const { x, y } = getDispatchSize(Math.ceil((sampleCount * samples.instanceCount) / 64), MaxWorkgroups);
  for (let z = 0; z < gridDim[2]; z += slabDepth) {
    // the clear ends the open compute pass
    if (z > 0) sums.clear();
    ValueCell.update(slab, Vec2.set(slab.ref.value, z, slabDepth));
    compute.beginComputePass({ label: 'color smoothing' });
    accumulateKernel.dispatch(x, y);
    normalizeKernel.dispatch(Math.ceil(width / 8), Math.ceil(height / 8));
    compute.endComputePass();
  }

  accumulateKernel.destroy();
  normalizeKernel.destroy();
  sums.destroy();
  transforms.destroy();
  instances.destroy();
}
