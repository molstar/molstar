/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
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
 * Color smoothing with compute kernels: the colors of sampled vertices accumulated into the cells of a grid, weighted
 * by distance (fixed-point atomic sums), then normalized into the 2d atlas of the grid (`color-smoothing-normalize`).
 * The vertices of texture meshes in the vertex arena, read through the page table (`color-smoothing-accumulate`, see
 * `computeChunkedIsosurface`), as the point sprite accumulation on WebGL; or of a mesh in storage buffers
 * (`color-smoothing-accumulate-mesh`, `GpuMeshVertices`), as the CPU's splat (`calcMeshColorSmoothing`). The sums span
 * a region of the grid (the reached instances' cells; all of it by default), in slabs of slices where they exceed the
 * storage buffer limit; the atlas outside the region zero.
 */

export const ColorAccumulateKernel = KernelCode('color-smoothing-accumulate', [64]);
export const ColorAccumulateMeshKernel = KernelCode('color-smoothing-accumulate-mesh', [64]);
export const ColorNormalizeKernel = KernelCode('color-smoothing-normalize', [8, 8]);

const AccumulateUniforms = {
  uSampleCount: UniformSpec('i'),
  uInstanceCount: UniformSpec('i'),
  uStride: UniformSpec('i'),
  uGroupCount: UniformSpec('i'),
  uGridDim: UniformSpec('iv3'),
  /** first slice and slice count of the slab, within the region */
  uSlab: UniformSpec('iv2'),
  /** the region of the grid the sums span: its first cell and the cell past its last */
  uRegionMin: UniformSpec('iv3'),
  uRegionMax: UniformSpec('iv3'),
  /** the grid's minimum */
  uBboxMin: UniformSpec('v3'),
  uResolution: UniformSpec('f'),
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  dColorType: DefineSpec('string', ['group', 'groupInstance']),
};

export const ColorAccumulateSchema = {
  ...AccumulateUniforms,
  /** an upper bound of the vertex count, the drawn count is in the page table */
  uVertexCount: UniformSpec('i'),
  /** the grid's size, for the point sprite's cells */
  uBboxSize: UniformSpec('v3'),
  isoTable: StorageBufferSpec('u32', 'read'),
  isoVertices: StorageBufferSpec('u32', 'read'),
  transforms: StorageBufferSpec('mat4x4<f32>', 'read'),
  instances: StorageBufferSpec('u32', 'read'),
  sums: StorageBufferSpec('atomic<u32>', 'read-write'),
};

export const ColorAccumulateMeshSchema = {
  ...AccumulateUniforms,
  uVertexCount: UniformSpec('i'),
  /** the mesh's first vertex in `positions` and `vertexGroups` */
  uVertexOffset: UniformSpec('i'),
  /** xyz per vertex */
  positions: StorageBufferSpec('f32', 'read'),
  /** the group per vertex */
  vertexGroups: StorageBufferSpec('u32', 'read'),
  transforms: StorageBufferSpec('mat4x4<f32>', 'read'),
  instances: StorageBufferSpec('u32', 'read'),
  sums: StorageBufferSpec('atomic<u32>', 'read-write'),
};

export const ColorNormalizeSchema = {
  uGridDim: UniformSpec('iv3'),
  uSlab: UniformSpec('iv2'),
  uRegionMin: UniformSpec('iv3'),
  uRegionMax: UniformSpec('iv3'),
  uTexCols: UniformSpec('i'),
  sums: StorageBufferSpec('u32', 'read'),
  outColor: StorageTextureSpec('rgba', 'ubyte', '2d', 'write'),
};

/** The kernels of `computeColorSmoothing` for `colorType`, see `ComputeContext.prepareKernels` */
export function getColorSmoothingKernels(colorType: 'group' | 'groupInstance'): KernelRequest[] {
  return [
    { code: ColorAccumulateKernel, schema: ColorAccumulateSchema, defines: { dColorType: colorType } },
    { code: ColorNormalizeKernel, schema: ColorNormalizeSchema },
  ];
}

/** The kernels of `computeMeshColorSmoothing` for `colorType`, see `ComputeContext.prepareKernels` */
export function getMeshColorSmoothingKernels(colorType: 'group' | 'groupInstance'): KernelRequest[] {
  return [
    { code: ColorAccumulateMeshKernel, schema: ColorAccumulateMeshSchema, defines: { dColorType: colorType } },
    { code: ColorNormalizeKernel, schema: ColorNormalizeSchema },
  ];
}

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

export function isComputeMeshColorSmoothingSupported(
  graphics: GraphicsContext,
): graphics is GraphicsContext & { compute: ComputeContext } {
  const { compute } = graphics;
  return (
    !!compute &&
    compute.isKernelSupported(ColorAccumulateMeshKernel, ColorAccumulateMeshSchema) &&
    compute.isKernelSupported(ColorNormalizeKernel)
  );
}

interface SampledInstances {
  /** every `stride`-th vertex is sampled */
  stride: number;
  /** the samples, `round(vertexCount / stride)` when not given (the CPU's: `ceil`) */
  sampleCount?: number;
  /** the instances splatted, their transforms and rows (`instanceBuffer`, of the colors) */
  instanceCount: number;
  groupCount: number;
  transformBuffer: Float32Array;
  instanceBuffer: Uint32Array;
  colorData: Texture;
  colorType: 'group' | 'groupInstance';
}

export interface ColorSmoothingSamples extends SampledInstances {
  /** an upper bound, the drawn count is in the page table */
  vertexCount: number;
  /** the page table and positions of the surface in the vertex arena */
  chunks: { table: StorageBuffer; vertices: StorageBuffer };
}

export interface MeshColorSmoothingSamples extends SampledInstances {
  vertexCount: number;
  /** xyz per vertex and the group per vertex (`GpuMeshVertices`), the mesh's from `vertexOffset` */
  positions: StorageBuffer;
  groups: StorageBuffer;
  vertexOffset: number;
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
  /**
   * the cells the sums span (min x, y, z, exclusive max x, y, z; `getSplatInstances`' region): the atlas zero outside
   * it; all of the grid when not given
   */
  region?: ArrayLike<number>;
}

/** Writes the smoothed colors of `samples` (texture meshes) on `grid` into the rgba `texture` (defined to the atlas size) */
export function computeColorSmoothing(
  compute: ComputeContext,
  samples: ColorSmoothingSamples,
  grid: ColorSmoothingGrid,
  texture: Texture,
) {
  smooth(compute, samples, grid, texture, (uniforms, buffers) =>
    compute.createKernel(ColorAccumulateKernel, ColorAccumulateSchema, {
      ...uniforms,
      ...buffers,
      uVertexCount: ValueCell.create(samples.vertexCount),
      uBboxSize: ValueCell.create(grid.bboxSize),
      isoTable: ValueCell.create(samples.chunks.table),
      isoVertices: ValueCell.create(samples.chunks.vertices),
    }),
  );
}

/** Writes the smoothed colors of `samples` (a mesh's) on `grid` into the rgba `texture` (defined to the atlas size) */
export function computeMeshColorSmoothing(
  compute: ComputeContext,
  samples: MeshColorSmoothingSamples,
  grid: ColorSmoothingGrid,
  texture: Texture,
) {
  smooth(compute, samples, grid, texture, (uniforms, buffers) =>
    compute.createKernel(ColorAccumulateMeshKernel, ColorAccumulateMeshSchema, {
      ...uniforms,
      ...buffers,
      uVertexCount: ValueCell.create(samples.vertexCount),
      uVertexOffset: ValueCell.create(samples.vertexOffset),
      positions: ValueCell.create(samples.positions),
      vertexGroups: ValueCell.create(samples.groups),
    }),
  );
}

type AccumulateValues = {
  [K in keyof typeof AccumulateUniforms]: ValueCell<any>;
};

function smooth(
  compute: ComputeContext,
  samples: SampledInstances & { vertexCount: number },
  grid: ColorSmoothingGrid,
  texture: Texture,
  createAccumulate: (
    uniforms: AccumulateValues,
    buffers: {
      transforms: ValueCell<StorageBuffer>;
      instances: ValueCell<StorageBuffer>;
      sums: ValueCell<StorageBuffer>;
    },
  ) => ReturnType<ComputeContext['createKernel']>,
) {
  const { gridDim, texCols, width, height } = grid;
  texture.define(width, height);

  const region = grid.region ?? [0, 0, 0, gridDim[0], gridDim[1], gridDim[2]];
  const regionMin = Vec3.create(region[0], region[1], region[2]);
  const regionMax = Vec3.create(region[3], region[4], region[5]);
  const depth = regionMax[2] - regionMin[2];
  const sliceSize = (regionMax[0] - regionMin[0]) * (regionMax[1] - regionMin[1]) * SumsPerCell * 4;
  const slabDepth = Math.max(0, Math.min(depth, Math.floor(compute.maxStorageBufferSize / sliceSize)));
  if (depth > 0 && slabDepth < 1) throw new Error(`a color smoothing grid slice exceeds the storage buffer limit`);
  const sums = compute.createBuffer(Math.max(4, slabDepth * sliceSize), 'color smoothing sums');
  const transforms = compute.createBuffer(
    Math.max(64, samples.transformBuffer.byteLength),
    'color smoothing transforms',
  );
  transforms.write(samples.transformBuffer);
  const instances = compute.createBuffer(Math.max(4, samples.instanceBuffer.byteLength), 'color smoothing instances');
  instances.write(samples.instanceBuffer);

  const sampleCount = samples.sampleCount ?? Math.round(samples.vertexCount / samples.stride);
  const gridDimValue = ValueCell.create(gridDim);
  const slab = ValueCell.create(Vec2.create(regionMin[2], slabDepth));
  const regionMinValue = ValueCell.create(regionMin);
  const regionMaxValue = ValueCell.create(regionMax);
  const accumulateKernel = createAccumulate(
    {
      uSampleCount: ValueCell.create(sampleCount),
      uInstanceCount: ValueCell.create(samples.instanceCount),
      uStride: ValueCell.create(samples.stride),
      uGroupCount: ValueCell.create(samples.groupCount),
      uGridDim: gridDimValue,
      uSlab: slab,
      uRegionMin: regionMinValue,
      uRegionMax: regionMaxValue,
      uBboxMin: ValueCell.create(grid.bboxMin),
      uResolution: ValueCell.create(grid.resolution),
      tColor: ValueCell.create(samples.colorData),
      dColorType: ValueCell.create(samples.colorType),
    },
    { transforms: ValueCell.create(transforms), instances: ValueCell.create(instances), sums: ValueCell.create(sums) },
  );
  const normalizeKernel = compute.createKernel(ColorNormalizeKernel, ColorNormalizeSchema, {
    uGridDim: gridDimValue,
    uSlab: slab,
    uRegionMin: regionMinValue,
    uRegionMax: regionMaxValue,
    uTexCols: ValueCell.create(texCols),
    sums: ValueCell.create(sums),
    outColor: ValueCell.create(texture),
  });

  const { x, y } = getDispatchSize(Math.ceil((sampleCount * samples.instanceCount) / 64), MaxWorkgroups);
  // an empty region: one pass, the atlas zeroed (the first slab's, the region's start)
  const end = Math.max(regionMax[2], regionMin[2] + 1);
  for (let z = regionMin[2]; z < end; z += Math.max(1, slabDepth)) {
    // the clear ends the open compute pass
    if (z > regionMin[2]) sums.clear();
    ValueCell.update(slab, Vec2.set(slab.ref.value, z, slabDepth));
    compute.beginComputePass({ label: 'color smoothing' });
    if (slabDepth > 0) accumulateKernel.dispatch(x, y);
    normalizeKernel.dispatch(Math.ceil(width / 8), Math.ceil(height / 8));
    compute.endComputePass();
  }

  accumulateKernel.destroy();
  normalizeKernel.destroy();
  sums.destroy();
  transforms.destroy();
  instances.destroy();
}
