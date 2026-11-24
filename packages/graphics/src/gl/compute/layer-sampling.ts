/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import {
  type ComputeContext,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
} from '../compute.js';
import type { GraphicsContext } from '../context.js';
import { DefineSpec, TextureSpec, UniformSpec } from '../renderable/schema.js';
import type { Texture } from '../texture.js';
import { getDispatchSize } from './prefix-sum.js';

/**
 * A smoothed layer sampled per vertex and instance on the GPU (`vertexInstance` data of a mesh): its grid (the 2d atlas
 * of `calcMeshColorSmoothing`, uploaded linear) sampled by the texture's filter at each vertex of each instance, as the
 * mesh shaders sample a layer's grid, into a storage buffer of `rgba8` or `r8` elements (`DataLayers`)
 */

export const LayerSampleKernel = KernelCode('layer-sample-vertices', [64]);

export const LayerSampleSchema = {
  uVertexCount: UniformSpec('i'),
  uInstanceCount: UniformSpec('i'),
  /** the mesh's first vertex in `positions` */
  uVertexOffset: UniformSpec('i'),
  /** the mesh's first instance in `transforms` */
  uInstanceOffset: UniformSpec('i'),
  /** the first element written in `outData` */
  uOutOffset: UniformSpec('i'),
  uGridDim: UniformSpec('v3'),
  uGridTexDim: UniformSpec('v2'),
  /** xyz: the grid's minimum, w: its scale (1 / resolution) */
  uGridTransform: UniformSpec('v4'),
  tGrid: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  dItemSize: DefineSpec('string', ['1', '4']),
  /** xyz per vertex */
  positions: StorageBufferSpec('f32', 'read'),
  transforms: StorageBufferSpec('mat4x4<f32>', 'read'),
  /** per instance: non-zero where the layer reaches it */
  reached: StorageBufferSpec('u32', 'read'),
  outData: StorageBufferSpec('u32', 'read-write'),
};

/** Per-dimension dispatch limit of the WebGPU core */
const MaxWorkgroups = 65535;

/** The kernels of `computeLayerSampling` for `itemSize`, see `ComputeContext.prepareKernels` */
export function getLayerSamplingKernels(itemSize: 1 | 4): KernelRequest[] {
  return [{ code: LayerSampleKernel, schema: LayerSampleSchema, defines: { dItemSize: `${itemSize}` } }];
}

/** The kernels of both item sizes */
export const LayerSamplingKernels: readonly KernelRequest[] = [
  ...getLayerSamplingKernels(1),
  ...getLayerSamplingKernels(4),
];

/**
 * Compiles the kernels of `computeLayerSampling` in the background where they run and are not ready (else the first
 * smoothed layer compiles them synchronously, blocking the queue); `undefined` where there is nothing to wait for
 */
export function prepareLayerSampling(graphics: GraphicsContext): Promise<void> | undefined {
  if (!isComputeLayerSamplingSupported(graphics) || graphics.compute.hasKernels(LayerSamplingKernels)) return;
  return graphics.compute.prepareKernels(LayerSamplingKernels);
}

export function isComputeLayerSamplingSupported(
  graphics: GraphicsContext,
): graphics is GraphicsContext & { compute: ComputeContext } {
  const { compute } = graphics;
  return !!compute && compute.isKernelSupported(LayerSampleKernel, LayerSampleSchema);
}

export interface LayerSamplingInput {
  vertexCount: number;
  instanceCount: number;
  /** xyz per vertex, the mesh's from `vertexOffset` */
  positions: StorageBuffer;
  vertexOffset: number;
  /** a matrix per instance, the mesh's from `instanceOffset` */
  transforms: StorageBuffer;
  instanceOffset: number;
  /** a `u32` per instance of the mesh, non-zero where the layer reaches it */
  reached: StorageBuffer;
  /** the grid (`calcMeshColorSmoothing`'s 2d atlas of slices), uploaded with linear filtering */
  grid: Texture;
  gridDim: Vec3;
  gridTexDim: Vec2;
  gridTransform: Vec4;
  /** 4: `rgba8` (overpaint, substance), 1: `r8` (transparency, emissive) */
  itemSize: 1 | 4;
}

/** Writes the layer of `input` into `out` from element `outOffset`, in a compute pass of its own */
export function computeLayerSampling(
  compute: ComputeContext,
  input: LayerSamplingInput,
  out: StorageBuffer,
  outOffset: number,
) {
  const kernel = compute.createKernel(LayerSampleKernel, LayerSampleSchema, {
    uVertexCount: ValueCell.create(input.vertexCount),
    uInstanceCount: ValueCell.create(input.instanceCount),
    uVertexOffset: ValueCell.create(input.vertexOffset),
    uInstanceOffset: ValueCell.create(input.instanceOffset),
    uOutOffset: ValueCell.create(outOffset),
    uGridDim: ValueCell.create(input.gridDim),
    uGridTexDim: ValueCell.create(input.gridTexDim),
    uGridTransform: ValueCell.create(input.gridTransform),
    tGrid: ValueCell.create(input.grid),
    dItemSize: ValueCell.create(`${input.itemSize}`),
    positions: ValueCell.create(input.positions),
    transforms: ValueCell.create(input.transforms),
    reached: ValueCell.create(input.reached),
    outData: ValueCell.create(out),
  });
  const count = input.vertexCount * input.instanceCount;
  // a thread an element (rgba8) or a word of four (r8)
  const threads = input.itemSize === 4 ? count : Math.ceil(((outOffset % 4) + count) / 4);
  const { x, y } = getDispatchSize(Math.ceil(threads / 64), MaxWorkgroups);
  compute.beginComputePass({ label: 'layer sampling' });
  kernel.dispatch(x, y);
  compute.endComputePass();
  kernel.destroy();
}
