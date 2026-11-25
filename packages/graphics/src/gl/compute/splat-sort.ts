/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
} from '../compute.js';
import { TextureSpec, UniformSpec } from '../renderable/schema.js';
import { createNullTexture, type Texture } from '../texture.js';
import { createPrefixSum, getDispatchSize, type PrefixSum, PrefixSumKernels } from './prefix-sum.js';

const WorkgroupSize = 256;
/** Keys of 16 bits: a bucket each */
export const SplatKeyCount = 65536;

export const SplatKeysKernel = KernelCode('splat-keys', [WorkgroupSize]);
export const SplatScatterKernel = KernelCode('splat-scatter', [WorkgroupSize]);
export const SplatCentersKernel = KernelCode('splat-centers', [WorkgroupSize]);

export const SplatKeysSchema = {
  uCount: UniformSpec('i'),
  uGroupsX: UniformSpec('i'),
  /** the view z row of the instance (model-view, object and instance transform): `dot(row.xyz, center)` */
  uRow: UniformSpec('v4'),
  uDepthMin: UniformSpec('f'),
  /** key steps per unit of depth */
  uDepthScale: UniformSpec('f'),
  /** center xyz and group bits per splat (`tPositionGroup`'s layout) */
  positions: StorageBufferSpec('vec4<f32>', 'read'),
  keys: StorageBufferSpec('u32', 'read-write'),
  histogram: StorageBufferSpec('atomic<u32>', 'read-write'),
};

export const SplatScatterSchema = {
  uCount: UniformSpec('i'),
  uGroupsX: UniformSpec('i'),
  /** the instance's block in the order */
  uOffset: UniformSpec('i'),
  keys: StorageBufferSpec('u32', 'read'),
  offsets: StorageBufferSpec('atomic<u32>', 'read-write'),
  order: StorageBufferSpec('u32', 'read-write'),
};

/** SOG's means read in place (the `sog` encoding), decoded into the centers the keys read */
export const SplatCentersSchema = {
  uCount: UniformSpec('i'),
  uGroupsX: UniformSpec('i'),
  /** the means' ranges, in the log domain */
  uMeansMin: UniformSpec('v3'),
  uMeansMax: UniformSpec('v3'),
  tMeansLow: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tMeansHigh: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  positions: StorageBufferSpec('vec4<f32>', 'read-write'),
};

/** The kernels of a splat sort, see `ComputeContext.prepareKernels` */
export const SplatSortKernels: readonly KernelRequest[] = [
  { code: SplatKeysKernel, schema: SplatKeysSchema },
  { code: SplatScatterKernel, schema: SplatScatterSchema },
  { code: SplatCentersKernel, schema: SplatCentersSchema },
  ...PrefixSumKernels,
];

/** An instance's sort: the order's block from `offset` (`count` splats), back to front along the view row */
export interface SplatSortJob {
  readonly offset: number;
  readonly count: number;
  /** the view z row: `a x + b y + c z` of a center */
  readonly a: number;
  readonly b: number;
  readonly c: number;
  /** the depths of the instance's splats, from its bounds */
  readonly depthMin: number;
  readonly depthMax: number;
}

/**
 * Sorts gaussian splats back to front on the GPU (a counting sort of 16-bit depth keys: the keys and their histogram,
 * its prefix sum, the scatter; no readback), each instance into its block of an order buffer the vertex shader reads
 * (`sSplatOrder`)
 */
export interface SplatSorter {
  /** Sorts the instances of `jobs` into `order` from the centers in `positions` (`vec4` a splat) */
  sort: (positions: StorageBuffer, order: StorageBuffer, jobs: readonly SplatSortJob[]) => void;
  /** Decodes the centers of `count` SOG splats (the `sog` encoding) from its means into `positions` */
  decodeSogCenters: (
    meansLow: Texture,
    meansHigh: Texture,
    meansMin: Vec3,
    meansMax: Vec3,
    count: number,
    positions: StorageBuffer,
  ) => void;
  destroy: () => void;
}

export function isSplatSortSupported(compute: ComputeContext) {
  return SplatSortKernels.every((k) => compute.isKernelSupported(k.code, k.schema));
}

export function createSplatSorter(compute: ComputeContext, maxWorkgroups = 65535): SplatSorter {
  const histogram = compute.createBuffer(SplatKeyCount * 4, 'splat sort histogram');
  const offsets = compute.createBuffer(SplatKeyCount * 4, 'splat sort offsets');
  let keys = compute.createBuffer(4, 'splat sort keys');
  const keysValues = {
    uCount: ValueCell.create(0),
    uGroupsX: ValueCell.create(0),
    uRow: ValueCell.create(Vec4()),
    uDepthMin: ValueCell.create(0),
    uDepthScale: ValueCell.create(0),
    positions: ValueCell.create(keys),
    keys: ValueCell.create(keys),
    histogram: ValueCell.create(histogram),
  };
  const scatterValues = {
    uCount: ValueCell.create(0),
    uGroupsX: ValueCell.create(0),
    uOffset: ValueCell.create(0),
    keys: ValueCell.create(keys),
    offsets: ValueCell.create(offsets),
    order: ValueCell.create(keys),
  };
  const keysKernel: ComputeKernel = compute.createKernel(SplatKeysKernel, SplatKeysSchema, keysValues);
  const scatterKernel: ComputeKernel = compute.createKernel(SplatScatterKernel, SplatScatterSchema, scatterValues);
  const prefixSum: PrefixSum = createPrefixSum(compute, maxWorkgroups);
  const centersValues = {
    uCount: ValueCell.create(0),
    uGroupsX: ValueCell.create(0),
    uMeansMin: ValueCell.create(Vec3()),
    uMeansMax: ValueCell.create(Vec3()),
    tMeansLow: ValueCell.create(createNullTexture()),
    tMeansHigh: ValueCell.create(createNullTexture()),
    positions: ValueCell.create(keys),
  };
  // made with its first SOG: textures of none can't be bound
  let centersKernel: ComputeKernel | undefined;

  return {
    sort: (positions, order, jobs) => {
      for (const job of jobs) {
        const { count } = job;
        if (count === 0) continue;
        if (keys.byteLength < count * 4) {
          keys.destroy();
          keys = compute.createBuffer(count * 4, 'splat sort keys');
          ValueCell.update(keysValues.keys, keys);
          ValueCell.update(scatterValues.keys, keys);
        }
        const { x, y } = getDispatchSize(Math.ceil(count / WorkgroupSize), maxWorkgroups);
        const range = job.depthMax - job.depthMin;

        ValueCell.updateIfChanged(keysValues.uCount, count);
        ValueCell.updateIfChanged(keysValues.uGroupsX, x);
        ValueCell.update(keysValues.uRow, Vec4.set(keysValues.uRow.ref.value, job.a, job.b, job.c, 0));
        ValueCell.updateIfChanged(keysValues.uDepthMin, job.depthMin);
        ValueCell.updateIfChanged(keysValues.uDepthScale, range > 0 ? (SplatKeyCount - 1) / range : 0);
        if (keysValues.positions.ref.value !== positions) ValueCell.update(keysValues.positions, positions);

        ValueCell.updateIfChanged(scatterValues.uCount, count);
        ValueCell.updateIfChanged(scatterValues.uGroupsX, x);
        ValueCell.updateIfChanged(scatterValues.uOffset, job.offset);
        if (scatterValues.order.ref.value !== order) ValueCell.update(scatterValues.order, order);

        // a clear ends a compute pass: each instance a pass of its own
        histogram.clear();
        compute.beginComputePass({ label: 'splat sort' });
        keysKernel.dispatch(x, y);
        // 65536 keys scan in two levels (256 blocks, one): no copy, the pass stays open
        prefixSum.run(histogram, offsets, SplatKeyCount);
        scatterKernel.dispatch(x, y);
        compute.endComputePass();
      }
    },
    decodeSogCenters: (meansLow, meansHigh, meansMin, meansMax, count, positions) => {
      if (count === 0) return;
      const { x, y } = getDispatchSize(Math.ceil(count / WorkgroupSize), maxWorkgroups);
      ValueCell.updateIfChanged(centersValues.uCount, count);
      ValueCell.updateIfChanged(centersValues.uGroupsX, x);
      ValueCell.update(centersValues.uMeansMin, Vec3.copy(centersValues.uMeansMin.ref.value, meansMin));
      ValueCell.update(centersValues.uMeansMax, Vec3.copy(centersValues.uMeansMax.ref.value, meansMax));
      if (centersValues.tMeansLow.ref.value !== meansLow) ValueCell.update(centersValues.tMeansLow, meansLow);
      if (centersValues.tMeansHigh.ref.value !== meansHigh) ValueCell.update(centersValues.tMeansHigh, meansHigh);
      if (centersValues.positions.ref.value !== positions) ValueCell.update(centersValues.positions, positions);
      centersKernel ??= compute.createKernel(SplatCentersKernel, SplatCentersSchema, centersValues);
      compute.beginComputePass({ label: 'splat centers' });
      centersKernel.dispatch(x, y);
      compute.endComputePass();
    },
    destroy: () => {
      centersKernel?.destroy();
      keysKernel.destroy();
      scatterKernel.destroy();
      prefixSum.destroy();
      histogram.destroy();
      offsets.destroy();
      keys.destroy();
    },
  };
}
