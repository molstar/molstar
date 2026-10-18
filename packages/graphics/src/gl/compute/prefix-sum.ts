/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type StorageBuffer,
  StorageBufferSpec,
  type KernelRequest,
} from '../compute.js';
import { UniformSpec } from '../renderable/schema.js';

/** Elements per workgroup of the scan kernels */
const BlockSize = 256;

/** Exclusive scan of each block of `input` into `output`, the block totals into `blockSums` */
export const ScanBlockKernel = KernelCode('scan-block', [BlockSize]);
/** Adds the scanned block totals to the elements of each block */
export const ScanAddKernel = KernelCode('scan-add', [BlockSize]);

export const ScanBlockSchema = {
  uCount: UniformSpec('i'),
  uBlockCount: UniformSpec('i'),
  uGroupsX: UniformSpec('i'),
  input: StorageBufferSpec('u32', 'read'),
  output: StorageBufferSpec('u32', 'read-write'),
  blockSums: StorageBufferSpec('u32', 'read-write'),
};

export const ScanAddSchema = {
  uCount: UniformSpec('i'),
  uBlockCount: UniformSpec('i'),
  uGroupsX: UniformSpec('i'),
  output: StorageBufferSpec('u32', 'read-write'),
  blockSums: StorageBufferSpec('u32', 'read'),
};

/** The kernels of a prefix sum, see `ComputeContext.prepareKernels` */
export const PrefixSumKernels: readonly KernelRequest[] = [
  { code: ScanBlockKernel, schema: ScanBlockSchema },
  { code: ScanAddKernel, schema: ScanAddSchema },
];

/** Workgroup counts for `count` workgroups, in 2d beyond the per-dimension limit; kernels linearize the id */
export function getDispatchSize(count: number, maxPerDimension: number) {
  if (count <= maxPerDimension) return { x: count, y: 1 };
  return { x: maxPerDimension, y: Math.ceil(count / maxPerDimension) };
}

/** Level buffers reused across runs of a prefix sum, by level */
type Level = { sums: StorageBuffer; scanned: StorageBuffer };

export interface PrefixSum {
  /**
   * Exclusive prefix sum of the first `count` u32 of `input` into `output` (in the open compute pass, which it ends
   * when it copies); returns a buffer whose first u32 is the total
   */
  run: (input: StorageBuffer, output: StorageBuffer, count: number) => StorageBuffer;
  destroy: () => void;
}

export function isPrefixSumSupported(compute: ComputeContext) {
  return compute.isKernelSupported(ScanBlockKernel) && compute.isKernelSupported(ScanAddKernel);
}

/**
 * Work-efficient exclusive prefix sum (scan) of u32: each workgroup scans a block of 256 elements in workgroup memory,
 * the block totals are scanned recursively and added back. `maxWorkgroups` is the per-dimension dispatch limit.
 */
export function createPrefixSum(compute: ComputeContext, maxWorkgroups = 65535): PrefixSum {
  const empty = compute.createBuffer(4, 'prefix sum empty');
  const scanValues = {
    uCount: ValueCell.create(0),
    uBlockCount: ValueCell.create(0),
    uGroupsX: ValueCell.create(0),
    input: ValueCell.create(empty),
    output: ValueCell.create(empty),
    blockSums: ValueCell.create(empty),
  };
  const addValues = {
    uCount: ValueCell.create(0),
    uBlockCount: ValueCell.create(0),
    uGroupsX: ValueCell.create(0),
    output: ValueCell.create(empty),
    blockSums: ValueCell.create(empty),
  };
  const scanKernel: ComputeKernel = compute.createKernel(ScanBlockKernel, ScanBlockSchema, scanValues);
  const addKernel: ComputeKernel = compute.createKernel(ScanAddKernel, ScanAddSchema, addValues);
  const levels: Level[] = [];

  function getLevel(index: number, blockCount: number) {
    let level = levels[index];
    if (!level || level.sums.byteLength < blockCount * 4) {
      level?.sums.destroy();
      level?.scanned.destroy();
      level = {
        sums: compute.createBuffer(blockCount * 4, `prefix sum level ${index}`),
        scanned: compute.createBuffer(blockCount * 4, `prefix sum level ${index} scanned`),
      };
      levels[index] = level;
    }
    return level;
  }

  function scan(input: StorageBuffer, output: StorageBuffer, count: number, levelIndex: number): StorageBuffer {
    const blockCount = Math.ceil(count / BlockSize);
    const { x, y } = getDispatchSize(blockCount, maxWorkgroups);
    const level = getLevel(levelIndex, blockCount);

    ValueCell.updateIfChanged(scanValues.uCount, count);
    ValueCell.updateIfChanged(scanValues.uBlockCount, blockCount);
    ValueCell.updateIfChanged(scanValues.uGroupsX, x);
    ValueCell.update(scanValues.input, input);
    ValueCell.update(scanValues.output, output);
    ValueCell.update(scanValues.blockSums, level.sums);
    scanKernel.dispatch(x, y);
    // the total of a single block is the total
    if (blockCount === 1) return level.sums;

    const total = scan(level.sums, level.scanned, blockCount, levelIndex + 1);

    ValueCell.updateIfChanged(addValues.uCount, count);
    ValueCell.updateIfChanged(addValues.uBlockCount, blockCount);
    ValueCell.updateIfChanged(addValues.uGroupsX, x);
    ValueCell.update(addValues.output, output);
    ValueCell.update(addValues.blockSums, level.scanned);
    addKernel.dispatch(x, y);
    return total;
  }

  return {
    run: (input, output, count) => {
      if (count <= 0) throw new Error('prefix sum of no elements');
      return scan(input, output, count, 0);
    },
    destroy: () => {
      for (const l of levels) {
        l.sums.destroy();
        l.scanned.destroy();
      }
      levels.length = 0;
      scanKernel.destroy();
      addKernel.destroy();
      empty.destroy();
    },
  };
}
