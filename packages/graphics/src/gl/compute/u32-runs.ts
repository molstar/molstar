/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import {
  type ComputeContext,
  type ComputeKernel,
  KernelCode,
  type KernelRequest,
  type StorageBuffer,
  StorageBufferSpec,
} from '../compute.js';
import { UniformSpec } from '../renderable/schema.js';

/**
 * Runs of u32 words written on the GPU: a run's `count` words from `word`, `(start + step i) | bits`. What a culled
 * object's data holds per instance that is serial or the same for a run of them (its identity order, slot words,
 * `aInstance`, segment indices: 4 B each an instance, 500 MB the yeast cell's 31 M) made where it is read, not uploaded
 */
export const U32RunsKernel = KernelCode('u32-runs', [64]);

/** vec4s of a run's record (`U32RunsSchema.runs`) */
export const U32RunRecordVec4s = 2;

export const U32RunsSchema = {
  uThreadCount: UniformSpec('i'),
  uRunCount: UniformSpec('i'),
  /** per run: its first word, count, start and step; its bits and first thread */
  runs: StorageBufferSpec('vec4<u32>', 'read'),
  dst: StorageBufferSpec('u32', 'read-write'),
};

export const U32RunsKernels: readonly KernelRequest[] = [{ code: U32RunsKernel, schema: U32RunsSchema }];

/** `count` words from `word`: `(start + step i) | bits` (u32, wrapping) */
export interface U32Run {
  readonly word: number;
  readonly count: number;
  readonly start: number;
  readonly step: number;
  readonly bits: number;
}

/** Writes runs of u32 words into a storage buffer on the GPU */
export interface U32RunWriter {
  /**
   * the runs into `dst`, submitted (later writes land after); `dispatchItems` a kernel over a count of items (in rows
   * past the limit). `submit` false: left to the next submit (a frame's: a submit costs), the caller submitting before
   * a write of its own to `dst`'s words would land ahead of it
   */
  write: (
    dst: StorageBuffer,
    runs: readonly U32Run[],
    dispatchItems: (kernel: ComputeKernel, n: number) => void,
    submit?: boolean,
  ) => void;
  destroy: () => void;
}

/** The runs' records and thread count, runs without words left out */
export function getU32RunRecords(runs: readonly U32Run[]) {
  const records = new Uint32Array(runs.length * U32RunRecordVec4s * 4);
  let threads = 0;
  let count = 0;
  for (const r of runs) {
    if (r.count <= 0) continue;
    const o = count * U32RunRecordVec4s * 4;
    records[o] = r.word;
    records[o + 1] = r.count;
    records[o + 2] = r.start >>> 0;
    records[o + 3] = r.step >>> 0;
    records[o + 4] = r.bits >>> 0;
    records[o + 5] = threads;
    threads += r.count;
    count += 1;
  }
  return { records: records.subarray(0, count * U32RunRecordVec4s * 4), runs: count, threads };
}

/**
 * The writer's kernel created at once, as the poses' expander's (`createCulledPoseExpander`): the data it writes is
 * complete before the first draw
 */
export function createU32RunWriter(compute: ComputeContext): U32RunWriter {
  const placeholder = compute.createBuffer(16, 'u32 runs placeholder');
  const values = {
    uThreadCount: ValueCell.create(0),
    uRunCount: ValueCell.create(0),
    runs: ValueCell.create(placeholder),
    dst: ValueCell.create(placeholder),
  };
  const kernel = compute.createKernel(U32RunsKernel, U32RunsSchema, values);
  return {
    write: (dst, runs, dispatchItems, submit = true) => {
      const { records, runs: runCount, threads } = getU32RunRecords(runs);
      if (threads === 0) return;
      const buffer = compute.createBuffer(Math.max(16, records.byteLength), 'u32 runs');
      buffer.write(records);
      ValueCell.updateIfChanged(values.runs, buffer);
      ValueCell.updateIfChanged(values.dst, dst);
      ValueCell.updateIfChanged(values.uThreadCount, threads);
      ValueCell.updateIfChanged(values.uRunCount, runCount);
      compute.beginComputePass({ label: 'u32 runs' });
      dispatchItems(kernel, threads);
      compute.endComputePass();
      if (submit) compute.submit();
      // destroyed once the submitted work is done (the device keeps it until then; a destroy before a submit deferred)
      buffer.destroy();
      ValueCell.updateIfChanged(values.runs, placeholder);
      ValueCell.updateIfChanged(values.dst, placeholder);
    },
    destroy: () => {
      kernel.destroy();
      placeholder.destroy();
    },
  };
}
