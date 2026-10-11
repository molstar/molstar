/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { alignTo } from '../uniform-layout.js';
import { BufferUsage } from './constants.js';

const ChunkSize = 1 << 20;

type Chunk = { buffer: GPUBuffer; data: Uint8Array; used: number };

export type UniformSlot = { readonly buffer: GPUBuffer; readonly offset: number };

/**
 * Snapshots of uniform block data, bound with dynamic offsets. All writes to GPU buffers happen
 * before the commands of a submission run, so values that change between the draws of one
 * submission (e.g. the global block between passes) need a snapshot per change.
 */
export interface UniformRing {
  /** Number of the current submission, snapshots of older epochs are gone */
  readonly epoch: number;
  /** Copies `data` into the ring */
  write: (data: Uint8Array) => UniformSlot;
  /** Uploads the snapshots of the current submission, call before submitting */
  upload: () => void;
  /** Starts the next submission, call after submitting */
  reset: () => void;
  /** Drops the buffers of a lost device, snapshots are written anew */
  restore: () => void;
  destroy: () => void;
}

/** `ref.device` is read on use, it changes when the context restores a lost device */
export function createUniformRing(ref: { readonly device: GPUDevice }): UniformRing {
  const chunks: Chunk[] = [];
  let current = 0;
  let epoch = 0;

  function createChunk(size: number): Chunk {
    const buffer = ref.device.createBuffer({ size, usage: BufferUsage.UNIFORM | BufferUsage.COPY_DST });
    return { buffer, data: new Uint8Array(size), used: 0 };
  }

  return {
    get epoch() {
      return epoch;
    },
    write: (data: Uint8Array) => {
      const alignment = ref.device.limits.minUniformBufferOffsetAlignment;
      let chunk = chunks[current];
      if (chunk && alignTo(chunk.used, alignment) + data.byteLength > chunk.data.byteLength) {
        current += 1;
        chunk = chunks[current];
      }
      if (!chunk) {
        chunk = createChunk(Math.max(ChunkSize, alignTo(data.byteLength, alignment)));
        chunks[current] = chunk;
      }
      const offset = alignTo(chunk.used, alignment);
      chunk.data.set(data, offset);
      chunk.used = offset + data.byteLength;
      return { buffer: chunk.buffer, offset };
    },
    upload: () => {
      for (let i = 0; i <= current && i < chunks.length; ++i) {
        const { buffer, data, used } = chunks[i];
        if (used > 0) ref.device.queue.writeBuffer(buffer, 0, data, 0, alignTo(used, 4));
      }
    },
    reset: () => {
      for (const chunk of chunks) chunk.used = 0;
      current = 0;
      epoch += 1;
    },
    restore: () => {
      chunks.length = 0;
      current = 0;
      epoch += 1;
    },
    destroy: () => {
      for (const chunk of chunks) chunk.buffer.destroy();
      chunks.length = 0;
    },
  };
}

/**
 * CPU-side data of a uniform block with a snapshot in the ring that is renewed when the data
 * changed or the submission changed
 */
export interface UniformSnapshot {
  readonly bytes: Uint8Array;
  markChanged: () => void;
  /** The current snapshot, written to the ring first if needed */
  getSlot: () => UniformSlot;
}

export function createUniformSnapshot(ring: UniformRing, size: number): UniformSnapshot {
  const bytes = new Uint8Array(size);
  let changed = true;
  let epoch = -1;
  let slot: UniformSlot | undefined;
  return {
    bytes,
    markChanged: () => {
      changed = true;
    },
    getSlot: () => {
      if (changed || epoch !== ring.epoch || !slot) {
        slot = ring.write(bytes);
        epoch = ring.epoch;
        changed = false;
      }
      return slot;
    },
  };
}
