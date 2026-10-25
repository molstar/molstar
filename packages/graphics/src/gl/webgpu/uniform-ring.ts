/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { alignTo } from '../uniform-layout.js';
import type { GraphicsStats } from '../stats.js';
import type { BufferPages, BufferRange } from './buffer-pages.js';
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
  readonly device: GPUDevice;
  /** bytes of the ring's buffers and of the arenas created with it */
  readonly byteCounts: { uniform: number };
  /** counts of the ring's buffers and of the arenas' (`resourceCounts.uniform`) */
  readonly stats?: GraphicsStats;
  /**
   * Buffers of the arenas created with it replaced so far (growth, the old one released): bind groups recorded before
   * (render bundles) may bind a released buffer
   */
  arenaReplacements: number;
  /** Copies `data` into the ring */
  write: (data: Uint8Array) => UniformSlot;
  /** Uploads the snapshots of the current submission, call before submitting */
  upload: () => void;
  /** Starts the next submission, call after submitting */
  reset: () => void;
  destroy: () => void;
}

export function createUniformRing(
  device: GPUDevice,
  byteCounts: { uniform: number } = { uniform: 0 },
  stats?: GraphicsStats,
): UniformRing {
  const chunks: Chunk[] = [];
  let current = 0;
  let epoch = 0;

  function createChunk(size: number): Chunk {
    const buffer = device.createBuffer({ size, usage: BufferUsage.UNIFORM | BufferUsage.COPY_DST });
    byteCounts.uniform += size;
    if (stats) stats.resourceCounts.uniform += 1;
    return { buffer, data: new Uint8Array(size), used: 0 };
  }

  return {
    get epoch() {
      return epoch;
    },
    device,
    byteCounts,
    arenaReplacements: 0,
    stats,
    write: (data: Uint8Array) => {
      const alignment = device.limits.minUniformBufferOffsetAlignment;
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
        if (used > 0) device.queue.writeBuffer(buffer, 0, data, 0, alignTo(used, 4));
      }
    },
    reset: () => {
      for (const chunk of chunks) chunk.used = 0;
      current = 0;
      epoch += 1;
    },
    destroy: () => {
      for (const chunk of chunks) {
        chunk.buffer.destroy();
        byteCounts.uniform -= chunk.data.byteLength;
        if (stats) stats.resourceCounts.uniform -= 1;
      }
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

/**
 * Uniform blocks of one owner (a render item, a shared material) in one buffer, a slot each, so that one bind group
 * with dynamic offsets serves all of them. Slots keep their data across submissions and are written with `writeBuffer`
 * when it changed. `writeBuffer` runs before the commands of the current submission, so a change after a draw of the
 * current submission read the slot goes to a ring snapshot for the rest of the submission, and to the slot in the next.
 * The buffer grows by replacing it (the old one is released after the submission, slots are written anew).
 *
 * With `pages`, the slots are ranges of shared uniform pages instead (`WebGPUDrawContext.uniformPages`), freed with
 * the arena: the owners' blocks in few buffers, so that their bind groups can be shared (`render-item.ts`)
 */
export interface UniformArena {
  createSnapshot: (size: number) => UniformSnapshot;
  destroy: () => void;
}

export function createUniformArena(
  ring: UniformRing,
  release: (buffer: GPUBuffer) => void,
  pages?: BufferPages,
): UniformArena {
  let buffer: GPUBuffer | undefined;
  const ranges: BufferRange[] = [];
  let capacity = 0;
  let used = 0;
  /** incremented when the buffer is replaced, slots of older generations are written anew */
  let generation = 0;

  function getBuffer() {
    if (!buffer || capacity < used) {
      if (buffer) {
        release(buffer);
        ring.arenaReplacements += 1;
      }
      const previous = buffer ? capacity : 0;
      capacity = Math.max(used, capacity * 2, 1024);
      buffer = ring.device.createBuffer({ size: capacity, usage: BufferUsage.UNIFORM | BufferUsage.COPY_DST });
      ring.byteCounts.uniform += capacity - previous;
      // a replaced buffer is released (above), so one buffer per arena
      if (ring.stats && !previous) ring.stats.resourceCounts.uniform += 1;
      generation += 1;
    }
    return buffer;
  }

  return {
    createSnapshot: (size: number) => {
      const range = pages?.allocate(alignTo(size, 16));
      let offset: number;
      if (range) {
        ranges.push(range);
        ring.byteCounts.uniform += range.size;
        offset = range.offset;
      } else {
        offset = alignTo(used, ring.device.limits.minUniformBufferOffsetAlignment);
        used = offset + alignTo(size, 16);
      }
      const bytes = new Uint8Array(size);
      let changed = true;
      /** the slot needs writing, after a change that went to the ring */
      let pending = false;
      let writtenGeneration = -1;
      let usedEpoch = -1;
      let usedGeneration = -1;
      let ringSlot: UniformSlot | undefined;
      let ringEpoch = -1;
      let slot: UniformSlot | undefined;
      return {
        bytes,
        markChanged: () => {
          changed = true;
        },
        getSlot: () => {
          const b = range ? range.buffer : getBuffer();
          const epoch = ring.epoch;
          if (ringEpoch === epoch && ringSlot) {
            if (changed) {
              ringSlot = ring.write(bytes);
              changed = false;
            }
            return ringSlot;
          }
          if (changed || pending || writtenGeneration !== generation) {
            if (usedEpoch === epoch && usedGeneration === generation) {
              ringSlot = ring.write(bytes);
              ringEpoch = epoch;
              changed = false;
              pending = true;
              return ringSlot;
            }
            ring.device.queue.writeBuffer(b, offset, bytes, 0, alignTo(size, 4));
            changed = false;
            pending = false;
            writtenGeneration = generation;
          }
          usedEpoch = epoch;
          usedGeneration = generation;
          if (slot?.buffer !== b) slot = { buffer: b, offset };
          return slot;
        },
      };
    },
    destroy: () => {
      for (const range of ranges) {
        pages!.free(range);
        ring.byteCounts.uniform -= range.size;
      }
      ranges.length = 0;
      if (buffer) {
        release(buffer);
        ring.byteCounts.uniform -= capacity;
        if (ring.stats) ring.stats.resourceCounts.uniform -= 1;
      }
      buffer = undefined;
    },
  };
}
