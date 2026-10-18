/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import {
  alignTo,
  createUniformBlockData,
  createUniformBlockWriter,
  UniformBlockInfo,
  type UniformBlockLayout,
} from '../uniform-layout.js';
import type { WebGLState } from './state.js';
import type { UniformBuffer, UniformsList } from '../uniform.js';

const getNextUniformSlabId = idFactory();
const getNextUniformBufferId = idFactory();

/** GPU buffer with a CPU mirror, uploads cover the dirty byte range */
export interface UniformSlab {
  readonly id: number;
  readonly data: ArrayBuffer;
  buffer: WebGLBuffer;
  dirtyStart: number;
  dirtyEnd: number;
}

export interface UniformSlice {
  readonly slab: UniformSlab;
  readonly offset: number;
  readonly size: number;
}

/**
 * Sub-allocates uniform block instances from a few large buffers, bound as ranges
 * . Slices are grouped into power-of-two size
 * classes, aligned to `UNIFORM_BUFFER_OFFSET_ALIGNMENT`.
 */
export interface UniformArena {
  readonly alignment: number;
  /** zero-initialized slice of at least `size` bytes */
  allocate: (size: number) => UniformSlice;
  free: (slice: UniformSlice) => void;
  /** dedicated slab, e.g. for ring buffers */
  createSlab: (byteSize: number) => UniformSlab;
  releaseSlab: (slab: UniformSlab) => void;
  markDirty: (slab: UniformSlab, start: number, end: number) => void;
  /** upload the dirty range of `slab` */
  flush: (slab: UniformSlab) => void;
  bind: (index: number, slab: UniformSlab, offset: number, size: number) => void;

  getByteCount: () => number;
  destroy: () => void;
}

const SlabByteSize = 64 * 1024;

function nextPowerOfTwo(n: number) {
  return Math.pow(2, Math.ceil(Math.log2(n)));
}

export function createUniformArena(gl: WebGL2RenderingContext, state: WebGLState): UniformArena {
  const alignment: number = gl.getParameter(gl.UNIFORM_BUFFER_OFFSET_ALIGNMENT);
  const slabs = new Set<UniformSlab>();
  const freeSlices = new Map<number, UniformSlice[]>();

  function getGlBuffer(byteSize: number) {
    const buffer = gl.createBuffer();
    if (buffer === null) throw new Error('Could not create WebGL buffer');
    gl.bindBuffer(gl.UNIFORM_BUFFER, buffer);
    gl.bufferData(gl.UNIFORM_BUFFER, byteSize, gl.DYNAMIC_DRAW);
    return buffer;
  }

  function createSlab(byteSize: number): UniformSlab {
    const slab: UniformSlab = {
      id: getNextUniformSlabId(),
      data: new ArrayBuffer(byteSize),
      buffer: getGlBuffer(byteSize),
      dirtyStart: byteSize,
      dirtyEnd: 0,
    };
    slabs.add(slab);
    return slab;
  }

  function markDirty(slab: UniformSlab, start: number, end: number) {
    if (start < slab.dirtyStart) slab.dirtyStart = start;
    if (end > slab.dirtyEnd) slab.dirtyEnd = end;
  }

  return {
    alignment,
    allocate: (size: number) => {
      const slotSize = alignTo(nextPowerOfTwo(Math.max(size, alignment)), alignment);
      let free = freeSlices.get(slotSize);
      if (free === undefined) {
        free = [];
        freeSlices.set(slotSize, free);
      }
      if (free.length === 0) {
        const slab = createSlab(Math.max(SlabByteSize, slotSize));
        // reversed so that slices are handed out in ascending offset order
        for (let o = slab.data.byteLength - slotSize; o >= 0; o -= slotSize) {
          free.push({ slab, offset: o, size: slotSize });
        }
      }
      const slice = free.pop()!;
      new Uint8Array(slice.slab.data, slice.offset, slotSize).fill(0);
      markDirty(slice.slab, slice.offset, slice.offset + slotSize);
      return slice;
    },
    free: (slice: UniformSlice) => {
      freeSlices.get(slice.size)!.push(slice);
    },
    createSlab,
    releaseSlab: (slab: UniformSlab) => {
      gl.deleteBuffer(slab.buffer);
      slabs.delete(slab);
    },
    markDirty,
    flush: (slab: UniformSlab) => {
      if (slab.dirtyEnd <= slab.dirtyStart) return;
      gl.bindBuffer(gl.UNIFORM_BUFFER, slab.buffer);
      gl.bufferSubData(
        gl.UNIFORM_BUFFER,
        slab.dirtyStart,
        new Uint8Array(slab.data),
        slab.dirtyStart,
        slab.dirtyEnd - slab.dirtyStart,
      );
      slab.dirtyStart = slab.data.byteLength;
      slab.dirtyEnd = 0;
    },
    bind: (index: number, slab: UniformSlab, offset: number, size: number) => {
      state.bindUniformBufferRange(index, slab.id, slab.buffer, offset, size);
    },

    getByteCount: () => {
      let count = 0;
      slabs.forEach((s) => (count += s.data.byteLength));
      return count;
    },
    destroy: () => {
      slabs.forEach((slab) => gl.deleteBuffer(slab.buffer));
      slabs.clear();
      freeSlices.clear();
    },
  };
}

//

export function createUniformBuffer(arena: UniformArena, layout: UniformBlockLayout): UniformBuffer {
  const id = getNextUniformBufferId();
  const { size } = layout;
  const { bindingIndex } = UniformBlockInfo[layout.block];

  const slice = arena.allocate(size);
  const { slab, offset } = slice;
  const data = createUniformBlockData(slab.data, offset, size);
  const markDirty = () => arena.markDirty(slab, offset, offset + size);

  let destroyed = false;

  return {
    id,

    getByteCount: () => size,
    createWriter: (uniformValues: UniformsList) => createUniformBlockWriter(layout, data, uniformValues, markDirty),
    sync: () => arena.flush(slab),
    bind: () => arena.bind(bindingIndex, slab, offset, slice.size),

    destroy: () => {
      if (destroyed) return;
      arena.free(slice);
      destroyed = true;
    },
  };
}

/**
 * Every `sync` after a change writes a new snapshot into the next slot instead of
 * updating data that previous draws use; wrapping around overwrites the oldest slot,
 * which is safe as WebGL commands are ordered.
 */
export function createUniformRingBuffer(
  arena: UniformArena,
  layout: UniformBlockLayout,
  slotCount = 64,
): UniformBuffer {
  const id = getNextUniformBufferId();
  const { size } = layout;
  const { bindingIndex } = UniformBlockInfo[layout.block];

  const slotSize = alignTo(size, arena.alignment);
  const slab = arena.createSlab(slotSize * slotCount);
  const slabU8 = new Uint8Array(slab.data);

  const staging = new ArrayBuffer(size);
  const stagingU8 = new Uint8Array(staging);
  const data = createUniformBlockData(staging, 0, size);

  let dirty = true;
  let slot = -1;
  let destroyed = false;

  function sync() {
    if (!dirty) return;
    slot = (slot + 1) % slotCount;
    const offset = slot * slotSize;
    slabU8.set(stagingU8, offset);
    arena.markDirty(slab, offset, offset + size);
    arena.flush(slab);
    dirty = false;
  }

  return {
    id,

    getByteCount: () => size,
    createWriter: (uniformValues: UniformsList) =>
      createUniformBlockWriter(layout, data, uniformValues, () => {
        dirty = true;
      }),
    sync,
    bind: () => {
      if (slot === -1) sync();
      arena.bind(bindingIndex, slab, slot * slotSize, slotSize);
    },

    destroy: () => {
      if (destroyed) return;
      arena.releaseSlab(slab);
      destroyed = true;
    },
  };
}
