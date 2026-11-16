/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { createRangeAllocator } from '../range-allocator.js';
import { createBufferPages } from '../webgpu/buffer-pages.js';

describe('range allocator', () => {
  it('allocates aligned ranges first fit, until full', () => {
    const a = createRangeAllocator(64, 16);
    expect(a.allocate(1)).toBe(0);
    expect(a.allocate(17)).toBe(16);
    expect(a.allocate(16)).toBe(48);
    expect(a.used).toBe(64);
    expect(a.allocate(1)).toBe(-1);
  });

  it('merges freed ranges with their free neighbors', () => {
    const a = createRangeAllocator(64, 16);
    const offsets = [a.allocate(16), a.allocate(16), a.allocate(16), a.allocate(16)];
    expect(offsets).toEqual([0, 16, 32, 48]);
    a.free(16, 16);
    a.free(48, 16);
    // not adjacent: no 32 bytes in one piece
    expect(a.allocate(32)).toBe(-1);
    a.free(32, 16);
    // 16..64 merged
    expect(a.allocate(48)).toBe(16);
    a.free(16, 48);
    a.free(0, 16);
    expect(a.used).toBe(0);
    expect(a.allocate(64)).toBe(0);
  });

  it('grows: the units past its capacity free, merged with a free range ending there', () => {
    const a = createRangeAllocator(64, 16);
    expect(a.allocate(48)).toBe(0);
    // 48..64 free: grown, 48..128 one range
    a.grow(128);
    expect(a.capacity).toBe(128);
    expect(a.allocate(80)).toBe(48);
    expect(a.allocate(1)).toBe(-1);
    // full: grown, a range of its own
    a.grow(160);
    expect(a.allocate(32)).toBe(128);
    // not smaller
    a.grow(16);
    expect(a.capacity).toBe(160);
  });

  it('reuses the first fitting gap', () => {
    const a = createRangeAllocator(128, 16);
    for (let i = 0; i < 8; ++i) a.allocate(16);
    a.free(32, 32);
    a.free(96, 16);
    expect(a.allocate(16)).toBe(32);
    expect(a.allocate(16)).toBe(48);
    expect(a.allocate(16)).toBe(96);
    expect(a.allocate(16)).toBe(-1);
  });
});

// Freed ranges of the WebGPU buffer pages: reusable after the next submission, or at once where nothing recorded is
// left to submit (an idle canvas); a page left empty dropped, its buffer destroyed, without waiting for an allocation.

/** a device of buffers that record their destroy */
function fakeDevice() {
  const buffers: { size: number; destroy: ReturnType<typeof vi.fn> }[] = [];
  const encoders: { copyBufferToBuffer: ReturnType<typeof vi.fn> }[] = [];
  const device = {
    createBuffer: (d: { size: number }) => {
      const buffer = { size: d.size, destroy: vi.fn() };
      buffers.push(buffer);
      return buffer;
    },
    queue: { writeBuffer: vi.fn(), submit: vi.fn() },
    createCommandEncoder: () => {
      const encoder = { copyBufferToBuffer: vi.fn(), finish: () => encoder };
      encoders.push(encoder);
      return encoder;
    },
  } as unknown as GPUDevice;
  return { device, buffers, encoders };
}

describe('WebGPU buffer pages, freed ranges', () => {
  const pageSize = 1024;

  it('drop a page left empty after the next submission, without an allocation', () => {
    const { device, buffers } = fakeDevice();
    let submit = 0;
    const onFree = vi.fn();
    const pages = createBufferPages(device, {
      usage: 0,
      pageSize,
      label: 'test',
      getSubmitCount: () => submit,
      onFree,
    });
    // larger than a quarter page: a page of its own
    const range = pages.allocate(pageSize);
    expect(pages.pageCount).toBe(1);
    pages.free(range);
    expect(onFree).toHaveBeenCalledTimes(1);

    // the commands of this submission may still read it
    pages.reclaim();
    expect(pages.pageCount).toBe(1);
    expect(buffers[0].destroy).not.toHaveBeenCalled();

    submit += 1;
    pages.reclaim();
    expect(pages.pageCount).toBe(0);
    expect(pages.pageBytes).toBe(0);
    expect(buffers[0].destroy).toHaveBeenCalledTimes(1);
  });

  it('drop it at once where nothing recorded is left to submit (idle)', () => {
    const { device, buffers } = fakeDevice();
    const pages = createBufferPages(device, { usage: 0, pageSize, label: 'test', getSubmitCount: () => 0 });
    pages.free(pages.allocate(pageSize));
    pages.reclaim(true);
    expect(pages.pageCount).toBe(0);
    expect(buffers[0].destroy).toHaveBeenCalledTimes(1);
  });
});

// Copies between ranges on the GPU (a render item's layers moved to a larger range, not uploaded again): submitted at
// once, so before the writes after; through a buffer of their own where a buffer is read and written
describe('WebGPU buffer pages, copies', () => {
  const pageSize = 1024;

  it('copy between buffers directly, in one submission', () => {
    const { device, buffers, encoders } = fakeDevice();
    const pages = createBufferPages(device, { usage: 0, pageSize, label: 'test', getSubmitCount: () => 0 });
    // larger than a quarter page: each a page of its own
    const from = pages.allocate(512);
    const to = pages.allocate(768);
    pages.copy([{ from, fromOffset: 16, to, toOffset: 32, size: 256 }]);
    expect(encoders).toHaveLength(1);
    expect(encoders[0].copyBufferToBuffer.mock.calls).toEqual([
      [from.buffer, from.offset + 16, to.buffer, to.offset + 32, 256],
    ]);
    expect((device.queue.submit as ReturnType<typeof vi.fn>).mock.calls).toHaveLength(1);
    expect(buffers).toHaveLength(2);
    // nothing to copy: no submission
    pages.copy([]);
    expect((device.queue.submit as ReturnType<typeof vi.fn>).mock.calls).toHaveLength(1);
  });

  it('copy through a staging buffer where a buffer is read and written: all read first, then written', () => {
    const { device, buffers, encoders } = fakeDevice();
    const pages = createBufferPages(device, { usage: 0, pageSize, label: 'test', getSubmitCount: () => 0 });
    // two ranges of a page: a layer moved behind another that grew, which moves too (overlapping its old place)
    const a = pages.allocate(64);
    const b = pages.allocate(64);
    expect(a.buffer).toBe(b.buffer);
    pages.copy([
      { from: a, fromOffset: 0, to: a, toOffset: 32, size: 32 },
      { from: a, fromOffset: 32, to: b, toOffset: 0, size: 16 },
    ]);
    const staging = buffers[1];
    expect(staging.size).toBe(48);
    expect(encoders[0].copyBufferToBuffer.mock.calls).toEqual([
      [a.buffer, a.offset, staging, 0, 32],
      [a.buffer, a.offset + 32, staging, 32, 16],
      [staging, 0, a.buffer, a.offset + 32, 32],
      [staging, 32, b.buffer, b.offset, 16],
    ]);
    expect(staging.destroy).toHaveBeenCalledTimes(1);
  });
});
