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
  const device = {
    createBuffer: (d: { size: number }) => {
      const buffer = { size: d.size, destroy: vi.fn() };
      buffers.push(buffer);
      return buffer;
    },
    queue: { writeBuffer: vi.fn() },
  } as unknown as GPUDevice;
  return { device, buffers };
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
