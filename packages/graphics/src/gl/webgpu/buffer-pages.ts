/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BufferUsage } from './constants.js';
import { createRangeAllocator, type RangeAllocator } from '../range-allocator.js';

/** A range of a page: bound with its offset and size (`setVertexBuffer`, `setIndexBuffer`) */
export interface BufferRange {
  readonly buffer: GPUBuffer;
  readonly offset: number;
  /** the bytes reserved (aligned), at least the bytes written */
  readonly size: number;
}

/**
 * Shared GPU buffers (pages) that data of many owners is sub-allocated in: draws binding ranges of one buffer cost the
 * replay of a render bundle far less than buffers of their own (see `shared-draw-resources.md`). A range freed becomes
 * reusable after the next submission: commands recorded before may still read it, and a queue write runs before them
 */
export interface BufferPages {
  /** a range of at least `size` bytes */
  allocate(size: number): BufferRange;
  free(range: BufferRange): void;
  /**
   * writes `array` into `range` at `byteOffset` (a multiple of 4; the default 0); a last partial word is padded with
   * zeros
   */
  write(range: BufferRange, array: ArrayBufferView, byteOffset?: number): void;
  /** pages in use and their bytes */
  readonly pageCount: number;
  readonly pageBytes: number;
  destroy(): void;
}

/** Offsets and sizes of ranges: vertex and index offsets need 4 bytes, 16 keeps them on vector boundaries */
const RangeAlignment = 16;

export function createBufferPages(
  device: GPUDevice,
  options: {
    readonly usage: number;
    /** the size of a page; larger ranges get a page of their own */
    readonly pageSize: number;
    /** of offsets and sizes, `RangeAlignment` by default (uniform ranges need the dynamic offset alignment) */
    readonly alignment?: number;
    readonly label: string;
    /** advances with each submission (`GraphicsContext.submitCount`) */
    readonly getSubmitCount: () => number;
  },
): BufferPages {
  type Page = { readonly buffer: GPUBuffer; readonly allocator: RangeAllocator; readonly dedicated: boolean };
  type Range = BufferRange & { readonly page: Page };
  const pages: Page[] = [];
  /** ranges freed, reusable once a submission followed (`submit`: the count when freed) */
  const pending: { range: Range; submit: number }[] = [];
  let pageBytes = 0;
  const alignment = options.alignment ?? RangeAlignment;

  const createPage = (size: number, dedicated: boolean): Page => {
    const buffer = device.createBuffer({
      size,
      usage: options.usage | BufferUsage.COPY_DST,
      label: `${options.label} ${pages.length}`,
    });
    const page = { buffer, allocator: createRangeAllocator(size, alignment), dedicated };
    pages.push(page);
    pageBytes += size;
    return page;
  };

  /** the submission pending ranges were last released in: none freed since is reusable before the next */
  let releasedAt = -1;
  const releasePending = () => {
    if (pending.length === 0) return;
    const submit = options.getSubmitCount();
    if (submit === releasedAt) return;
    releasedAt = submit;
    let kept = 0;
    for (let i = 0, il = pending.length; i < il; ++i) {
      const p = pending[i];
      if (p.submit === submit) {
        pending[kept++] = p;
        continue;
      }
      const { page, offset, size } = p.range;
      page.allocator.free(offset, size);
      if (page.dedicated || page.allocator.used === 0) dropIfEmpty(page);
    }
    pending.length = kept;
  };

  const dropIfEmpty = (page: Page) => {
    if (page.allocator.used !== 0) return;
    // the last shared page stays for the next allocations
    if (!page.dedicated && pages.filter((p) => !p.dedicated).length <= 1) return;
    pages.splice(pages.indexOf(page), 1);
    pageBytes -= page.allocator.capacity;
    page.buffer.destroy();
  };

  return {
    allocate(size) {
      releasePending();
      const n = Math.ceil(Math.max(size, 4) / alignment) * alignment;
      let page: Page | undefined;
      let offset = -1;
      if (n > options.pageSize / 4) {
        page = createPage(n, true);
        offset = page.allocator.allocate(n);
      } else {
        for (const p of pages) {
          if (p.dedicated) continue;
          offset = p.allocator.allocate(n);
          if (offset >= 0) {
            page = p;
            break;
          }
        }
        if (!page) {
          page = createPage(options.pageSize, false);
          offset = page.allocator.allocate(n);
        }
      }
      const range: Range = { buffer: page.buffer, offset, size: n, page };
      return range;
    },
    free(range) {
      pending.push({ range: range as Range, submit: options.getSubmitCount() });
    },
    write(range, array, byteOffset = 0) {
      const offset = range.offset + byteOffset;
      const bytes = Math.floor(array.byteLength / 4) * 4;
      if (bytes > 0) device.queue.writeBuffer(range.buffer, offset, array.buffer, array.byteOffset, bytes);
      const rest = array.byteLength - bytes;
      if (rest > 0) {
        const tail = new Uint8Array(4);
        tail.set(new Uint8Array(array.buffer, array.byteOffset + bytes, rest));
        device.queue.writeBuffer(range.buffer, offset + bytes, tail);
      }
    },
    get pageCount() {
      return pages.length;
    },
    get pageBytes() {
      return pageBytes;
    },
    destroy() {
      for (const p of pages) p.buffer.destroy();
      pages.length = 0;
      pending.length = 0;
      pageBytes = 0;
    },
  };
}
