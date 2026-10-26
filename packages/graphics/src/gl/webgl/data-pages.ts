/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { DataFormat } from '../renderable/schema.js';
import type { DataArray } from '../renderable/util.js';
import { createRangeAllocator, type RangeAllocator } from '../range-allocator.js';
import type { GLTexture } from './texture.js';
import { getDataTextureSize } from './data-textures.js';

/**
 * WebGL's pages of data arrays (`DataSpec`), as WebGPU's storage pages: textures of one format shared by the data of
 * many items, each in a range of elements (texels, row by row), its offset an object uniform (`u<Name>Offset`). Draws of
 * items in the same pages bind the same textures, which the programs skip (`WebGLState.boundTextures`)
 */

/** A range of elements of a page: the texture, of `width` × `height` texels, and the elements from `offset` */
export interface DataTextureRange {
  readonly texture: GLTexture;
  readonly width: number;
  readonly height: number;
  readonly offset: number;
  /** the elements reserved */
  readonly size: number;
}

export interface DataTexturePages {
  /** a range of at least `count` elements */
  allocate(count: number): DataTextureRange;
  free(range: DataTextureRange): void;
  /** writes the element `ranges` of `data` (`[start, end)` pairs, all where undefined) into `range` */
  write(range: DataTextureRange, data: DataArray<any>, ranges?: readonly number[]): void;
  /** bytes of the pages */
  readonly pageBytes: number;
  readonly bytesPerElement: number;
  destroy(): void;
}

const DataPageTextures = new WeakSet<object>();

/** Whether `texture` is a page of data arrays: never a render target, so programs skip binding it again */
export function isDataPageTexture(texture: object) {
  return DataPageTextures.has(texture);
}

/** Bytes per element of the formats */
export const DataFormatBytes: { readonly [F in DataFormat]: number } = { rgba8: 4, r8: 1, rgba32f: 16 };

/** Elements per page: 1M of the 8-bit formats (4 MB of `rgba8`), 256K of `rgba32f` (4 MB) */
const PageElements: { readonly [F in DataFormat]: number } = { rgba8: 1 << 20, r8: 1 << 20, rgba32f: 1 << 18 };
/** Element offsets and sizes: 4 keeps small ranges from splitting rows into many pieces */
const Alignment = 4;

export function createDataTexturePages(
  format: DataFormat,
  createTexture: (format: DataFormat) => GLTexture,
  maxTextureSize: number,
): DataTexturePages {
  type Page = { readonly texture: GLTexture; readonly allocator: RangeAllocator; readonly dedicated: boolean };
  type Range = DataTextureRange & { readonly page: Page };
  const pages: Page[] = [];
  // a page as wide as allowed up to 2048 texels, its rows the elements per page
  const pageWidth = Math.min(2048, maxTextureSize);
  const pageHeight = Math.min(maxTextureSize, Math.ceil(PageElements[format] / pageWidth));
  const pageElements = pageWidth * pageHeight;
  let pageBytes = 0;

  const createPage = (width: number, height: number, dedicated: boolean): Page => {
    const texture = createTexture(format);
    DataPageTextures.add(texture);
    // allocated at its size on the first write (`loadData`)
    texture.define(width, height);
    const page = { texture, allocator: createRangeAllocator(width * height, Alignment), dedicated };
    pages.push(page);
    pageBytes += width * height * DataFormatBytes[format];
    return page;
  };

  const drop = (page: Page) => {
    pages.splice(pages.indexOf(page), 1);
    const { capacity } = page.allocator;
    pageBytes -= capacity * DataFormatBytes[format];
    page.texture.destroy();
  };

  return {
    allocate(count) {
      const n = Math.max(1, count);
      if (n > pageElements / 4) {
        // a texture of its own, sized as a data texture of the aligned count: the allocator rounds the range up to
        // `Alignment`, which a texture of exactly `n` texels (a square count) does not hold
        const aligned = Math.ceil(n / Alignment) * Alignment;
        const { width, height } = getDataTextureSize(aligned);
        const page = createPage(width, height, true);
        const offset = page.allocator.allocate(n);
        if (offset < 0) throw new Error(`data page of ${width}x${height} texels too small for ${n} elements`);
        return { texture: page.texture, width, height, offset, size: n, page } as Range;
      }
      for (const page of pages) {
        if (page.dedicated) continue;
        const offset = page.allocator.allocate(n);
        if (offset >= 0) {
          return { texture: page.texture, width: pageWidth, height: pageHeight, offset, size: n, page } as Range;
        }
      }
      const page = createPage(pageWidth, pageHeight, false);
      const offset = page.allocator.allocate(n);
      return { texture: page.texture, width: pageWidth, height: pageHeight, offset, size: n, page } as Range;
    },
    free(range) {
      const { page, offset, size } = range as Range;
      // in order with the draws before: GL commands run in the order issued, so the range is reusable at once
      page.allocator.free(offset, size);
      if (page.allocator.used !== 0) return;
      // the last shared page stays for the next allocations
      if (!page.dedicated && pages.filter((p) => !p.dedicated).length <= 1) return;
      drop(page);
    },
    write(range, data, ranges) {
      range.texture.loadData!(data, range.width, range.height, range.offset, ranges);
    },
    get pageBytes() {
      return pageBytes;
    },
    bytesPerElement: DataFormatBytes[format],
    destroy() {
      for (const p of pages) p.texture.destroy();
      pages.length = 0;
      pageBytes = 0;
    },
  };
}
