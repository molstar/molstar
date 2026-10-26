/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { DataFormat } from '../renderable/schema.js';
import { createDataTexturePages } from '../webgl/data-pages.js';

/** A texture that only records its size: allocating defines it, nothing is written */
function stubTexture() {
  return { define: () => {}, destroy: () => {} } as any;
}

/** Elements above which an array gets a texture of its own (a quarter of a page) */
const DedicatedAbove: { [F in DataFormat]: number } = { rgba8: 1 << 18, r8: 1 << 18, rgba32f: 1 << 16 };

describe('WebGL data pages', () => {
  // the page's allocator rounds a range up to a multiple of 4: a texture of exactly n texels (513², a square count)
  // did not hold it, and the range was used at offset -1 (every element read the one before it)
  it('give a dedicated range offset 0 for square counts', () => {
    const pages = createDataTexturePages('rgba8', stubTexture, 4096);
    const range = pages.allocate(513 * 513);
    expect(range.offset).toBe(0);
    expect(range.width * range.height).toBeGreaterThanOrEqual(513 * 513);
  });

  it('fit every dedicated range in its texture, around squares and exact sizes', () => {
    for (const format of ['rgba8', 'r8', 'rgba32f'] as const) {
      const pages = createDataTexturePages(format, stubTexture, 4096);
      const counts: number[] = [];
      for (let side = Math.ceil(Math.sqrt(DedicatedAbove[format])) + 1; side < 600; side += 7) {
        for (const d of [-2, -1, 0, 1, 2]) counts.push(side * side + d, side * (side + 1) + d);
      }
      for (const n of counts) {
        const range = pages.allocate(n);
        expect(range.offset, `${format} ${n}`).toBe(0);
        expect(range.width * range.height, `${format} ${n}`).toBeGreaterThanOrEqual(Math.ceil(n / 4) * 4);
        pages.free(range);
      }
      pages.destroy();
    }
  });
});
