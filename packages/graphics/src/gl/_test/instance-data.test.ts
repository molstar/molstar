/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import type { RenderableValues } from '../renderable/schema.js';
import { getInstanceDataSegmentCount, InstanceDataVec4s, packInstanceData } from '../renderable/merged.js';
import { writeAffineRowsOf, writeAffineRowsRange } from '../transform-parts.js';

/** 10 instances translated by (i, 2i, 3i), ids 100 + i; segment 0 (groups 3, 10, 20) the first 6, segment 1 (5, 30, 40) the rest */
function createValues() {
  const count = 10;
  const transforms = new Float32Array(count * 16);
  const instances = new Uint32Array(count);
  const segments = new Int32Array(count * 4);
  for (let i = 0; i < count; ++i) {
    Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(i, 2 * i, 3 * i)), transforms, i * 16);
    instances[i] = 100 + i;
    segments.set(i < 6 ? [3, 10, 20, 0] : [5, 30, 40, 1], i * 4);
  }
  return {
    aTransform: ValueCell.create(transforms),
    instanceCount: ValueCell.create(count),
    aInstance: ValueCell.create(instances),
    aSegment: ValueCell.create(segments),
  } as unknown as RenderableValues;
}

describe('merged instance data (WebGL2 pages)', () => {
  it('writes the rows of a range of instances as the rows of all from its first', () => {
    const values = createValues() as any;
    const all = writeAffineRowsOf(values, 10, new Float32Array(10 * 12));
    const range = writeAffineRowsRange(values, 4, 3, new Float32Array(3 * 12));
    expect(Array.from(range)).toEqual(Array.from(all.subarray(4 * 12, 7 * 12)));
  });

  it('packs a page: its instances, then its segments’ groups after a fixed slot', () => {
    const values = createValues();
    expect(getInstanceDataSegmentCount(values)).toBe(2);
    const stride = InstanceDataVec4s * 4;
    // instances 4 to 7 (segments 0 and 1), the segments from slot 4
    const out = new Float32Array((4 + 2) * stride);
    packInstanceData(values, out, 4, 4, 4);
    const bits = new Uint32Array(out.buffer);
    const signed = new Int32Array(out.buffer);
    for (let k = 0; k < 4; ++k) {
      const i = 4 + k;
      // the translation in the rows' w
      expect([out[k * stride + 3], out[k * stride + 7], out[k * stride + 11]]).toEqual([i, 2 * i, 3 * i]);
      expect(bits[k * stride + 12]).toBe(100 + i);
      expect(signed[k * stride + 13]).toBe(i < 6 ? 0 : 1);
    }
    expect(Array.from(signed.subarray(4 * stride, 4 * stride + 4))).toEqual([3, 10, 20, 0]);
    expect(Array.from(signed.subarray(5 * stride, 5 * stride + 4))).toEqual([5, 30, 40, 0]);
  });
});
