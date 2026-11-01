/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util/value-cell';
import {
  ClipStorageVec4s,
  getInstanceClipCount,
  getPixelClipCount,
  getStorageClipCount,
  isStorageClipSupported,
  packStorageClipObjects,
} from '../clip-storage.js';

function clipValues(variant: string, clipping = false) {
  const frame = new Float32Array(32);
  const shape = new Float32Array(8);
  for (let i = 0; i < 32; ++i) frame[i] = i;
  for (let i = 0; i < 8; ++i) shape[i] = 100 + i;
  return {
    dClipObjects: ValueCell.create(true),
    dClipVariant: ValueCell.create(variant),
    dClipping: ValueCell.create(clipping),
    uClipObjectCount: ValueCell.create(2),
    uClipObjectFrame: ValueCell.create(frame),
    uClipObjectShape: ValueCell.create(shape),
  };
}

describe('clip objects in storage', () => {
  it('are those of an object, counted for its variant; instance clipping not with a clipping layer', () => {
    expect(getStorageClipCount(clipValues('instance'))).toBe(2);
    expect(getStorageClipCount(clipValues('pixel'))).toBe(2);
    expect(getInstanceClipCount(clipValues('instance'))).toBe(2);
    expect(getInstanceClipCount(clipValues('pixel'))).toBe(0);
    expect(getPixelClipCount(clipValues('pixel'))).toBe(2);
    expect(getPixelClipCount(clipValues('instance'))).toBe(0);
    // a clipping layer: its flags per sphere for pixel clipping, none per instance
    expect(getStorageClipCount(clipValues('instance', true))).toBe(2);
    expect(getInstanceClipCount(clipValues('instance', true))).toBe(0);
    expect(getPixelClipCount(clipValues('pixel', true))).toBe(2);
    expect(isStorageClipSupported(clipValues('instance', true))).toBe(false);
    expect(isStorageClipSupported(clipValues('pixel', true))).toBe(true);
    const none = clipValues('instance');
    ValueCell.update(none.dClipObjects, false);
    expect(getStorageClipCount(none)).toBe(0);
  });

  it('pack per object its frame columns, then its shape', () => {
    const out = new Float32Array(2 * ClipStorageVec4s * 4 + 4).fill(-1);
    expect(packStorageClipObjects(clipValues('instance'), out, 1)).toBe(2);
    expect(Array.from(out.subarray(0, 4))).toEqual([-1, -1, -1, -1]);
    expect(Array.from(out.subarray(4, 24))).toEqual([...Array.from({ length: 16 }, (_, i) => i), 100, 101, 102, 103]);
    expect(Array.from(out.subarray(24, 44))).toEqual([
      ...Array.from({ length: 16 }, (_, i) => 16 + i),
      104,
      105,
      106,
      107,
    ]);
  });
});
