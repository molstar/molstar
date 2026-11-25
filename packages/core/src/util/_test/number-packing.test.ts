/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { arrayMaxPackedIntToRGBA, packIntToRGBAArray, packIntToRGBArray, unpackRGBAToInt } from '../number-packing.js';

describe('number packing', () => {
  it('packs integers into four bytes past 2^24, zero for none', () => {
    const bytes = new Uint8Array(4);
    for (const v of [0, 1, 255, 65535, 16777214, 16777215, 16777216, 16777217, 16777219, 2 ** 31, 2 ** 32 - 2]) {
      packIntToRGBAArray(v, bytes, 0);
      expect(unpackRGBAToInt(bytes[0], bytes[1], bytes[2], bytes[3])).toBe(v);
    }
    // all zero: none
    expect(unpackRGBAToInt(0, 0, 0, 0)).toBe(-1);
  });

  it('packs below 2^24 as the three bytes did, the fourth zero', () => {
    const rgb = new Uint8Array(4);
    const rgba = new Uint8Array(4);
    for (const v of [0, 7, 70000, 16777214]) {
      packIntToRGBArray(v, rgb, 0);
      packIntToRGBAArray(v, rgba, 0);
      expect(Array.from(rgba)).toEqual([rgb[0], rgb[1], rgb[2], 0]);
    }
  });

  it('finds the largest packed integer', () => {
    const bytes = new Uint8Array(12);
    packIntToRGBAArray(5, bytes, 0);
    packIntToRGBAArray(16777219, bytes, 4);
    packIntToRGBAArray(9, bytes, 8);
    expect(arrayMaxPackedIntToRGBA(bytes, 4)).toBe(16777219);
  });
});
