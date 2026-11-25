/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import {
  floatBitsToInt32,
  getDataWords,
  getPositionGroupGroupCount,
  int32ToFloatBits,
  interleavePositionGroup,
} from '../center-group.js';

const Groups = [0, 1, 7, 16777215, 16777216, 16777217, 16777219, 2 ** 30 + 3, -1];

describe('position groups', () => {
  it('keep groups past 2^24 exact as their bits, through copies of the array', () => {
    const array = new Float32Array(Groups.length);
    for (let i = 0; i < Groups.length; i++) array[i] = int32ToFloatBits(Groups[i]);
    // copied as typed arrays are (merging, levels of detail) and element by element through numbers
    const copy = new Float32Array(array.length);
    for (let i = 0; i < array.length; i++) copy[i] = array[i];
    for (const a of [array, array.slice(), copy]) {
      expect(Array.from(getDataWords(a))).toEqual(Groups);
      expect(Array.from(a, floatBitsToInt32)).toEqual(Groups);
    }
  });

  it('interleave centers and groups, the groups as their bits, and count the groups', () => {
    const n = 3;
    const centers = Float32Array.of(1, 2, 3, 4, 5, 6, 7, 8, 9);
    const groups = [5, 16777219, 2];
    const array = interleavePositionGroup(centers, groups, n);
    const words = getDataWords(array);
    for (let i = 0; i < n; i++) {
      expect(Array.from(array.subarray(i * 4, i * 4 + 3))).toEqual(Array.from(centers.subarray(i * 3, i * 3 + 3)));
      expect(words[i * 4 + 3]).toBe(groups[i]);
    }
    expect(getPositionGroupGroupCount(createDataArray(n, 4, Float32Array, array), n)).toBe(16777220);
  });
});
