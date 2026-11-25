/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import {
  createSiteTypeSize,
  createTextureSize,
  createValueSize,
  getMaxSize,
  getSizeTypeDefine,
  type SizeType,
  SizeTypeIndex,
} from '../size-data.js';

describe('size-data', () => {
  it('shares one define for the types read by index, which the index selects', () => {
    const textureTypes: SizeType[] = ['instance', 'group', 'groupInstance', 'vertex', 'vertexInstance'];
    for (const type of textureTypes) expect(getSizeTypeDefine(type)).toBe('texture');
    expect(new Set(textureTypes.map((t) => SizeTypeIndex[t])).size).toBe(textureTypes.length);
    expect(getSizeTypeDefine('uniform')).toBe('uniform');
  });

  it('keeps the size type, its define and index in step', () => {
    const data = createTextureSize(createDataArray(4, 4, Uint8Array), 'group');
    expect(data.sizeType.ref.value).toBe('group');
    expect(data.dSizeType.ref.value).toBe('texture');
    expect(data.uSizeIndexType.ref.value).toBe(SizeTypeIndex.group);

    createTextureSize(createDataArray(4, 4, Uint8Array), 'vertexInstance', data);
    expect(data.sizeType.ref.value).toBe('vertexInstance');
    expect(data.uSizeIndexType.ref.value).toBe(SizeTypeIndex.vertexInstance);

    createValueSize(2, data);
    expect(data.sizeType.ref.value).toBe('uniform');
    expect(data.dSizeType.ref.value).toBe('uniform');
    expect(getMaxSize(data)).toBe(2);
  });

  it("makes a uniform size's table of one entry its own array, not the empty one every uniform size shares", () => {
    const a = createValueSize(1);
    const b = createValueSize(1);
    expect(a.tSize.ref.value).toBe(b.tSize.ref.value);
    createSiteTypeSize([3], a);
    expect(a.tSize.ref.value.array).not.toBe(b.tSize.ref.value.array);
    expect(Array.from(b.tSize.ref.value.array)).toEqual([0, 0, 0, 0]);
  });
});
