/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import {
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
});
