/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import { Color } from '@molstar/core/util/color';
import {
  type ColorType,
  ColorTypeIndex,
  createTextureColor,
  createValueColor,
  getColorTypeDefine,
  setColorType,
} from '../color-data.js';

describe('color-data', () => {
  it('shares one define for the types read by index and one for the grid types, which the index selects', () => {
    const textureTypes: ColorType[] = ['instance', 'group', 'groupInstance', 'vertex', 'vertexInstance'];
    for (const type of textureTypes) expect(getColorTypeDefine(type)).toBe('texture');
    expect(new Set(textureTypes.map((t) => ColorTypeIndex[t])).size).toBe(textureTypes.length);
    const gridTypes: ColorType[] = ['volume', 'volumeInstance'];
    for (const type of gridTypes) expect(getColorTypeDefine(type)).toBe('grid');
    expect(new Set(gridTypes.map((t) => ColorTypeIndex[t])).size).toBe(gridTypes.length);
    for (const type of ['uniform', 'direct'] as ColorType[]) expect(getColorTypeDefine(type)).toBe(type);
  });

  it('keeps the color type, its define and index in step', () => {
    const data = createTextureColor(createDataArray(4, 4, Uint8Array), 'groupInstance');
    expect(data.colorType.ref.value).toBe('groupInstance');
    expect(data.dColorType.ref.value).toBe('texture');
    expect(data.uColorIndexType.ref.value).toBe(ColorTypeIndex.groupInstance);

    setColorType(data, 'volumeInstance');
    expect(data.colorType.ref.value).toBe('volumeInstance');
    expect(data.dColorType.ref.value).toBe('grid');
    expect(data.uColorIndexType.ref.value).toBe(ColorTypeIndex.volumeInstance);

    createValueColor(Color(0xff0000), data);
    expect(data.colorType.ref.value).toBe('uniform');
    expect(data.dColorType.ref.value).toBe('uniform');
  });
});
