/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { type DataArray, createDataArray } from '@molstar/graphics/gl/renderable/util';
import { Color } from '@molstar/core/util/color';
import { createNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import { getLayerTypeDefine, type LayerType } from './layer-type.js';

export type OverpaintType = LayerType;

export type OverpaintData = {
  tOverpaint: ValueCell<DataArray<Uint8Array>>;
  dOverpaint: ValueCell<boolean>;

  tOverpaintGrid: ValueCell<Texture>;
  uOverpaintGridDim: ValueCell<Vec3>;
  uOverpaintGridTexDim: ValueCell<Vec2>;
  uOverpaintGridTransform: ValueCell<Vec4>;
  /** the layer type (`OverpaintType`) */
  overpaintType: ValueCell<string>;
  /** its define (`getLayerTypeDefine`) */
  dOverpaintType: ValueCell<string>;
  uOverpaintStrength: ValueCell<number>;
};

export function applyOverpaintColor(array: Uint8Array, start: number, end: number, color: Color) {
  for (let i = start; i < end; ++i) {
    Color.toArray(color, array, i * 4);
    array[i * 4 + 3] = 255;
  }
  return true;
}

export function clearOverpaint(array: Uint8Array, start: number, end: number) {
  array.fill(0, start * 4, end * 4);
  return true;
}

/** Sets the layer type `type` and its define (`getLayerTypeDefine`) of `values` */
export function setOverpaintType(
  values: Pick<OverpaintData, 'overpaintType' | 'dOverpaintType'>,
  type: OverpaintType | string,
) {
  ValueCell.updateIfChanged(values.overpaintType, type);
  ValueCell.updateIfChanged(values.dOverpaintType, getLayerTypeDefine(type));
}

export function createOverpaint(count: number, type: OverpaintType, overpaintData?: OverpaintData): OverpaintData {
  const overpaint = createDataArray(
    Math.max(1, count),
    4,
    Uint8Array,
    overpaintData && overpaintData.tOverpaint.ref.value.array,
  );
  if (overpaintData) {
    ValueCell.update(overpaintData.tOverpaint, overpaint);
    ValueCell.updateIfChanged(overpaintData.dOverpaint, count > 0);
    setOverpaintType(overpaintData, type);
    return overpaintData;
  } else {
    return {
      tOverpaint: ValueCell.create(overpaint),
      dOverpaint: ValueCell.create(count > 0),

      tOverpaintGrid: ValueCell.create(createNullTexture()),
      uOverpaintGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uOverpaintGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uOverpaintGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      overpaintType: ValueCell.create<string>(type),
      dOverpaintType: ValueCell.create(getLayerTypeDefine(type)),
      uOverpaintStrength: ValueCell.create(1),
    };
  }
}

const emptyOverpaintTexture = { array: new Uint8Array(4), count: 1 };
export function createEmptyOverpaint(overpaintData?: OverpaintData): OverpaintData {
  if (overpaintData) {
    ValueCell.update(overpaintData.tOverpaint, emptyOverpaintTexture);
    return overpaintData;
  } else {
    return {
      tOverpaint: ValueCell.create(emptyOverpaintTexture),
      dOverpaint: ValueCell.create(false),

      tOverpaintGrid: ValueCell.create(createNullTexture()),
      uOverpaintGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uOverpaintGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uOverpaintGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      overpaintType: ValueCell.create<string>('groupInstance'),
      dOverpaintType: ValueCell.create(getLayerTypeDefine('groupInstance')),
      uOverpaintStrength: ValueCell.create(1),
    };
  }
}
