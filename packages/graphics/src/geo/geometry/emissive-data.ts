/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { type DataArray, createDataArray } from '@molstar/graphics/gl/renderable/util';
import { createNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import { getLayerTypeDefine, type LayerType } from './layer-type.js';

export type EmissiveType = LayerType;

export type EmissiveData = {
  tEmissive: ValueCell<DataArray<Uint8Array>>;
  dEmissive: ValueCell<boolean>;
  emissiveAverage: ValueCell<number>;

  tEmissiveGrid: ValueCell<Texture>;
  uEmissiveGridDim: ValueCell<Vec3>;
  uEmissiveGridTexDim: ValueCell<Vec2>;
  uEmissiveGridTransform: ValueCell<Vec4>;
  /** the layer type (`EmissiveType`) */
  emissiveType: ValueCell<string>;
  /** its define (`getLayerTypeDefine`) */
  dEmissiveType: ValueCell<string>;
  uEmissiveStrength: ValueCell<number>;
};

export function applyEmissiveValue(array: Uint8Array, start: number, end: number, value: number) {
  for (let i = start; i < end; ++i) {
    array[i] = value * 255;
  }
  return true;
}

export function getEmissiveAverage(array: Uint8Array, count: number): number {
  if (count === 0 || array.length < count) return 0;
  let sum = 0;
  for (let i = 0; i < count; ++i) {
    sum += array[i];
  }
  return sum / (255 * count);
}

export function clearEmissive(array: Uint8Array, start: number, end: number) {
  array.fill(0, start, end);
}

/** Sets the layer type `type` and its define (`getLayerTypeDefine`) of `values` */
export function setEmissiveType(
  values: Pick<EmissiveData, 'emissiveType' | 'dEmissiveType'>,
  type: EmissiveType | string,
) {
  ValueCell.updateIfChanged(values.emissiveType, type);
  ValueCell.updateIfChanged(values.dEmissiveType, getLayerTypeDefine(type));
}

export function createEmissive(count: number, type: EmissiveType, emissiveData?: EmissiveData): EmissiveData {
  const emissive = createDataArray(
    Math.max(1, count),
    1,
    Uint8Array,
    emissiveData && emissiveData.tEmissive.ref.value.array,
  );
  if (emissiveData) {
    ValueCell.update(emissiveData.tEmissive, emissive);
    ValueCell.updateIfChanged(emissiveData.dEmissive, count > 0);
    ValueCell.updateIfChanged(emissiveData.emissiveAverage, getEmissiveAverage(emissive.array, count));
    setEmissiveType(emissiveData, type);
    return emissiveData;
  } else {
    return {
      tEmissive: ValueCell.create(emissive),
      dEmissive: ValueCell.create(count > 0),
      emissiveAverage: ValueCell.create(0),

      tEmissiveGrid: ValueCell.create(createNullTexture()),
      uEmissiveGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uEmissiveGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uEmissiveGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      emissiveType: ValueCell.create<string>(type),
      dEmissiveType: ValueCell.create(getLayerTypeDefine(type)),
      uEmissiveStrength: ValueCell.create(1),
    };
  }
}

const emptyEmissiveTexture = { array: new Uint8Array(4), count: 1 };
export function createEmptyEmissive(emissiveData?: EmissiveData): EmissiveData {
  if (emissiveData) {
    ValueCell.update(emissiveData.tEmissive, emptyEmissiveTexture);
    return emissiveData;
  } else {
    return {
      tEmissive: ValueCell.create(emptyEmissiveTexture),
      dEmissive: ValueCell.create(false),
      emissiveAverage: ValueCell.create(0),

      tEmissiveGrid: ValueCell.create(createNullTexture()),
      uEmissiveGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uEmissiveGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uEmissiveGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      emissiveType: ValueCell.create<string>('groupInstance'),
      dEmissiveType: ValueCell.create(getLayerTypeDefine('groupInstance')),
      uEmissiveStrength: ValueCell.create(1),
    };
  }
}
