/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { type DataArray, createDataArray } from '@molstar/graphics/gl/renderable/util';
import { createNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import { Material } from '@molstar/core/util/material';
import { getLayerTypeDefine, type LayerType } from './layer-type.js';

export type SubstanceType = LayerType;

export type SubstanceData = {
  tSubstance: ValueCell<DataArray<Uint8Array>>;
  dSubstance: ValueCell<boolean>;

  tSubstanceGrid: ValueCell<Texture>;
  uSubstanceGridDim: ValueCell<Vec3>;
  uSubstanceGridTexDim: ValueCell<Vec2>;
  uSubstanceGridTransform: ValueCell<Vec4>;
  /** the layer type (`SubstanceType`) */
  substanceType: ValueCell<string>;
  /** its define (`getLayerTypeDefine`) */
  dSubstanceType: ValueCell<string>;
  uSubstanceStrength: ValueCell<number>;
};

export function applySubstanceMaterial(array: Uint8Array, start: number, end: number, material: Material) {
  for (let i = start; i < end; ++i) {
    Material.toArray(material, array, i * 4);
    array[i * 4 + 3] = 255;
  }
  return true;
}

export function clearSubstance(array: Uint8Array, start: number, end: number) {
  array.fill(0, start * 4, end * 4);
  return true;
}

/** Sets the layer type `type` and its define (`getLayerTypeDefine`) of `values` */
export function setSubstanceType(
  values: Pick<SubstanceData, 'substanceType' | 'dSubstanceType'>,
  type: SubstanceType | string,
) {
  ValueCell.updateIfChanged(values.substanceType, type);
  ValueCell.updateIfChanged(values.dSubstanceType, getLayerTypeDefine(type));
}

export function createSubstance(count: number, type: SubstanceType, substanceData?: SubstanceData): SubstanceData {
  const substance = createDataArray(
    Math.max(1, count),
    4,
    Uint8Array,
    substanceData && substanceData.tSubstance.ref.value.array,
  );
  if (substanceData) {
    ValueCell.update(substanceData.tSubstance, substance);
    ValueCell.updateIfChanged(substanceData.dSubstance, count > 0);
    setSubstanceType(substanceData, type);
    return substanceData;
  } else {
    return {
      tSubstance: ValueCell.create(substance),
      dSubstance: ValueCell.create(count > 0),

      tSubstanceGrid: ValueCell.create(createNullTexture()),
      uSubstanceGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uSubstanceGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uSubstanceGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      substanceType: ValueCell.create<string>(type),
      dSubstanceType: ValueCell.create(getLayerTypeDefine(type)),
      uSubstanceStrength: ValueCell.create(1),
    };
  }
}

const emptySubstanceTexture = { array: new Uint8Array(4), count: 1 };
export function createEmptySubstance(substanceData?: SubstanceData): SubstanceData {
  if (substanceData) {
    ValueCell.update(substanceData.tSubstance, emptySubstanceTexture);
    return substanceData;
  } else {
    return {
      tSubstance: ValueCell.create(emptySubstanceTexture),
      dSubstance: ValueCell.create(false),

      tSubstanceGrid: ValueCell.create(createNullTexture()),
      uSubstanceGridDim: ValueCell.create(Vec3.create(1, 1, 1)),
      uSubstanceGridTexDim: ValueCell.create(Vec2.create(1, 1)),
      uSubstanceGridTransform: ValueCell.create(Vec4.create(0, 0, 0, 1)),
      substanceType: ValueCell.create<string>('groupInstance'),
      dSubstanceType: ValueCell.create(getLayerTypeDefine('groupInstance')),
      uSubstanceStrength: ValueCell.create(1),
    };
  }
}
