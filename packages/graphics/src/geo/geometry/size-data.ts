/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { type DataArray, createDataArray } from '@molstar/graphics/gl/renderable/util';
import type { LocationIterator } from '../util/location-iterator.js';
import { type Location, NullLocation } from '@molstar/model/model/location';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import { Geometry } from './geometry.js';
import { unpackRGBToInt, packIntToRGBArray } from '@molstar/core/util/number-packing';

export type SizeType = 'uniform' | 'instance' | 'group' | 'groupInstance' | 'vertex' | 'vertexInstance';

export type SizeData = {
  uSize: ValueCell<number>;
  tSize: ValueCell<DataArray<Uint8Array>>;
  /** the index of `tSize` by the size type (`SizeTypeIndex`) */
  uSizeIndexType: ValueCell<number>;
  /** the size type's define (`getSizeTypeDefine`) */
  dSizeType: ValueCell<string>;
  /** the size type (`SizeType`) */
  sizeType: ValueCell<string>;
};

/**
 * How the size types read from `tSize` by index select it (`uSizeIndexType`): they share the define value 'texture',
 * one program, the shaders select at runtime
 */
export const SizeTypeIndex: { readonly [T in SizeType]?: number } = {
  instance: 0,
  group: 1,
  groupInstance: 2,
  vertex: 3,
  vertexInstance: 4,
};

/** The define value of the size type `type`: the types read from `tSize` by index share 'texture' */
export function getSizeTypeDefine(type: SizeType | string): string {
  return type === 'uniform' ? 'uniform' : 'texture';
}

/** Sets the size type `type`, its define and index (`SizeTypeIndex`) of `values` */
export function setSizeType(values: Pick<SizeData, 'sizeType' | 'dSizeType' | 'uSizeIndexType'>, type: SizeType) {
  ValueCell.updateIfChanged(values.sizeType, type);
  ValueCell.updateIfChanged(values.dSizeType, getSizeTypeDefine(type));
  ValueCell.updateIfChanged(values.uSizeIndexType, SizeTypeIndex[type] ?? 0);
}

/** The values of the size type `type` (its define and index, see `setSizeType`) */
function createSizeType(type: SizeType) {
  return {
    uSizeIndexType: ValueCell.create(SizeTypeIndex[type] ?? 0),
    dSizeType: ValueCell.create(getSizeTypeDefine(type)),
    sizeType: ValueCell.create<string>(type),
  };
}

/** Whether sizes of `type` are sampled at vertex positions (the position iterator of `createSizes`) */
export function isPositionSizeType(type: SizeType) {
  return type === 'vertex' || type === 'vertexInstance';
}

export function createSizes(
  locationIt: LocationIterator,
  positionIt: LocationIterator,
  sizeTheme: SizeTheme<any>,
  sizeData?: SizeData,
): SizeData {
  switch (Geometry.getGranularity(locationIt, sizeTheme.granularity)) {
    case 'uniform':
      return createUniformSize(locationIt, sizeTheme.size, sizeData);
    case 'instance':
      return createInstanceSize(locationIt, sizeTheme.size, sizeData);
    case 'group':
      return createGroupSize(locationIt, sizeTheme.size, sizeData);
    case 'groupInstance':
      return createGroupInstanceSize(locationIt, sizeTheme.size, sizeData);
    case 'vertex':
      return createVertexSize(positionIt, sizeTheme.size, sizeData);
    case 'vertexInstance':
      return createVertexInstanceSize(positionIt, sizeTheme.size, sizeData);
  }
}

export const sizeDataFactor = 100; // NOTE same factor is set in shaders

export function getMaxSize(sizeData: SizeData): number {
  const type = sizeData.sizeType.ref.value as SizeType;
  switch (type) {
    case 'uniform':
      return sizeData.uSize.ref.value;
    case 'instance':
    case 'group':
    case 'groupInstance':
    case 'vertex':
    case 'vertexInstance':
      let maxSize = 0;
      // the elements only: a reused array may be longer
      const { array, count } = sizeData.tSize.ref.value;
      for (let i = 0, il = count * 4; i < il; i += 4) {
        const value = unpackRGBToInt(array[i], array[i + 1], array[i + 2]);
        if (maxSize < value) maxSize = value;
      }
      return maxSize / sizeDataFactor;
  }
}

export type LocationSize = (location: Location) => number;

const emptySizeTexture = { array: new Uint8Array(4), count: 1 };
function createEmptySizeTexture() {
  return { tSize: ValueCell.create(emptySizeTexture) };
}

export function createValueSize(value: number, sizeData?: SizeData): SizeData {
  if (sizeData) {
    ValueCell.update(sizeData.uSize, value);
    setSizeType(sizeData, 'uniform');
    return sizeData;
  } else {
    return {
      uSize: ValueCell.create(value),
      ...createEmptySizeTexture(),
      ...createSizeType('uniform'),
    };
  }
}

/** Creates size uniform */
export function createUniformSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
  locationIt.reset();
  const location = locationIt.hasNext ? locationIt.move().location : NullLocation;
  return createValueSize(sizeFn(location), sizeData);
}

export function createTextureSize(sizes: DataArray<Uint8Array>, type: SizeType, sizeData?: SizeData): SizeData {
  if (sizeData) {
    ValueCell.update(sizeData.tSize, sizes);
    setSizeType(sizeData, type);
    return sizeData;
  } else {
    return {
      uSize: ValueCell.create(0),
      tSize: ValueCell.create(sizes),
      ...createSizeType(type),
    };
  }
}

/** Creates size texture with size for each instance/unit */
export function createInstanceSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
  const { instanceCount } = locationIt;
  const sizes = createDataArray(Math.max(1, instanceCount), 4, Uint8Array, sizeData && sizeData.tSize.ref.value.array);
  locationIt.reset();
  while (locationIt.hasNext && !locationIt.isNextNewInstance) {
    const v = locationIt.move();
    packIntToRGBArray(sizeFn(v.location) * sizeDataFactor, sizes.array, v.instanceIndex * 4);
    locationIt.skipInstance();
  }
  return createTextureSize(sizes, 'instance', sizeData);
}

/** Creates size texture with size for each group (i.e. shared across instances/units) */
export function createGroupSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
  const { groupCount } = locationIt;
  const sizes = createDataArray(Math.max(1, groupCount), 4, Uint8Array, sizeData && sizeData.tSize.ref.value.array);
  locationIt.reset();
  while (locationIt.hasNext && !locationIt.isNextNewInstance) {
    const v = locationIt.move();
    packIntToRGBArray(sizeFn(v.location) * sizeDataFactor, sizes.array, v.groupIndex * 4);
  }
  return createTextureSize(sizes, 'group', sizeData);
}

/** Creates size texture with size for each group in each instance (i.e. for each unit) */
export function createGroupInstanceSize(
  locationIt: LocationIterator,
  sizeFn: LocationSize,
  sizeData?: SizeData,
): SizeData {
  const { groupCount, instanceCount } = locationIt;
  const count = instanceCount * groupCount;
  const sizes = createDataArray(Math.max(1, count), 4, Uint8Array, sizeData && sizeData.tSize.ref.value.array);
  locationIt.reset();
  while (locationIt.hasNext) {
    const v = locationIt.move();
    packIntToRGBArray(sizeFn(v.location) * sizeDataFactor, sizes.array, v.index * 4);
  }
  return createTextureSize(sizes, 'groupInstance', sizeData);
}

/** Creates size texture with size for each vertex */
export function createVertexSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
  const { groupCount } = locationIt;
  const sizes = createDataArray(Math.max(1, groupCount), 4, Uint8Array, sizeData && sizeData.tSize.ref.value.array);
  locationIt.reset();
  while (locationIt.hasNext) {
    const v = locationIt.move();
    packIntToRGBArray(sizeFn(v.location) * sizeDataFactor, sizes.array, v.index * 4);
  }
  return createTextureSize(sizes, 'vertex', sizeData);
}

/** Creates size texture with size for each vertex instance */
export function createVertexInstanceSize(
  locationIt: LocationIterator,
  sizeFn: LocationSize,
  sizeData?: SizeData,
): SizeData {
  const { groupCount, instanceCount } = locationIt;
  const count = instanceCount * groupCount;
  const sizes = createDataArray(Math.max(1, count), 4, Uint8Array, sizeData && sizeData.tSize.ref.value.array);
  locationIt.reset();
  while (locationIt.hasNext) {
    const v = locationIt.move();
    packIntToRGBArray(sizeFn(v.location) * sizeDataFactor, sizes.array, v.index * 4);
  }
  return createTextureSize(sizes, 'vertexInstance', sizeData);
}
