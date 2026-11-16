/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Sphere3D } from '@molstar/core/math/geometry';
import { hashFnv32a } from '@molstar/core/data/util';
import {
  type DataArray,
  calculateInvariantBoundingSphere,
  createDataArray,
} from '@molstar/graphics/gl/renderable/util';
import { createGroupMapping, type GroupMapping } from '../util.js';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { LocationIterator, PositionLocation } from '../util/location-iterator.js';
import { type TransformData, getInstanceRows, getInstanceTransforms } from './transform-data.js';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';

/** Words per item of a position group: center xyz (floats) and group (an int32 as its bits, `int32ToFloatBits`) */
export const PositionGroupStride = 4;

const bitsFloat = new Float32Array(1);
const bitsInt = new Int32Array(bitsFloat.buffer);

/**
 * The float whose bits are the int32 `value`: a group in a data array of groups (`rgb32f-i32`: position groups, the
 * cylinders', lines' and text's data), exact past 2^24 (a float's value only to 2^24); written as is, never computed
 * with (small values are denormals). Negative groups (-1, none) and those from 2^31 - 2^23 are quiet NaNs, whose bits
 * typed-array copies and copies through numbers keep (`center-group.test.ts`)
 */
export function int32ToFloatBits(value: number) {
  bitsInt[0] = value;
  return bitsFloat[0];
}

/** The int32 of the bits of float `value` (`int32ToFloatBits`'s inverse) */
export function floatBitsToInt32(value: number) {
  bitsFloat[0] = value;
  return bitsInt[0];
}

/**
 * The words of a data array of groups (`rgb32f-i32`: position groups, the cylinders', lines' and text's data) as
 * int32s, the groups among them (`int32ToFloatBits`)
 */
export function getDataWords(array: Float32Array) {
  return new Int32Array(array.buffer, array.byteOffset, array.length);
}

/**
 * Center (xyz) and group per item, shared by the geometries drawn as one quad of 6 vertices per item (spheres,
 * points): one array (`positionGroup`, the group an int32 as its bits: `int32ToFloatBits`, `getDataWords`),
 * which the shaders read the items from (`tPositionGroup`, `rgb32f-i32`, vertex pulling), the
 * derived bounding sphere and group mapping. The items keep their order with levels of detail: the shaders map a level's
 * position to its item (`getLodItem`)
 */
export interface CenterGroupData {
  /** center xyz and group per item (`PositionGroupStride` floats), in item order */
  readonly positionGroup: ValueCell<DataArray<Float32Array>>;

  getBoundingSphere(): Sphere3D;
  /** group to item indices, times `PositionGroupStride` (offsets into `positionGroup`'s array) */
  getGroupMapping(): GroupMapping;
  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;
}

export function createCenterGroupData(array: Float32Array, count: number, getCount: () => number): CenterGroupData {
  const positionGroup = ValueCell.create(createDataArray(count, PositionGroupStride, Float32Array, array));

  const boundingSphere = Sphere3D();
  let groupMapping: GroupMapping;
  let currentHash = -1;
  let currentGroup = -1;

  const hashCode = () => hashFnv32a([getCount(), positionGroup.ref.version]);

  return {
    positionGroup,
    getBoundingSphere() {
      const newHash = hashCode();
      if (newHash !== currentHash) {
        Sphere3D.copy(
          boundingSphere,
          calculateInvariantBoundingSphere(positionGroup.ref.value.array, getCount(), 1, PositionGroupStride),
        );
        currentHash = newHash;
      }
      return boundingSphere;
    },
    getGroupMapping() {
      if (positionGroup.ref.version !== currentGroup) {
        const words = getDataWords(positionGroup.ref.value.array);
        groupMapping = createGroupMapping(words, getCount(), PositionGroupStride, 3);
        currentGroup = positionGroup.ref.version;
      }
      return groupMapping;
    },
    setBoundingSphere(sphere: Sphere3D) {
      Sphere3D.copy(boundingSphere, sphere);
      currentHash = hashCode();
    },
    hasBoundingSphere() {
      return currentHash === hashCode();
    },
  };
}

/** Sets the `count` items of `array` (`PositionGroupStride` floats each) to `positionGroup`, in their order */
export function setPositionGroup(
  positionGroup: ValueCell<DataArray<Float32Array>>,
  array: Float32Array,
  count: number,
) {
  ValueCell.update(positionGroup, createDataArray(count, PositionGroupStride, Float32Array, array));
}

/**
 * The positions of the `count` items of `positionGroup` (in its order, as the shaders draw them; the item as the
 * iterator's group), per instance of `transform`
 */
export function createPositionGroupIterator(
  positionGroup: ValueCell<DataArray<Float32Array>>,
  count: number,
  transform: TransformData,
): LocationIterator {
  const instanceCount = transform.instanceCount.ref.value;
  const location = PositionLocation();
  const p = location.position;
  const v = positionGroup.ref.value.array;
  const t = getInstanceTransforms(transform);
  // the instance's row: the grid stores them in its order (`getInstanceRows`)
  const rows = getInstanceRows(transform);
  const local = Vec3();
  const getLocation = (groupIndex: number, instanceIndex: number) => {
    if (instanceIndex < 0) {
      Vec3.fromArray(p, v, groupIndex * PositionGroupStride);
    } else {
      const row = rows ? rows[instanceIndex] : instanceIndex;
      InstanceTransforms.getPoint(t, row, Vec3.fromArray(local, v, groupIndex * PositionGroupStride), p);
    }
    return location;
  };
  return LocationIterator(count, instanceCount, 1, getLocation);
}

/**
 * The largest group of the `count` items of `positionGroup` plus one (0 without items); `stride` the floats per item,
 * the group the 4th, an int32 as its bits (also of the cylinders', lines' and text's data, `CylinderStride`)
 */
export function getPositionGroupGroupCount(
  positionGroup: DataArray<Float32Array>,
  count: number,
  stride = PositionGroupStride,
) {
  const array = getDataWords(positionGroup.array);
  let max = -1;
  for (let i = 3, il = count * stride; i < il; i += stride) if (array[i] > max) max = array[i];
  return max + 1;
}

/**
 * Interleaves `centers` (xyz) and `groups` of `count` items into `out` (reused where long enough), the groups as their
 * bits (`int32ToFloatBits`)
 */
export function interleavePositionGroup(
  centers: ArrayLike<number>,
  groups: ArrayLike<number>,
  count: number,
  out?: Float32Array,
): Float32Array {
  const length = count * PositionGroupStride;
  const array = out && out.length >= length ? out : new Float32Array(length);
  for (let i = 0; i < count; ++i) {
    array[i * 4] = centers[i * 3];
    array[i * 4 + 1] = centers[i * 3 + 1];
    array[i * 4 + 2] = centers[i * 3 + 2];
    array[i * 4 + 3] = int32ToFloatBits(groups[i]);
  }
  return array;
}
