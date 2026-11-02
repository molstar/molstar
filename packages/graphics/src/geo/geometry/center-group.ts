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
import type { TransformData } from './transform-data.js';

/** Floats per item of a position group: center xyz and group */
export const PositionGroupStride = 4;

/**
 * Center (xyz) and group per item, shared by the geometries drawn as one quad of 6 vertices per item (spheres,
 * points): one array (`positionGroup`), which the shaders read the items from (`tPositionGroup`, vertex pulling), the
 * derived bounding sphere and group mapping. Its item order has no meaning beyond the geometry (groups name the items):
 * with levels of detail it is the levels' order (`orderPositionGroup`)
 */
export interface CenterGroupData {
  /** center xyz and group per item (`PositionGroupStride` floats), in level of detail order once strides are applied */
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

  // the items' set, not their order: a bounding sphere computed in another order differs slightly
  const hashCode = () => hashFnv32a([getCount(), getContentVersion(positionGroup)]);

  return {
    positionGroup,
    getBoundingSphere() {
      const newHash = hashCode();
      if (newHash !== currentHash) {
        Sphere3D.copy(boundingSphere, calculateItemOrderBoundingSphere(positionGroup, getCount()));
        currentHash = newHash;
      }
      return boundingSphere;
    },
    getGroupMapping() {
      if (positionGroup.ref.version !== currentGroup) {
        groupMapping = createGroupMapping(positionGroup.ref.value.array, getCount(), PositionGroupStride, 3);
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

/**
 * The strides `positionGroup`'s items are ordered by (none: in their order) and the levels' offsets, by its cell (the
 * geometry's and the values')
 */
const appliedStrides = new WeakMap<
  ValueCell<DataArray<Float32Array>>,
  { readonly strides: readonly number[]; readonly offsets: number[] }
>();

/**
 * The bounding sphere of the `count` items of `positionGroup` in their order (not the levels'): computed in another
 * order it differs slightly, and a geometry would get another sphere when built with levels than without
 */
function calculateItemOrderBoundingSphere(positionGroup: ValueCell<DataArray<Float32Array>>, count: number) {
  const { array } = positionGroup.ref.value;
  const applied = appliedStrides.get(positionGroup);
  if (!applied) return calculateInvariantBoundingSphere(array, count, 1, PositionGroupStride);
  // the centers in the items' order, while computing
  const centers = new Float32Array(count * 3);
  forEachInLevelOrder(count, applied.strides, (item, position) => {
    centers[item * 3] = array[position * PositionGroupStride];
    centers[item * 3 + 1] = array[position * PositionGroupStride + 1];
    centers[item * 3 + 2] = array[position * PositionGroupStride + 2];
  });
  return calculateInvariantBoundingSphere(centers, count, 1);
}

/** The version a reorder of `positionGroup` (`orderPositionGroup`) left, and the version of the content it reordered */
const reorderedVersions = new WeakMap<
  ValueCell<DataArray<Float32Array>>,
  { readonly order: number; readonly content: number }
>();

/** The version of `positionGroup`'s items, which a reorder keeps (it changes the cell's version) */
function getContentVersion(positionGroup: ValueCell<DataArray<Float32Array>>) {
  const reordered = reorderedVersions.get(positionGroup);
  return reordered && reordered.order === positionGroup.ref.version ? reordered.content : positionGroup.ref.version;
}

/**
 * Calls `f(item, position)` for the `count` items in the order of the levels of `strides`: per level (coarsest first)
 * the items its stride takes (`item % stride === 0`) that no level before took, then the items of no level (drawn by
 * none, kept to be ordered again); returns the vertex offsets (6 per item) of the levels, undefined without strides
 */
function forEachInLevelOrder(count: number, strides: readonly number[], f: (item: number, position: number) => void) {
  if (strides.length === 0) {
    for (let j = 0; j < count; ++j) f(j, j);
    return;
  }
  // each item's level, computed once (a byte each, while ordering)
  const levels = new Uint8Array(count);
  for (let j = 0; j < count; ++j) {
    let k = 0;
    while (k < strides.length && j % strides[k] !== 0) ++k;
    levels[j] = k;
  }
  const offsets = [0];
  let o = 0;
  for (let i = 0, il = strides.length; i <= il; ++i) {
    for (let j = 0; j < count; ++j) {
      if (levels[j] === i) f(j, o++);
    }
    if (i < il) offsets.push(o * 6);
  }
  return offsets;
}

function stridesEqual(a: readonly number[], b: readonly number[]) {
  if (a.length !== b.length) return false;
  for (let i = 0, il = a.length; i < il; ++i) if (a[i] !== b[i]) return false;
  return true;
}

/**
 * Orders the `count` items of `positionGroup` (`PositionGroupStride` floats each) into the levels of `strides`, in
 * place: from the order of the strides applied before, so that the order depends on `strides` only (not on the levels
 * set before); returns the vertex offsets (6 per item) of the levels, undefined without strides
 */
export function orderPositionGroup(
  positionGroup: ValueCell<DataArray<Float32Array>>,
  count: number,
  strides: readonly number[],
): number[] | undefined {
  const applied = appliedStrides.get(positionGroup);
  const previous = applied?.strides ?? [];
  if (stridesEqual(previous, strides)) return applied?.offsets.slice();

  const { array } = positionGroup.ref.value;
  const source = array.slice(0, count * PositionGroupStride);
  // the position of each item in the previous order
  let positionOf: Uint32Array | undefined;
  if (previous.length > 0) {
    const p = new Uint32Array(count);
    forEachInLevelOrder(count, previous, (item, position) => {
      p[item] = position;
    });
    positionOf = p;
  }
  const offsets = forEachInLevelOrder(count, strides, (item, position) => {
    const s = (positionOf ? positionOf[item] : item) * PositionGroupStride;
    const t = position * PositionGroupStride;
    array[t] = source[s];
    array[t + 1] = source[s + 1];
    array[t + 2] = source[s + 2];
    array[t + 3] = source[s + 3];
  });
  if (offsets) appliedStrides.set(positionGroup, { strides: strides.slice(), offsets: offsets.slice() });
  else appliedStrides.delete(positionGroup);
  const content = getContentVersion(positionGroup);
  ValueCell.update(positionGroup, positionGroup.ref.value);
  reorderedVersions.set(positionGroup, { order: positionGroup.ref.version, content });
  return offsets;
}

/** Sets the `count` items of `array` (`PositionGroupStride` floats each) to `positionGroup`, in their order */
export function setPositionGroup(
  positionGroup: ValueCell<DataArray<Float32Array>>,
  array: Float32Array,
  count: number,
) {
  appliedStrides.delete(positionGroup);
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
  const m = transform.aTransform.ref.value;
  const getLocation = (groupIndex: number, instanceIndex: number) => {
    if (instanceIndex < 0) {
      Vec3.fromArray(p, v, groupIndex * PositionGroupStride);
    } else {
      Vec3.transformMat4Offset(p, v, m, 0, groupIndex * PositionGroupStride, instanceIndex * 16);
    }
    return location;
  };
  return LocationIterator(count, instanceCount, 1, getLocation);
}

/**
 * The largest group of the `count` items of `positionGroup` plus one (0 without items); `stride` the floats per item,
 * the group the 4th (also of the cylinders' data, `CylinderStride`)
 */
export function getPositionGroupGroupCount(
  positionGroup: DataArray<Float32Array>,
  count: number,
  stride = PositionGroupStride,
) {
  const { array } = positionGroup;
  let max = -1;
  for (let i = 3, il = count * stride; i < il; i += stride) if (array[i] > max) max = array[i];
  return max + 1;
}

/** Interleaves `centers` (xyz) and `groups` of `count` items into `out` (reused where long enough) */
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
    array[i * 4 + 3] = groups[i];
  }
  return array;
}
