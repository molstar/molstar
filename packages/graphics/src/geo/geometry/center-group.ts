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

/**
 * Center (xyz) and group per item, shared by the geometries drawn as one quad of 6 vertices per item
 * (spheres, points): the buffers, the derived bounding sphere and group mapping, and the
 * `tPositionGroup` data the shaders read the item data from (vertex pulling)
 */
export interface CenterGroupData {
  readonly centerBuffer: ValueCell<Float32Array>;
  readonly groupBuffer: ValueCell<Float32Array>;
  /** center xyz and group per item, ordered by level of detail when strides are given */
  readonly positionGroup: ValueCell<DataArray<Float32Array>>;

  getBoundingSphere(): Sphere3D;
  getGroupMapping(): GroupMapping;
  setBoundingSphere(boundingSphere: Sphere3D): void;
  hasBoundingSphere(): boolean;
  /**
   * Writes the items to `positionGroup`; with `strides`, items are reordered into levels of detail
   * (items of a coarser level first) and the vertex offsets of the levels are returned
   */
  updatePositionGroup(strides?: number[]): number[] | undefined;
}

export function createCenterGroupData(
  centers: Float32Array,
  groups: Float32Array,
  getCount: () => number,
): CenterGroupData {
  const centerBuffer = ValueCell.create(centers);
  const groupBuffer = ValueCell.create(groups);
  const positionGroup = ValueCell.create(createDataArray(1, 4, Float32Array));

  const boundingSphere = Sphere3D();
  let groupMapping: GroupMapping;
  let currentHash = -1;
  let currentGroup = -1;

  const hashCode = () => hashFnv32a([getCount(), centerBuffer.ref.version, groupBuffer.ref.version]);

  return {
    centerBuffer,
    groupBuffer,
    positionGroup,
    getBoundingSphere() {
      const newHash = hashCode();
      if (newHash !== currentHash) {
        const b = calculateInvariantBoundingSphere(centerBuffer.ref.value, getCount(), 1);
        Sphere3D.copy(boundingSphere, b);
        currentHash = newHash;
      }
      return boundingSphere;
    },
    getGroupMapping() {
      if (groupBuffer.ref.version !== currentGroup) {
        groupMapping = createGroupMapping(groupBuffer.ref.value, getCount());
        currentGroup = groupBuffer.ref.version;
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
    updatePositionGroup(strides: number[] = []) {
      const count = getCount();
      const pgt = createDataArray(count, 4, Float32Array, positionGroup.ref.value.array);
      const offsets = setPositionGroup(pgt, centerBuffer.ref.value, groupBuffer.ref.value, count, strides);
      ValueCell.update(positionGroup, pgt);
      return offsets;
    },
  };
}

/**
 * Writes center and group per item to `out`; returns the vertex offsets (6 per item) of the levels
 * of detail when `strides` are given
 */
export function setPositionGroup(
  out: DataArray<Float32Array>,
  centers: Float32Array,
  groups: Float32Array,
  count: number,
  strides: number[],
) {
  const { array } = out;
  if (strides.length === 0) {
    for (let i = 0; i < count; ++i) {
      array[i * 4 + 0] = centers[i * 3 + 0];
      array[i * 4 + 1] = centers[i * 3 + 1];
      array[i * 4 + 2] = centers[i * 3 + 2];
      array[i * 4 + 3] = groups[i];
    }
    return;
  }

  const offsets = [0];

  let o = 0;
  for (let i = 0, il = strides.length; i < il; ++i) {
    const s = strides[i];
    for (let j = 0; j < count; ++j) {
      let handled = false;
      for (let k = 0; k < i; ++k) {
        if (j % strides[k] === 0) {
          handled = true;
          break;
        }
      }
      if (!handled && j % s === 0) {
        array[o * 4 + 0] = centers[j * 3 + 0];
        array[o * 4 + 1] = centers[j * 3 + 1];
        array[o * 4 + 2] = centers[j * 3 + 2];
        array[o * 4 + 3] = groups[j];
        o += 1;
      }
    }
    offsets.push(o * 6);
  }

  return offsets;
}
