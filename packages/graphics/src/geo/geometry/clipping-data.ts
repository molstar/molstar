/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { type DataArray, createDataArray } from '@molstar/graphics/gl/renderable/util';
import type { Clipping } from '@molstar/graphics/theme/clipping';

export type ClippingType = 'instance' | 'groupInstance';

export type ClippingData = {
  tClipping: ValueCell<DataArray<Uint8Array>>;
  dClipping: ValueCell<boolean>;
  /** the type (`ClippingType`); the object uniform `uInstanceGranularity` selects the index */
  clippingType: ValueCell<string>;
};

export function applyClippingGroups(array: Uint8Array, start: number, end: number, groups: Clipping.Groups) {
  array.fill(groups, start, end);
  return true;
}

export function clearClipping(array: Uint8Array, start: number, end: number) {
  array.fill(0, start, end);
}

export function createClipping(count: number, type: ClippingType, clippingData?: ClippingData): ClippingData {
  const clipping = createDataArray(
    Math.max(1, count),
    1,
    Uint8Array,
    clippingData && clippingData.tClipping.ref.value.array,
  );
  if (clippingData) {
    ValueCell.update(clippingData.tClipping, clipping);
    ValueCell.updateIfChanged(clippingData.dClipping, count > 0);
    ValueCell.updateIfChanged(clippingData.clippingType, type);
    return clippingData;
  } else {
    return {
      tClipping: ValueCell.create(clipping),
      dClipping: ValueCell.create(count > 0),
      clippingType: ValueCell.create<string>(type),
    };
  }
}

const emptyClippingTexture = { array: new Uint8Array(4), count: 1 };
export function createEmptyClipping(clippingData?: ClippingData): ClippingData {
  if (clippingData) {
    ValueCell.update(clippingData.tClipping, emptyClippingTexture);
    ValueCell.updateIfChanged(clippingData.dClipping, false);
    return clippingData;
  } else {
    return {
      tClipping: ValueCell.create(emptyClippingTexture),
      dClipping: ValueCell.create(false),
      clippingType: ValueCell.create<string>('groupInstance'),
    };
  }
}
