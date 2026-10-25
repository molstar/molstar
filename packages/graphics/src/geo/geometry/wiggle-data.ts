/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { type DataArray, createDataArray } from '@molstar/graphics/gl/renderable/util';

export type WiggleType = 'instance' | 'groupInstance';

export type WiggleData = {
  tWiggle: ValueCell<DataArray<Uint8Array>>;
  dWiggle: ValueCell<boolean>;
  wiggleAverage: ValueCell<number>;
  /** the type (`WiggleType`); the object uniform `uInstanceGranularity` selects the index */
  wiggleType: ValueCell<string>;
  uWiggleStrength: ValueCell<number>;
};

export function applyWiggleValue(array: Uint8Array, start: number, end: number, value: number) {
  for (let i = start; i < end; ++i) {
    array[i] = value * 255;
  }
  return true;
}

export function getWiggleAverage(array: Uint8Array, count: number): number {
  if (count === 0 || array.length < count) return 0;
  let sum = 0;
  for (let i = 0; i < count; ++i) {
    sum += array[i];
  }
  return sum / (255 * count);
}

export function clearWiggle(array: Uint8Array, start: number, end: number) {
  array.fill(0, start, end);
}

export function createWiggle(count: number, type: WiggleType, wiggleData?: WiggleData): WiggleData {
  const wiggle = createDataArray(Math.max(1, count), 1, Uint8Array, wiggleData && wiggleData.tWiggle.ref.value.array);
  if (wiggleData) {
    ValueCell.update(wiggleData.tWiggle, wiggle);
    ValueCell.updateIfChanged(wiggleData.dWiggle, count > 0);
    ValueCell.updateIfChanged(wiggleData.wiggleAverage, getWiggleAverage(wiggle.array, count));
    ValueCell.updateIfChanged(wiggleData.wiggleType, type);
    return wiggleData;
  } else {
    return {
      tWiggle: ValueCell.create(wiggle),
      dWiggle: ValueCell.create(count > 0),
      wiggleAverage: ValueCell.create(0),
      wiggleType: ValueCell.create<string>(type),
      uWiggleStrength: ValueCell.create(1),
    };
  }
}

const emptyWiggleTexture = { array: new Uint8Array(4), count: 1 };
export function createEmptyWiggle(wiggleData?: WiggleData): WiggleData {
  if (wiggleData) {
    ValueCell.update(wiggleData.tWiggle, emptyWiggleTexture);
    return wiggleData;
  } else {
    return {
      tWiggle: ValueCell.create(emptyWiggleTexture),
      dWiggle: ValueCell.create(false),
      wiggleAverage: ValueCell.create(0),
      wiggleType: ValueCell.create<string>('groupInstance'),
      uWiggleStrength: ValueCell.create(1),
    };
  }
}
