/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import { type DataArray, createDataArray } from '@molstar/graphics/gl/renderable/util';

export type MarkerType = 'instance' | 'groupInstance';

export type MarkerData = {
  uMarker: ValueCell<number>;
  tMarker: ValueCell<DataArray<Uint8Array>>;
  markerAverage: ValueCell<number>;
  markerStatus: ValueCell<number>;
  /** the type (`MarkerType`); the object uniform `uInstanceGranularity` selects the index */
  markerType: ValueCell<string>;
  /**
   * per instance, how many of its entries are marked (0 or 1 with 'instance' markers), kept by `Visual.mark` and sized
   * by it (shorter: nothing marked yet); written in place, its cell updated where replaced
   */
  markerInstanceCounts: ValueCell<Uint32Array>;
  /** the instances with a marked entry (`markerInstanceCounts` above 0) */
  markedInstanceCount: ValueCell<number>;
};

const MarkerCountLut = new Uint8Array(0x0303 + 1);
MarkerCountLut[0x0001] = 1;
MarkerCountLut[0x0002] = 1;
MarkerCountLut[0x0003] = 1;
MarkerCountLut[0x0100] = 1;
MarkerCountLut[0x0200] = 1;
MarkerCountLut[0x0300] = 1;
MarkerCountLut[0x0101] = 2;
MarkerCountLut[0x0201] = 2;
MarkerCountLut[0x0301] = 2;
MarkerCountLut[0x0102] = 2;
MarkerCountLut[0x0202] = 2;
MarkerCountLut[0x0302] = 2;
MarkerCountLut[0x0103] = 2;
MarkerCountLut[0x0203] = 2;
MarkerCountLut[0x0303] = 2;

/**
 * Calculates the average number of entries that have any marker flag set.
 *
 * For alternative implementations and performance tests see
 * `src\perf-tests\markers-average.ts`.
 */
export function getMarkersAverage(array: Uint8Array, count: number): number {
  if (count === 0) return 0;

  const view = new Uint32Array(array.buffer, 0, array.buffer.byteLength >> 2);
  const viewEnd = (count - 4) >> 2;
  const backStart = 4 * viewEnd;

  let sum = 0;
  if (viewEnd < 0) {
    // avoid edge cases with small arrays
    for (let i = 0; i < count; ++i) {
      sum += array[i] && 1;
    }
  } else {
    for (let i = 0; i < viewEnd; ++i) {
      const v = view[i];
      sum += MarkerCountLut[v & 0xffff] + MarkerCountLut[v >> 16];
    }
    for (let i = backStart; i < count; ++i) {
      sum += array[i] && 1;
    }
  }
  return sum / count;
}

/**
 * Adds `sign` times the marked entries of `[start, end)` of `array` to the counts of their instances (`groupCount`
 * entries an instance), `sign` -1 before an action changes the range and 1 after; returns the change of the number of
 * instances with a marked entry. Reads the range alone: a mark costs what it marks
 */
export function addMarkerInstanceCounts(
  array: Uint8Array,
  start: number,
  end: number,
  groupCount: number,
  counts: Uint32Array,
  sign: 1 | -1,
): number {
  let changed = 0;
  let i = start;
  while (i < end) {
    const instance = Math.floor(i / groupCount);
    const e = Math.min(end, (instance + 1) * groupCount);
    let n = 0;
    for (let k = i; k < e; ++k) if (array[k] !== 0) n += 1;
    if (n > 0) {
      const before = counts[instance];
      const after = before + sign * n;
      counts[instance] = after;
      if (before === 0) changed += 1;
      else if (after === 0) changed -= 1;
    }
    i = e;
  }
  return changed;
}

/** The counts of all `instanceCount` instances from `array`; returns the instances with a marked entry */
export function setMarkerInstanceCounts(
  array: Uint8Array,
  instanceCount: number,
  groupCount: number,
  counts: Uint32Array,
): number {
  counts.fill(0, 0, instanceCount);
  return addMarkerInstanceCounts(array, 0, instanceCount * groupCount, groupCount, counts, 1);
}

/**
 * The counts of the instances counted above 0 read again from `array`, after an action that only removes markers (no
 * other instance gains one); returns the instances with a marked entry
 */
export function updateMarkerInstanceCounts(
  array: Uint8Array,
  instanceCount: number,
  groupCount: number,
  counts: Uint32Array,
): number {
  let marked = 0;
  for (let i = 0; i < instanceCount; ++i) {
    if (counts[i] === 0) continue;
    let n = 0;
    for (let k = i * groupCount, kl = k + groupCount; k < kl; ++k) if (array[k] !== 0) n += 1;
    counts[i] = n;
    if (n > 0) marked += 1;
  }
  return marked;
}

const EmptyCounts = new Uint32Array(0);

function resetInstanceCounts(markerData: MarkerData) {
  if (!markerData.markerInstanceCounts) {
    markerData.markerInstanceCounts = ValueCell.create(EmptyCounts);
    markerData.markedInstanceCount = ValueCell.create(0);
    return;
  }
  markerData.markerInstanceCounts.ref.value.fill(0);
  ValueCell.updateIfChanged(markerData.markedInstanceCount, 0);
}

export function createMarkers(count: number, type: MarkerType, markerData?: MarkerData): MarkerData {
  if (markerData) {
    const markers = createDataArray(Math.max(1, count), 1, Uint8Array, markerData.tMarker.ref.value.array);
    markers.array.fill(0, 0, count);
    ValueCell.updateIfChanged(markerData.uMarker, 0);
    ValueCell.update(markerData.tMarker, markers);
    ValueCell.updateIfChanged(markerData.markerAverage, 0);
    ValueCell.updateIfChanged(markerData.markerStatus, 0);
    ValueCell.updateIfChanged(markerData.markerType, type);
    resetInstanceCounts(markerData);
    return markerData;
  } else {
    const markers = createDataArray(Math.max(1, count), 1, Uint8Array);
    return {
      uMarker: ValueCell.create(0),
      tMarker: ValueCell.create(markers),
      markerAverage: ValueCell.create(0),
      markerStatus: ValueCell.create(0),
      markerType: ValueCell.create<string>(type),
      markerInstanceCounts: ValueCell.create(EmptyCounts),
      markedInstanceCount: ValueCell.create(0),
    };
  }
}

const emptyMarkerTexture = { array: new Uint8Array(4), count: 1 };
export function createEmptyMarkers(markerData?: MarkerData): MarkerData {
  if (markerData) {
    ValueCell.updateIfChanged(markerData.uMarker, 0);
    ValueCell.update(markerData.tMarker, emptyMarkerTexture);
    ValueCell.updateIfChanged(markerData.markerAverage, 0);
    ValueCell.updateIfChanged(markerData.markerStatus, 0);
    resetInstanceCounts(markerData);
    return markerData;
  } else {
    return {
      uMarker: ValueCell.create(0),
      tMarker: ValueCell.create(emptyMarkerTexture),
      markerAverage: ValueCell.create(0),
      markerStatus: ValueCell.create(0),
      markerType: ValueCell.create<string>('groupInstance'),
      markerInstanceCounts: ValueCell.create(EmptyCounts),
      markedInstanceCount: ValueCell.create(0),
    };
  }
}
