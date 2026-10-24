/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sort, arraySwap } from './sort.js';
import type { AssignableArrayLike } from '@molstar/core/util/type-helpers';

type Bucket = {
  key: any;
  count: number;
  offset: number;
  /** position in the bucket list before any sort */
  ordinal: number;
};

function sortAsc(bs: Bucket[], i: number, j: number) {
  return bs[i].key < bs[j].key ? -1 : 1;
}

function _makeBuckets(
  indices: AssignableArrayLike<number>,
  getKey: (i: number) => any,
  sortBuckets: boolean,
  start: number,
  end: number,
) {
  const buckets = new Map<any, Bucket>();
  const bucketList: Bucket[] = [];
  // the bucket ordinal of each index, so that the reorder pass needs no getKey or map lookup; a run of one key reuses
  // its bucket without a lookup (perf-audit-2 01/P4)
  const ordinals = new Int32Array(end - start);

  let prev: Bucket | undefined;
  let isBucketed = true;
  for (let i = start; i < end; i++) {
    const key = getKey(indices[i]);
    let bucket = prev !== undefined && prev.key === key ? prev : buckets.get(key);
    if (bucket === undefined) {
      bucket = { key, count: 0, offset: i, ordinal: bucketList.length };
      buckets.set(key, bucket);
      bucketList[bucketList.length] = bucket;
    } else if (bucket !== prev) {
      isBucketed = false;
    }
    bucket.count++;
    ordinals[i - start] = bucket.ordinal;
    prev = bucket;
  }

  const bucketOffsets = new Int32Array(bucketList.length + 1);
  bucketOffsets[bucketList.length] = end;

  let sorted = true;
  if (sortBuckets) {
    for (let i = 1, _i = bucketList.length; i < _i; i++) {
      if (bucketList[i - 1].key > bucketList[i].key) {
        sorted = false;
        break;
      }
    }
  }

  if (isBucketed && sorted) {
    for (let i = 0; i < bucketList.length; i++) bucketOffsets[i] = bucketList[i].offset;
    return bucketOffsets;
  }

  const byOrdinal = sortBuckets && !sorted ? bucketList.slice() : bucketList;
  if (sortBuckets && !sorted) {
    sort(bucketList, 0, bucketList.length, sortAsc, arraySwap);
  }

  let offset = 0;
  for (let i = 0; i < bucketList.length; i++) {
    const b = bucketList[i];
    b.offset = offset;
    offset += b.count;
  }

  const reorderedIndices = new Int32Array(end - start);
  for (let i = start; i < end; i++) {
    const bucket = byOrdinal[ordinals[i - start]];
    reorderedIndices[bucket.offset++] = indices[i];
  }

  for (let i = 0, _i = reorderedIndices.length; i < _i; i++) {
    indices[i + start] = reorderedIndices[i];
  }

  bucketOffsets[0] = start;
  for (let i = 1; i < bucketList.length; i++) bucketOffsets[i] = bucketList[i - 1].offset + start;

  return bucketOffsets;
}

export interface MakeBucketsOptions<K> {
  // If specified, will be sorted
  sort?: boolean;
  // inclusive start indidex
  start?: number;
  // exclusive end index
  end?: number;
}

/**
 * Reorders indices so that the same keys are next to each other, [start, end)
 * Returns the offsets of buckets. So that [offsets[i], offsets[i + 1]) determines the range.
 */
export function makeBuckets<K extends string | number>(
  indices: AssignableArrayLike<number>,
  getKey: (i: number) => K,
  options?: MakeBucketsOptions<K>,
): ArrayLike<number> {
  const s = (options && options.start) || 0;
  const e = (options && options.end) || indices.length;
  if (e - s <= 0) throw new Error('Can only bucket non-empty collections.');

  return _makeBuckets(indices, getKey, !!(options && options.sort), s, e);
}
