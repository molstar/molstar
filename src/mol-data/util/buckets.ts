/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sort, arraySwap } from './sort';
import { AssignableArrayLike } from '../../mol-util/type-helpers';

type Bucket = {
    key: any,
    count: number,
    offset: number
}

function sortAsc(bs: Bucket[], i: number, j: number) { return bs[i].key < bs[j].key ? -1 : 1; }

function _makeBuckets(indices: AssignableArrayLike<number>,
    getKey: (i: number) => any, sortBuckets: boolean, start: number, end: number) {

    const buckets = new Map<any, Bucket>();
    const bucketList: Bucket[] = [];

    let prevKey = getKey(indices[0]);
    let isBucketed = true;
    for (let i = start; i < end; i++) {
        const key = getKey(indices[i]);
        if (buckets.has(key)) {
            buckets.get(key)!.count++;
            if (prevKey !== key) isBucketed = false;
        } else {
            const bucket: Bucket = { key, count: 1, offset: i };
            buckets.set(key, bucket);
            bucketList[bucketList.length] = bucket;
        }
        prevKey = key;
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
        const key = getKey(indices[i]);
        const bucket = buckets.get(key)!;
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
    sort?: boolean,
    // inclusive start indidex
    start?: number,
    // exclusive end index
    end?: number
}

/**
 * Reorders indices so that the same keys are next to each other, [start, end)
 * Returns the offsets of buckets. So that [offsets[i], offsets[i + 1]) determines the range.
 */
export function makeBuckets<K extends string | number>(
    indices: AssignableArrayLike<number>, getKey: (i: number) => K, options?: MakeBucketsOptions<K>): ArrayLike<number> {
    const s = (options && options.start) || 0;
    const e = (options && options.end) || indices.length;
    if (e - s <= 0) throw new Error('Can only bucket non-empty collections.');

    return _makeBuckets(indices, getKey, !!(options && options.sort), s, e);
}