/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

type Bucket = {
    count: number,
    offset: number
}

function _makeBuckets(indices: Helpers.ArrayLike<number>, getKey: (i: number) => any, start: number, end: number) {
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
            const bucket: Bucket = { count: 1, offset: i };
            buckets.set(key, bucket);
            bucketList[bucketList.length] = bucket;
        }
        prevKey = key;
    }

    const bucketOffsets = new Int32Array(bucketList.length + 1);
    bucketOffsets[bucketList.length] = end;

    if (isBucketed) {
        for (let i = 0; i < bucketList.length; i++) bucketOffsets[i] = bucketList[i].offset;
        return bucketOffsets;
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

/**
 * Reorders indices so that the same keys are next to each other, [start, end)
 * Returns the offsets of buckets. So that [offsets[i], offsets[i + 1]) determines the range.
 */
export function makeBuckets<T>(indices: Helpers.ArrayLike<number>, getKey: (i: number) => string | number, start?: number, end?: number): ArrayLike<number> {
    const s = start || 0;
    const e = typeof end === 'undefined' ? indices.length : end;

    if (e - s <= 0) throw new Error('Can only bucket non-empty collections.');

    return _makeBuckets(indices, getKey, s, e);
}