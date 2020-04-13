/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../../iterator';
import OrderedSet from '../ordered-set';
import Interval from '../interval';
import SortedArray from '../sorted-array';
import Segs from '../segmentation';

interface Segmentation {
    /** Segments stored as a sorted array */
    offsets: SortedArray,
    /** Mapping of values to segments */
    index: Int32Array,
    /** Number of segments */
    count: number
}

export function create(values: ArrayLike<number>): Segmentation {
    const offsets = SortedArray.ofSortedArray(values);
    const max = SortedArray.max(offsets);
    const index = new Int32Array(max);
    for (let i = 0, _i = values.length - 1; i < _i; i++) {
        for (let j = values[i], _j = values[i + 1]; j < _j; j++) {
            index[j] = i;
        }
    }
    return { offsets, index, count: values.length - 1 };
}

export function ofOffsets(offsets: ArrayLike<number>, bounds: Interval): Segmentation {
    const s = Interval.start(bounds);
    const segments = new Int32Array(offsets.length + 1);
    for (let i = 0, _i = offsets.length; i < _i; i++) {
        segments[i] = offsets[i] - s;
    }
    segments[offsets.length] = Interval.end(bounds) - s;
    return create(segments);
}

/** Get number of segments in a segmentation */
export function count({ count }: Segmentation) { return count; }
export function getSegment({ index }: Segmentation, value: number) { return index[value]; }

export function projectValue({ offsets }: Segmentation, set: OrderedSet, value: number): Interval {
    const last = OrderedSet.max(offsets);
    const idx = value >= last ? -1 : OrderedSet.findPredecessorIndex(offsets, value - 1);
    return OrderedSet.findRange(set, OrderedSet.getAt(offsets, idx), OrderedSet.getAt(offsets, idx + 1) - 1);
}

export class SegmentIterator<I extends number = number> implements Iterator<Segs.Segment<I>> {
    private segmentMin = 0;
    private segmentMax = 0;
    private setRange = Interval.Empty;
    private value: Segs.Segment<I> = { index: 0 as I, start: 0, end: 0 };

    hasNext: boolean = false;

    move() {
        while (this.hasNext) {
            if (this.updateValue()) {
                this.value.index = this.segmentMin++ as I;
                this.hasNext = this.segmentMax >= this.segmentMin && Interval.size(this.setRange) > 0;
                break;
            } else {
                this.updateSegmentRange();
            }
        }
        return this.value;
    }

    private updateValue() {
        const segmentEnd = this.segments[this.segmentMin + 1];
        // TODO: add optimized version for interval and array?
        const setEnd = OrderedSet.findPredecessorIndexInInterval(this.set, segmentEnd, this.setRange);
        this.value.start = Interval.start(this.setRange);
        this.value.end = setEnd;
        this.setRange = Interval.ofBounds(setEnd, Interval.end(this.setRange));
        return setEnd > this.value.start;
    }

    private updateSegmentRange() {
        const sMin = Interval.min(this.setRange), sMax = Interval.max(this.setRange);
        if (sMax < sMin) {
            this.hasNext = false;
            return;
        }
        this.segmentMin = this.segmentMap[OrderedSet.getAt(this.set, sMin)];
        this.segmentMax = this.segmentMap[OrderedSet.getAt(this.set, sMax)];
        this.hasNext = this.segmentMax >= this.segmentMin;
    }

    setSegment(segment: Segs.Segment<number>) {
        this.setRange = Interval.ofBounds(segment.start, segment.end);
        this.updateSegmentRange();
    }

    constructor(private segments: SortedArray, private segmentMap: Int32Array, private set: OrderedSet, inputRange: Interval) {
        this.setRange = inputRange;
        this.updateSegmentRange();
    }
}

export function segments(segs: Segmentation, set: OrderedSet, segment?: Segs.Segment) {
    const int = typeof segment !== 'undefined' ? Interval.ofBounds(segment.start, segment.end) : Interval.ofBounds(0, OrderedSet.size(set));
    return new SegmentIterator(segs.offsets, segs.index, set, int);
}