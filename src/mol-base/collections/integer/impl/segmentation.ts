/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../../iterator'
import OrderedSet from '../ordered-set'
import Interval from '../interval'
import SortedArray from '../sorted-array'
import Segs from '../segmentation'

type Segmentation = { segments: SortedArray, segmentMap: Int32Array, count: number }

export function create(values: ArrayLike<number>): Segmentation {
    const segments = SortedArray.ofSortedArray(values);
    const max = SortedArray.max(segments);
    const segmentMap = new Int32Array(max);
    for (let i = 0, _i = values.length - 1; i < _i; i++) {
        for (let j = values[i], _j = values[i + 1]; j < _j; j++) {
            segmentMap[j] = i;
        }
    }
    return { segments, segmentMap, count: values.length - 1 };
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

export function count({ count }: Segmentation) { return count; }
export function getSegment({ segmentMap }: Segmentation, value: number) { return segmentMap[value]; }

export function projectValue({ segments }: Segmentation, set: OrderedSet, value: number): Interval {
    const last = OrderedSet.max(segments);
    const idx = value >= last ? -1 : OrderedSet.findPredecessorIndex(segments, value - 1);
    return OrderedSet.findRange(set, OrderedSet.getAt(segments, idx), OrderedSet.getAt(segments, idx + 1) - 1);
}

class SegmentIterator implements Iterator<Segs.Segment> {
    private segmentStart = 0;
    private segmentEnd = 0;
    private setRange = Interval.Empty;
    private value: Segs.Segment = { index: 0, start: 0, end: 0 };
    private last: number = 0;

    hasNext: boolean = false;

    move() {
        while (this.hasNext) {
            if (this.updateValue()) {
                this.value.index = this.segmentStart++;
                this.hasNext = this.segmentEnd >= this.segmentStart && Interval.size(this.setRange) > 0;
                break;
            } else {
                this.updateSegmentRange();
            }
        }
        return this.value;
    }

    private getSegmentIndex(value: number) {
        if (value >= this.last) return -1;
        return OrderedSet.findPredecessorIndex(this.segments, value - 1);
    }

    private updateValue() {
        const segmentEnd = OrderedSet.getAt(this.segments, this.segmentStart + 1);
        const setEnd = OrderedSet.findPredecessorIndexInRange(this.set, segmentEnd, this.setRange);
        this.value.start = Interval.start(this.setRange);
        this.value.end = setEnd;
        const rEnd = Interval.end(this.setRange);
        this.setRange = Interval.ofBounds(setEnd, rEnd);
        return setEnd > this.value.start;
    }

    private updateSegmentRange() {
        const sMin = Interval.min(this.setRange), sMax = Interval.max(this.setRange);
        if (sMax < sMin) {
            this.hasNext = false;
            return;
        }

        const min = OrderedSet.getAt(this.set, Interval.min(this.setRange));
        const max = OrderedSet.getAt(this.set, Interval.max(this.setRange));

        this.segmentStart = this.getSegmentIndex(min);
        this.segmentEnd = this.getSegmentIndex(max) + 1;
        this.hasNext = this.segmentEnd >= this.segmentStart && Interval.size(this.setRange) > 0;
    }

    constructor(private segments: OrderedSet, private set: OrderedSet, inputRange: Interval) {
        this.last = OrderedSet.max(segments);
        this.setRange = inputRange;
        this.updateSegmentRange();
    }
}

export function segments(segs: Segmentation, set: OrderedSet, segment?: Segs.Segment) {
    const int = typeof segment !== 'undefined' ? Interval.ofBounds(segment.start, segment.end) : Interval.ofBounds(0, OrderedSet.size(set));
    return new SegmentIterator(segs.segments, set, int);
}