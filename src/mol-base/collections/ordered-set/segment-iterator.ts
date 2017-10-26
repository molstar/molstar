/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../iterator'
import OrderedSet from '../ordered-set'

class SegmentIterator implements Iterator<{ segment: number } & OrderedSet.IndexRange> {
    private segmentRange = OrderedSet.IndexRange();
    private setRange = OrderedSet.IndexRange();
    private value = { segment: 0, start: 0, end: 0 };
    private last: number = 0;

    [Symbol.iterator]() { return new SegmentIterator(this.segments, this.set, this.start, this.end); };
    done: boolean = false;

    next() {
        const value = this.move();
        return { value, done: this.done }
    }

    move() {
        this.done = this.segmentRange.end <= this.segmentRange.start;
        while (!this.done) {
            if (!this.updateValue()) {
                this.updateSegmentRange();
            } else {
                this.value.segment = this.segmentRange.start++;
                break;
            }
        }
        return this.value;
    }

    private getSegmentIndex(value: number) {
        if (value >= this.last) return -1;
        return OrderedSet.findPredecessorIndex(this.segments, value - 1);
    }

    private updateValue() {
        const segmentEnd = OrderedSet.getAt(this.segments, this.segmentRange.start + 1);
        const setEnd = OrderedSet.findPredecessorIndexInRange(this.set, segmentEnd, this.setRange);
        this.value.start = this.setRange.start;
        this.value.end = setEnd;
        this.setRange.start = setEnd;
        return setEnd > this.value.start;
    }

    private updateSegmentRange() {
        const min = OrderedSet.getAt(this.set, this.setRange.start), max = OrderedSet.getAt(this.set, this.setRange.end - 1);
        this.segmentRange.start = this.getSegmentIndex(min);
        this.segmentRange.end = this.getSegmentIndex(max) + 1;
        this.done = this.segmentRange.end <= this.segmentRange.start;
    }

    constructor(private segments: OrderedSet, private set: OrderedSet, private start: number, private end: number) {
        this.last = OrderedSet.max(segments);
        this.setRange.start = start;
        this.setRange.end = end;
        this.updateSegmentRange();
    }
}

function createIterator(segments: OrderedSet, set: OrderedSet, range?: OrderedSet.IndexRange) {
    const start = !!range ? range.start : 0;
    const end = !!range ? range.end : OrderedSet.size(set);
    return new SegmentIterator(segments, set, start, end);
}

export default createIterator