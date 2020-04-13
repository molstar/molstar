/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Iterator from '../iterator';
import { OrderedSet, Interval, Segmentation } from '../int';

/** Emits a segment of length one for each element in the interval that is also in the set */
export class IntervalIterator<I extends number = number> implements Iterator<Segmentation.Segment<I>> {
    private value: Segmentation.Segment<I> = { index: 0 as I, start: 0, end: 0 }

    private curIndex = 0
    private maxIndex = 0

    hasNext: boolean = false;

    updateValue() {
        this.value.index = this.curIndex as I;
        this.value.start = OrderedSet.findPredecessorIndex(this.set, Interval.getAt(this.interval, this.curIndex));
        this.value.end = this.value.start + 1;
    }

    move() {
        if (this.hasNext) {
            this.updateValue();
            while (this.curIndex <= this.maxIndex) {
                ++this.curIndex;
                if (OrderedSet.has(this.set, this.curIndex)) break;
            }
            this.hasNext = this.curIndex <= this.maxIndex;
        }
        return this.value;
    }

    constructor(private interval: Interval<I>, private set: OrderedSet<I>) {
        if (Interval.size(interval)) {
            this.curIndex = Interval.findPredecessorIndex(interval, OrderedSet.min(set));
            this.maxIndex = Interval.findPredecessorIndex(interval, OrderedSet.max(set));
        }

        this.hasNext = OrderedSet.areIntersecting(this.interval, this.set);
    }
}