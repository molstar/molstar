/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

 import Iterator from '../iterator'
import { OrderedSet, Interval, Segmentation } from '../int';

/** Emits a segment of length one for each element in the interval that is also in the set */
export class IntervalIterator<T extends number = number> implements Iterator<Segmentation.Segment<T>> {
    private value: Segmentation.Segment<T> = { index: 0, start: 0 as T, end: 0 as T }

    private curIndex = 0
    private maxIndex = 0

    hasNext: boolean = false;

    updateValue() {
        this.value.index = this.curIndex
        this.value.start = OrderedSet.findPredecessorIndex(this.set, Interval.getAt(this.interval, this.curIndex)) as  T
        this.value.end = this.value.start + 1 as T 
    }

    move() {
        if (this.hasNext) {
            this.updateValue()
            while (this.curIndex <= this.maxIndex) {
                ++this.curIndex
                if (OrderedSet.has(this.set, this.curIndex)) break
            }
            this.hasNext = this.curIndex <= this.maxIndex
        }
        return this.value;
    }

    constructor(private interval: Interval<T>, private set: OrderedSet<T>) {
        if (Interval.size(interval)) {
            this.curIndex = Interval.findPredecessorIndex(interval, OrderedSet.min(set))
            this.maxIndex = Interval.findPredecessorIndex(interval, OrderedSet.max(set))
        }

        this.hasNext = OrderedSet.areIntersecting(this.interval, this.set)
    }
}