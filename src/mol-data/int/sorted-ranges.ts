/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Segmentation, OrderedSet, SortedArray, Interval } from '../int';
import _Iterator from '../iterator';

/** Pairs of min and max indices of sorted, non-overlapping ranges */
type SortedRanges<T extends number = number> = SortedArray<T>

namespace SortedRanges {
    export function ofSortedRanges<T extends number = number>(array: ArrayLike<T>) { return SortedArray.ofSortedArray<T>(array) }
    export function start<T extends number = number>(ranges: SortedRanges<T>) { return ranges[0] }
    export function end<T extends number = number>(ranges: SortedRanges<T>) { return ranges[ranges.length - 1] + 1 }
    export function min<T extends number = number>(ranges: SortedRanges<T>) { return ranges[0] }
    export function max<T extends number = number>(ranges: SortedRanges<T>) { return ranges[ranges.length - 1] }
    export function size<T extends number = number>(ranges: SortedRanges<T>) {
        let size = 0
        for (let i = 0, il = ranges.length; i < il; i += 2) {
            size += ranges[i + 1] - ranges[i] + 1
        }
        return size
    }

    export function transientSegments<T extends number = number, I extends number = number>(ranges: SortedRanges<T>, set: OrderedSet<T>) {
        return new Iterator<T, I>(ranges, set)
    }

    export class Iterator<T extends number = number, I extends number = number> implements _Iterator<Segmentation.Segment<I>> {
        private value: Segmentation.Segment<I> = { index: 0 as I, start: 0 as T, end: 0 as T }

        private curIndex = 0
        private maxIndex = 0
        private interval: Interval<T>
        private curMin: T = 0 as T

        hasNext: boolean = false;

        private updateInterval() {
            this.interval = Interval.ofRange(this.ranges[this.curIndex], this.ranges[this.curIndex + 1])
        }

        private updateValue() {
            this.value.index = this.curIndex / 2 as I
            this.value.start = OrderedSet.findPredecessorIndex(this.set, this.ranges[this.curIndex])
            this.value.end = OrderedSet.findPredecessorIndex(this.set, this.ranges[this.curIndex + 1])
        }

        move() {
            if (this.hasNext) {
                this.updateValue()
                while (this.curIndex <= this.maxIndex) {
                    this.curIndex += 2
                    this.curMin = Interval.end(this.interval)
                    this.updateInterval()
                    if (Interval.min(this.interval) >= this.curMin && OrderedSet.areIntersecting(this.interval, this.set)) break
                }
                this.hasNext = this.curIndex <= this.maxIndex
            }
            return this.value;
        }

        private getRangeIndex(value: number) {
            const index = SortedArray.findPredecessorIndex(this.ranges, value)
            return (index % 2 === 1) ? index - 1 : index
        }

        constructor(private ranges: SortedRanges<T>, private set: OrderedSet<T>) {
            const min = OrderedSet.min(set)
            const max = OrderedSet.max(set)
            const b = ranges.length > 0 && ranges[0] <= max && ranges[ranges.length - 1] >= min
            if (b) {
                this.curIndex = this.getRangeIndex(min)
                this.maxIndex = Math.min(ranges.length - 2, this.getRangeIndex(max) + 2)
                this.curMin = ranges[this.curIndex]
                this.updateInterval()
            }
            this.hasNext = b && this.curIndex <= this.maxIndex
        }
    }
}

export default SortedRanges