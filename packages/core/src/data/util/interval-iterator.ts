/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Interval } from '../int/interval.js';
import { OrderedSet } from '../int/ordered-set.js';
import type { Segmentation } from '../int/segmentation.js';
import type { Iterator } from '../iterator.js';

/** Emits a segment of length one for each element in the interval that is also in the set */
export class IntervalIterator<I extends number = number> implements Iterator<Segmentation.Segment<I>> {
  private value: Segmentation.Segment<I> = { index: 0 as I, start: 0, end: 0 };

  private curIndex = 0;
  private maxIndex = 0;

  hasNext: boolean = false;

  updateValue() {
    this.value.index = this.curIndex as I;
    this.value.start = OrderedSet.findPredecessorIndex(this.set, Interval.getAt(this.interval, this.curIndex));
    this.value.end = this.value.start + 1;
  }

  /** Moves `curIndex` to the first interval index at or after it whose value is in the set */
  private skipToMember() {
    while (this.curIndex <= this.maxIndex && !OrderedSet.has(this.set, Interval.getAt(this.interval, this.curIndex))) {
      ++this.curIndex;
    }
    this.hasNext = this.curIndex <= this.maxIndex;
  }

  move() {
    if (this.hasNext) {
      this.updateValue();
      ++this.curIndex;
      this.skipToMember();
    }
    return this.value;
  }

  constructor(
    private interval: Interval<I>,
    private set: OrderedSet<I>,
  ) {
    if (Interval.size(interval) && OrderedSet.size(set)) {
      // interval indices of the first value >= min(set) and the last value <= max(set)
      this.curIndex = Interval.findPredecessorIndex(interval, OrderedSet.min(set));
      this.maxIndex = Interval.findPredecessorIndex(interval, OrderedSet.max(set) + 1) - 1;
      this.skipToMember();
    }
  }
}
