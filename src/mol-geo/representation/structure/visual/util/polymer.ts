/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Element, StructureProperties } from 'mol-model/structure';
import { Segmentation } from 'mol-data/int';
// import Iterator from 'mol-data/iterator';

export function getPolymerElementCount(unit: Unit) {
    let count = 0
    const { elements } = unit
    const l = Element.Location(unit)
    if (Unit.isAtomic(unit)) {
        const { polymerSegments, residueSegments } = unit.model.atomicHierarchy
        const polymerIt = Segmentation.transientSegments(polymerSegments, elements);
        const residuesIt = Segmentation.transientSegments(residueSegments, elements);
        while (polymerIt.hasNext) {
            residuesIt.setSegment(polymerIt.move());
            while (residuesIt.hasNext) {
                residuesIt.move();
                count++
            }
        }
    } else if (Unit.isSpheres(unit)) {
        for (let i = 0, il = elements.length; i < il; ++i) {
            l.element = elements[i]
            if (StructureProperties.entity.type(l) === 'polymer') count++
        }
    }
    return count
}

// interface PolymerTrace {
//     center: Element
//     direction: Element
// }

// export class PolymerTraceIterator<T extends number = number> implements Iterator<PolymerTrace> {
//     private segmentMin = 0;
//     private segmentMax = 0;
//     private setRange = Interval.Empty;
//     private value: TraceSegment = { index: 0, start: 0 as T, end: 0 as T };

//     hasNext: boolean = false;

//     move() {
//         while (this.hasNext) {
//             if (this.updateValue()) {
//                 this.value.index = this.segmentMin++;
//                 this.hasNext = this.segmentMax >= this.segmentMin && Interval.size(this.setRange) > 0;
//                 break;
//             } else {
//                 this.updateSegmentRange();
//             }
//         }
//         return this.value;
//     }

//     private updateValue() {
//         const segmentEnd = this.segments[this.segmentMin + 1];
//         // TODO: add optimized version for interval and array?
//         const setEnd = OrderedSet.findPredecessorIndexInInterval(this.set, segmentEnd, this.setRange);
//         this.value.start = Interval.start(this.setRange) as T;
//         this.value.end = setEnd as T;
//         this.setRange = Interval.ofBounds(setEnd, Interval.end(this.setRange));
//         return setEnd > this.value.start;
//     }

//     private updateSegmentRange() {
//         const sMin = Interval.min(this.setRange), sMax = Interval.max(this.setRange);
//         if (sMax < sMin) {
//             this.hasNext = false;
//             return;
//         }

//         this.segmentMin = this.segmentMap[OrderedSet.getAt(this.set, sMin)];
//         this.segmentMax = this.segmentMap[OrderedSet.getAt(this.set, sMax)];

//         this.hasNext = this.segmentMax >= this.segmentMin;
//     }

//     setSegment(segment: Segs.Segment<T>) {
//         this.setRange = Interval.ofBounds(segment.start, segment.end);
//         this.updateSegmentRange();
//     }

//     constructor(private segments: SortedArray, private segmentMap: Int32Array, private set: OrderedSet, inputRange: Interval) {
//         this.setRange = inputRange;
//         this.updateSegmentRange();
//     }
// }