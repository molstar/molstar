/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CoarseRanges, CoarseElementData } from '../coarse/hierarchy';
import { Segmentation, Interval } from 'mol-data/int';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { Element } from '../../../../structure';
import { ChemicalComponent } from '../chemical-component';

// TODO assumes all coarse elements are part of a polymer

export function getCoarseRanges(data: CoarseElementData, chemicalComponentMap: Map<string, ChemicalComponent>): CoarseRanges {
    const polymerRanges: number[] = []
    const gapRanges: number[] = []
    const chainIt = Segmentation.transientSegments(data.chainSegments, Interval.ofBounds(0, data.count))

    const { seq_id_begin, seq_id_end } = data

    while (chainIt.hasNext) {
        const { start, end } = chainIt.move();
        console.log('chain', start, end)

        let startIndex = -1
        let prevSeqEnd = -1
        for (let i = start; i < end; ++i) {
            const seqEnd = seq_id_end.value(i)
            if (i === start) {
                startIndex = i
                prevSeqEnd = seq_id_end.value(i)
            } else {
                if (prevSeqEnd !== seq_id_begin.value(i) - 1) {
                    polymerRanges.push(startIndex, i - 1)
                    startIndex = i
                }
            }
            if (i === end - 1) {
                polymerRanges.push(startIndex, i)
            }
            prevSeqEnd = seqEnd
        }
    }

    console.log(polymerRanges, gapRanges)

    return {
        polymerRanges: SortedRanges.ofSortedRanges(polymerRanges as Element[]),
        gapRanges: SortedRanges.ofSortedRanges(gapRanges as Element[])
    }
}