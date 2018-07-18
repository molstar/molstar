/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomicSegments } from '../atomic';
import { AtomicData, AtomicRanges } from '../atomic/hierarchy';
import { Segmentation, Interval } from 'mol-data/int';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { ChemicalComponent } from '../chemical-component';
import { MoleculeType, isPolymer } from '../../types';
import { ElementIndex } from '../../indexing';

// TODO add gaps at the ends of the chains by comparing to the polymer sequence data

export function getAtomicRanges(data: AtomicData, segments: AtomicSegments, chemicalComponentMap: Map<string, ChemicalComponent>): AtomicRanges {
    const polymerRanges: number[] = []
    const gapRanges: number[] = []
    const chainIt = Segmentation.transientSegments(segments.chainAtomSegments, Interval.ofBounds(0, data.atoms._rowCount))
    const residueIt = Segmentation.transientSegments(segments.residueAtomSegments, Interval.ofBounds(0, data.atoms._rowCount))
    const { label_seq_id, label_comp_id } = data.residues

    let prevSeqId: number
    let prevStart: number
    let prevEnd: number
    let startIndex: number

    while (chainIt.hasNext) {
        const chainSegment = chainIt.move();
        residueIt.setSegment(chainSegment);
        prevSeqId = -1
        prevStart = -1
        prevEnd = -1
        startIndex = -1

        while (residueIt.hasNext) {
            const residueSegment = residueIt.move();
            const residueIndex = residueSegment.index
            const cc = chemicalComponentMap.get(label_comp_id.value(residueIndex))
            const moleculeType = cc ? cc.moleculeType : MoleculeType.unknown
            const seqId = label_seq_id.value(residueIndex)
            if (isPolymer(moleculeType)) {
                if (startIndex !== -1) {
                    if (seqId !== prevSeqId + 1) {
                        polymerRanges.push(startIndex, prevEnd - 1)
                        gapRanges.push(prevStart, residueSegment.end - 1)
                        startIndex = residueSegment.start
                    } else if (!residueIt.hasNext) {
                        polymerRanges.push(startIndex, residueSegment.end - 1)
                    }
                } else {
                    startIndex = residueSegment.start // start polymer
                }
            } else {
                if (startIndex !== -1) {
                    polymerRanges.push(startIndex, prevEnd - 1)
                    startIndex = -1
                }
            }

            prevStart = residueSegment.start
            prevEnd = residueSegment.end
            prevSeqId = seqId
        }
    }

    console.log('polymerRanges', polymerRanges)
    console.log('gapRanges', gapRanges)

    return {
        polymerRanges: SortedRanges.ofSortedRanges(polymerRanges as ElementIndex[]),
        gapRanges: SortedRanges.ofSortedRanges(gapRanges as ElementIndex[])
    }
}