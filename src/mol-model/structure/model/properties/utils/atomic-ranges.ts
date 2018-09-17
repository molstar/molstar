/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomicSegments } from '../atomic';
import { AtomicData, AtomicRanges } from '../atomic/hierarchy';
import { Segmentation, Interval } from 'mol-data/int';
import SortedRanges from 'mol-data/int/sorted-ranges';
import { ChemicalComponentMap } from '../chemical-component';
import { MoleculeType, isPolymer } from '../../types';
import { ElementIndex, ResidueIndex } from '../../indexing';
import { getAtomIdForAtomRole } from '../../../util';
import { AtomicConformation } from '../atomic/conformation';
import { Vec3 } from 'mol-math/linear-algebra';

// TODO add gaps at the ends of the chains by comparing to the polymer sequence data

function getMoleculeType(compId: string, chemicalComponentMap: ChemicalComponentMap) {
    const cc = chemicalComponentMap.get(compId)
    return cc ? cc.moleculeType : MoleculeType.unknown
}

function getElementIndexForAtomId(rI: ResidueIndex, atomId: string, data: AtomicData, segments: AtomicSegments): ElementIndex {
    const { offsets } = segments.residueAtomSegments
    const { label_atom_id } = data.atoms
    for (let j = offsets[rI], _j = offsets[rI + 1]; j < _j; j++) {
        if (label_atom_id.value(j) === atomId) return j
    }
    return offsets[rI]
}

function areBackboneConnected(riStart: ResidueIndex, riEnd: ResidueIndex, data: AtomicData, segments: AtomicSegments, conformation: AtomicConformation, chemicalComponentMap: ChemicalComponentMap) {
    const { label_comp_id } = data.residues
    const mtStart = getMoleculeType(label_comp_id.value(riStart), chemicalComponentMap)
    const mtEnd = getMoleculeType(label_comp_id.value(riEnd), chemicalComponentMap)
    if (!isPolymer(mtStart) || !isPolymer(mtEnd)) return false

    const startId = getAtomIdForAtomRole(mtStart, 'backboneStart')
    const endId = getAtomIdForAtomRole(mtEnd, 'backboneEnd')

    const eiStart = getElementIndexForAtomId(riStart, startId, data, segments)
    const eiEnd = getElementIndexForAtomId(riEnd, endId, data, segments)

    const { x, y, z } = conformation
    const pStart = Vec3.create(x[eiStart], y[eiStart], z[eiStart])
    const pEnd = Vec3.create(x[eiEnd], y[eiEnd], z[eiEnd])
    return Vec3.distance(pStart, pEnd) < 10
}

export function getAtomicRanges(data: AtomicData, segments: AtomicSegments, conformation: AtomicConformation, chemicalComponentMap: ChemicalComponentMap): AtomicRanges {
    const polymerRanges: number[] = []
    const gapRanges: number[] = []
    const cyclicPolymerMap = new Map<ResidueIndex, ResidueIndex>()
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

        const riStart = segments.residueAtomSegments.index[chainSegment.start]
        const riEnd = segments.residueAtomSegments.index[chainSegment.end - 1]
        if (areBackboneConnected(riStart, riEnd, data, segments, conformation, chemicalComponentMap)) {
            cyclicPolymerMap.set(riStart, riEnd)
            cyclicPolymerMap.set(riEnd, riStart)
        }

        while (residueIt.hasNext) {
            const residueSegment = residueIt.move();
            const residueIndex = residueSegment.index
            const moleculeType = getMoleculeType(label_comp_id.value(residueIndex), chemicalComponentMap)
            const seqId = label_seq_id.value(residueIndex)
            if (isPolymer(moleculeType)) {
                if (startIndex !== -1) {
                    if (seqId !== prevSeqId + 1) {
                        polymerRanges.push(startIndex, prevEnd - 1)
                        gapRanges.push(prevStart, residueSegment.end - 1)
                        startIndex = residueSegment.start
                    } else if (!residueIt.hasNext) {
                        polymerRanges.push(startIndex, residueSegment.end - 1)
                    } else {
                        const riStart = segments.residueAtomSegments.index[residueSegment.start]
                        const riEnd = segments.residueAtomSegments.index[prevEnd - 1]
                        if (!areBackboneConnected(riStart, riEnd, data, segments, conformation, chemicalComponentMap)) {
                            polymerRanges.push(startIndex, prevEnd - 1)
                            startIndex = residueSegment.start
                        }
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

    return {
        polymerRanges: SortedRanges.ofSortedRanges(polymerRanges as ElementIndex[]),
        gapRanges: SortedRanges.ofSortedRanges(gapRanges as ElementIndex[]),
        cyclicPolymerMap
    }
}