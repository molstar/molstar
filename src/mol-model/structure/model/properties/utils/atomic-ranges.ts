/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomicRanges, AtomicIndex, AtomicHierarchy, AtomicDerivedData } from '../atomic/hierarchy';
import { Segmentation, Interval } from '../../../../../mol-data/int';
import SortedRanges from '../../../../../mol-data/int/sorted-ranges';
import { isPolymer, PolymerType } from '../../types';
import { ElementIndex, ResidueIndex } from '../../indexing';
import { getAtomIdForAtomRole } from '../../../util';
import { AtomicConformation } from '../atomic/conformation';
import { Vec3 } from '../../../../../mol-math/linear-algebra';
import { Entities } from '../common';
import StructureSequence from '../sequence';

function areBackboneConnected(riStart: ResidueIndex, riEnd: ResidueIndex, conformation: AtomicConformation, index: AtomicIndex, derived: AtomicDerivedData) {
    const { polymerType, traceElementIndex, directionFromElementIndex, directionToElementIndex } = derived.residue;
    const ptStart = polymerType[riStart];
    const ptEnd = polymerType[riEnd];
    if (ptStart === PolymerType.NA || ptEnd === PolymerType.NA) return false;
    if (traceElementIndex[riStart] === -1 || traceElementIndex[riEnd] === -1) return false;

    let eiStart = index.findAtomsOnResidue(riStart, getAtomIdForAtomRole(ptStart, 'backboneStart'));
    let eiEnd = index.findAtomsOnResidue(riEnd, getAtomIdForAtomRole(ptEnd, 'backboneEnd'));

    if (eiStart === -1 || eiEnd === -1) {
        eiStart = index.findAtomsOnResidue(riStart, getAtomIdForAtomRole(ptStart, 'coarseBackbone'));
        eiEnd = index.findAtomsOnResidue(riEnd, getAtomIdForAtomRole(ptEnd, 'coarseBackbone'));
    }

    const { x, y, z } = conformation;
    const pStart = Vec3.create(x[eiStart], y[eiStart], z[eiStart]);
    const pEnd = Vec3.create(x[eiEnd], y[eiEnd], z[eiEnd]);
    const isCoarse = directionFromElementIndex[riStart] === -1 || directionToElementIndex[riStart] === -1 || directionFromElementIndex[riEnd] === -1 || directionToElementIndex[riEnd] === -1;
    return Vec3.distance(pStart, pEnd) < (isCoarse ? 10 : 3);
}

export function getAtomicRanges(hierarchy: AtomicHierarchy, entities: Entities, conformation: AtomicConformation, sequence: StructureSequence): AtomicRanges {
    const polymerRanges: number[] = [];
    const gapRanges: number[] = [];
    const cyclicPolymerMap = new Map<ResidueIndex, ResidueIndex>();
    const chainIt = Segmentation.transientSegments(hierarchy.chainAtomSegments, Interval.ofBounds(0, hierarchy.atoms._rowCount));
    const residueIt = Segmentation.transientSegments(hierarchy.residueAtomSegments, Interval.ofBounds(0, hierarchy.atoms._rowCount));
    const { index, derived } = hierarchy;
    const { label_seq_id } = hierarchy.residues;
    const { label_entity_id } = hierarchy.chains;
    const { moleculeType, traceElementIndex } = derived.residue;

    let prevSeqId: number;
    let prevStart: number;
    let prevEnd: number;
    let startIndex: number;

    while (chainIt.hasNext) {
        const chainSegment = chainIt.move();
        residueIt.setSegment(chainSegment);
        prevSeqId = -1;
        prevStart = -1;
        prevEnd = -1;
        startIndex = -1;

        const eI = entities.getEntityIndex(label_entity_id.value(chainSegment.index));
        const seq = sequence.byEntityKey[eI];
        const maxSeqId = seq ? seq.sequence.seqId.value(seq.sequence.seqId.rowCount - 1) : -1;

        // check cyclic peptides, seqIds and distance must be compatible
        const riStart = hierarchy.residueAtomSegments.index[chainSegment.start];
        const riEnd = hierarchy.residueAtomSegments.index[chainSegment.end - 1];
        const seqIdStart = label_seq_id.value(riStart);
        const seqIdEnd = label_seq_id.value(riEnd);
        if (seqIdStart === 1 && seqIdEnd === maxSeqId && conformation.xyzDefined && areBackboneConnected(riStart, riEnd, conformation, index, derived)) {
            cyclicPolymerMap.set(riStart, riEnd);
            cyclicPolymerMap.set(riEnd, riStart);
        }

        while (residueIt.hasNext) {
            const residueSegment = residueIt.move();
            const residueIndex = residueSegment.index;
            const seqId = label_seq_id.value(residueIndex);
            // treat polymers residues that don't have a trace element resolved as gaps
            if (isPolymer(moleculeType[residueIndex]) && traceElementIndex[residueIndex] !== -1) {
                if (startIndex !== -1) {
                    if (seqId !== prevSeqId + 1) {
                        polymerRanges.push(startIndex, prevEnd - 1);
                        gapRanges.push(prevStart, residueSegment.end - 1);
                        startIndex = residueSegment.start;
                    } else if (!residueIt.hasNext) {
                        polymerRanges.push(startIndex, residueSegment.end - 1);
                        // TODO store terminal gaps
                    } else {
                        const riStart = hierarchy.residueAtomSegments.index[residueSegment.start];
                        const riEnd = hierarchy.residueAtomSegments.index[prevEnd - 1];
                        if (conformation.xyzDefined && !areBackboneConnected(riStart, riEnd, conformation, hierarchy.index, derived)) {
                            polymerRanges.push(startIndex, prevEnd - 1);
                            // add gap even for consecutive residues if they are not connected
                            gapRanges.push(prevStart, residueSegment.end - 1);
                            startIndex = residueSegment.start;
                        }
                    }
                } else {
                    startIndex = residueSegment.start; // start polymer
                    // TODO store terminal gaps
                }
            } else {
                if (startIndex !== -1) {
                    polymerRanges.push(startIndex, prevEnd - 1);
                    startIndex = -1;
                }
            }

            prevStart = residueSegment.start;
            prevEnd = residueSegment.end;
            prevSeqId = seqId;
        }
    }

    return {
        polymerRanges: SortedRanges.ofSortedRanges(polymerRanges as ElementIndex[]),
        gapRanges: SortedRanges.ofSortedRanges(gapRanges as ElementIndex[]),
        cyclicPolymerMap
    };
}