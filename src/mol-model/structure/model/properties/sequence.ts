/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Column } from '../../../../mol-data/db';
import { AtomicHierarchy } from './atomic/hierarchy';
import { Entities } from './common';
import { Sequence } from '../../../sequence';
import { ChainIndex } from '../indexing';
import { CoarseHierarchy } from './coarse';
import { CoarseElements } from './coarse/hierarchy';

interface StructureSequence {
    readonly sequences: ReadonlyArray<StructureSequence.Entity>,
    readonly byEntityKey: { [key: number]: StructureSequence.Entity }
}

namespace StructureSequence {
    export interface Entity {
        readonly entityId: string,
        readonly sequence: Sequence
    }

    const Empty: StructureSequence = { byEntityKey: {}, sequences: [] };

    function merge(...entitySeqs: StructureSequence[]): StructureSequence {
        const sequences: StructureSequence.Entity[] = [];
        const byEntityKey: { [key: number]: StructureSequence.Entity } = {};

        for (let i = 0, il = entitySeqs.length; i < il; ++i) {
            sequences.push(...entitySeqs[i].sequences);
            Object.assign(byEntityKey, entitySeqs[i].byEntityKey);
        }
        return { sequences, byEntityKey };
    }

    export function fromHierarchy(entities: Entities, atomicHierarchy: AtomicHierarchy, coarseHierarchy: CoarseHierarchy): StructureSequence {
        const atomic = fromAtomicHierarchy(entities, atomicHierarchy);
        const coarse = coarseHierarchy.isDefined ? fromCoarseHierarchy(entities, coarseHierarchy) : Empty;
        return merge(atomic, coarse);
    }

    export function fromAtomicHierarchy(entities: Entities, hierarchy: AtomicHierarchy): StructureSequence {
        const { label_comp_id } = hierarchy.atoms;
        const { label_seq_id } = hierarchy.residues;
        const { chainAtomSegments, residueAtomSegments } = hierarchy;
        const { count, offsets } = chainAtomSegments;

        const byEntityKey: StructureSequence['byEntityKey'] = { };
        const sequences: StructureSequence.Entity[] = [];

        // check if chain segments are empty
        if (count === 1 && offsets[0] === 0 && offsets[1] === 0) {
            return { byEntityKey, sequences };
        }

        for (let cI = 0 as ChainIndex, _cI = hierarchy.chains._rowCount; cI < _cI; cI++) {
            const entityKey = hierarchy.index.getEntityFromChain(cI);
            // Only for polymers, trying to mirror _entity_poly_seq
            if (byEntityKey[entityKey] !== void 0 || entities.data.type.value(entityKey) !== 'polymer') continue;

            let start = cI;
            cI++;
            while (cI < _cI && entityKey === hierarchy.index.getEntityFromChain(cI) && entities.data.type.value(entityKey) !== 'polymer') {
                cI++;
            }
            cI--;

            const rStart = residueAtomSegments.index[offsets[start]];
            const rEnd = residueAtomSegments.index[offsets[cI + 1] - 1] + 1;
            const seqId = Column.window(label_seq_id, rStart, rEnd);

            const _compId: string[] = [];
            for (let rI = rStart; rI < rEnd; ++rI) {
                _compId.push(label_comp_id.value(residueAtomSegments.offsets[rI]));
            }
            const compId = Column.ofStringArray(_compId);

            byEntityKey[entityKey] = {
                entityId: entities.data.id.value(entityKey),
                sequence: Sequence.ofResidueNames(compId, seqId)
            };

            sequences.push(byEntityKey[entityKey]);
        }

        return { byEntityKey, sequences };
    }

    export function fromCoarseHierarchy(entities: Entities, hierarchy: CoarseHierarchy): StructureSequence {
        const spheres = fromCoarseElements(entities, hierarchy.spheres);
        const gaussians = fromCoarseElements(entities, hierarchy.gaussians);
        return merge(spheres, gaussians);
    }

    export function fromCoarseElements(entities: Entities, elements: CoarseElements): StructureSequence {
        const { chainElementSegments, seq_id_begin, seq_id_end } = elements;
        const { count, offsets } = chainElementSegments;

        const byEntityKey: StructureSequence['byEntityKey'] = { };
        const sequences: StructureSequence.Entity[] = [];

        // check if chain segments are empty
        if (count === 1 && offsets[0] === 0 && offsets[1] === 0) {
            return { byEntityKey, sequences };
        }

        for (let cI = 0 as ChainIndex, _cI = count; cI < _cI; cI++) {
            const eK = elements.getEntityFromChain(cI);
            if (byEntityKey[eK] !== void 0) continue;

            let start = cI;
            cI++;
            while (cI < _cI && eK === elements.getEntityFromChain(cI)) {
                cI++;
            }
            cI--;

            const eStart = offsets[start];
            const eEnd = offsets[cI + 1] - 1;

            const seqIdBegin = Column.window(seq_id_begin, eStart, eEnd);
            const seqIdEnd = Column.window(seq_id_end, eStart, eEnd);

            byEntityKey[eK] = {
                entityId: entities.data.id.value(eK),
                sequence: Sequence.ofSequenceRanges(seqIdBegin, seqIdEnd)
            };

            sequences.push(byEntityKey[eK]);
        }

        return { byEntityKey, sequences };
    }
}

export default StructureSequence;