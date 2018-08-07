/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db'
import { AtomicHierarchy } from './atomic/hierarchy';
import { Entities } from './common';
import { Sequence } from '../../../sequence';
import { ChainIndex } from '../indexing';

interface StructureSequence {
    readonly sequences: ReadonlyArray<StructureSequence.Entity>,
    readonly byEntityKey: { [key: number]: StructureSequence.Entity }
}

namespace StructureSequence {
    export interface Entity {
        readonly entityId: string,
        readonly num: Column<number>,
        // Corresponds to _entity_poly_seq.mon_id
        readonly compId: Column<string>,
        readonly sequence: Sequence
    }

    export function fromAtomicHierarchy(entities: Entities, hierarchy: AtomicHierarchy, modResMap?: ReadonlyMap<string, string>): StructureSequence {
        const { label_comp_id, label_seq_id } = hierarchy.residues
        const { chainAtomSegments, residueAtomSegments } = hierarchy

        const byEntityKey: StructureSequence['byEntityKey'] = { };
        const sequences: StructureSequence.Entity[] = [];

        for (let cI = 0 as ChainIndex, _cI = hierarchy.chains._rowCount; cI < _cI; cI++) {
            const entityKey = hierarchy.getEntityKey(cI);
            // Only for polymers, trying to mirror _entity_poly_seq
            if (byEntityKey[entityKey] !== void 0 || entities.data.type.value(entityKey) !== 'polymer') continue;

            let start = cI;
            cI++;
            while (cI < _cI && entityKey === hierarchy.getEntityKey(cI) && entities.data.type.value(entityKey) !== 'polymer') {
                cI++;
            }
            cI--;

            const rStart = residueAtomSegments.index[chainAtomSegments.offsets[start]];
            const rEnd = residueAtomSegments.index[chainAtomSegments.offsets[cI + 1]];

            const compId = Column.window(label_comp_id, rStart, rEnd);
            const num = Column.window(label_seq_id, rStart, rEnd);

            byEntityKey[entityKey] = {
                entityId: entities.data.id.value(entityKey),
                compId,
                num,
                sequence: Sequence.ofResidueNames(compId, num, modResMap)
            };

            sequences.push(byEntityKey[entityKey]);
        }

        return { byEntityKey, sequences };
    }
}

export default StructureSequence