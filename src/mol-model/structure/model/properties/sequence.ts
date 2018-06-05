/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db'
import { AtomicHierarchy } from './atomic/hierarchy';
import { Entities } from './common';

interface Sequence {
    readonly byEntityKey: { [key: number]: Sequence.Entity }
}

// TODO lift to model/sequence/ folder
// TODO add one letter code sequence string
// TODO add mapping support to other sequence spaces, e.g. uniprot

namespace Sequence {
    export interface Entity {
        readonly entityId: string,
        readonly num: Column<number>
        // _entity_poly_seq.mon_id
        readonly compId: Column<string>
    }

    export function fromAtomicHierarchy(entities: Entities, hierarchy: AtomicHierarchy): Sequence {
        const { label_entity_id } = hierarchy.chains
        const { label_comp_id, label_seq_id } = hierarchy.residues
        const { chainSegments, residueSegments } = hierarchy

        const byEntityKey: Sequence['byEntityKey'] = {};

        // TODO get min/max of label_seq_id to handle missing residues at start and in between
        //   note that this assumes label_seq_id is monotonically increasing

        const chainCount = hierarchy.chains._rowCount
        for (let i = 0; i < chainCount; ++i) {
            const entityId = label_entity_id.value(i)
            const entityIndex = entities.getEntityIndex(entityId)
            // TODO only for polymers, mirroring _entity_poly_seq, ok???
            if (entities.data.type.value(i) !== 'polymer') continue

            const entityKey = hierarchy.entityKey[entityIndex]
            if (byEntityKey[entityKey] !== undefined) continue

            const start = residueSegments.segmentMap[chainSegments.segments[i]]
            let end = residueSegments.segmentMap[chainSegments.segments[i + 1]]
            // TODO better way to handle end???
            if (end === undefined) end = hierarchy.residues._rowCount

            byEntityKey[entityKey] = {
                entityId,
                compId: Column.window(label_comp_id, start, end),
                num: Column.window(label_seq_id, start, end)
            }
        }

        return { byEntityKey }
    }
}

export default Sequence