/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database } from 'mol-io/reader/cif/schema/mmcif'
import { AtomicHierarchy } from '../../properties/atomic';
import { Entities } from '../../properties/common';

class ModifiedResidues {
    private residueIndexMap = new Map<number, number>();

    /** All residue indices within the given model that are modified. */
    readonly residueIndices: ReadonlyArray<number>;
    /** Indexed same as residueIndex and is key to ModifiedResidues.data table */
    readonly keys: ReadonlyArray<number>;

    /** Index into the data table. -1 if the residue is not modified. */
    getKey(residueIndex: number): number {
        return this.residueIndexMap.has(residueIndex) ? this.residueIndexMap.get(residueIndex)! : -1;
    }

    constructor(public data: mmCIF_Database['pdbx_struct_mod_residue'], hierarchy: AtomicHierarchy, entities: Entities) {
        if (data._rowCount === 0) {
            this.residueIndices = [];
            this.keys = [];
            return;
        }

        const { PDB_ins_code, auth_seq_id, _rowCount } = data;

        const asym_id = data.label_asym_id.isDefined ? data.label_asym_id : data.auth_asym_id;
        const comp_id = data.label_comp_id.isDefined ? data.label_comp_id : data.auth_comp_id;

        const entityIds = entities.data.id.toArray();

        const residueIndices: number[] = [];
        const keys: number[] = [];

        for (let i = 0; i < _rowCount; i++) {
            const aId = asym_id.value(i);
            const eIdx = getEntityId(hierarchy, entityIds, aId);
            if (eIdx < 0) continue;
            const key = hierarchy.findResidueKey(entityIds[eIdx], aId, comp_id.value(i), auth_seq_id.value(i), PDB_ins_code.value(i));
            if (key >= 0) {
                this.residueIndexMap.set(key, i);
                residueIndices[residueIndices.length] = key;
                keys[keys.length] = i;
            }
        }

        this.residueIndices = residueIndices;
        this.keys = keys;
    }
}

function getEntityId(hierarchy: AtomicHierarchy, ids: ArrayLike<string>, asym_id: string) {
    for (let i = 0, _i = ids.length; i < _i; i++) {
        if (hierarchy.findChainKey(ids[i], asym_id) >= 0) return i;
    }
    return -1;
}

export { ModifiedResidues }