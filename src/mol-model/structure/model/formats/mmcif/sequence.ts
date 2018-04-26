/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database as mmCIF } from 'mol-io/reader/cif/schema/mmcif'
import Sequence from '../../properties/sequence'
import { Column } from 'mol-data/db';
import { Hierarchy } from '../../properties/hierarchy';

export function getSequence(cif: mmCIF, hierarchy: Hierarchy): Sequence {
    if (!cif.entity_poly_seq._rowCount) return Sequence.fromHierarchy(hierarchy);

    const { entity_id, num, mon_id } = cif.entity_poly_seq;

    const byEntityKey: Sequence['byEntityKey'] = {};
    const count = entity_id.rowCount;

    let i = 0;
    while (i < count) {
        const start = i;
        while (i < count - 1 && entity_id.areValuesEqual(i, i + 1)) i++;
        i++;

        const id = entity_id.value(start);
        byEntityKey[hierarchy.findEntityKey(id)] = { entityId: id, compId: Column.window(mon_id, start, i), num: Column.window(num, start, i)  }
    }

    return { byEntityKey };
}