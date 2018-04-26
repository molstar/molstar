/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Column } from 'mol-data/db'
import { Hierarchy } from './hierarchy';

interface Sequence {
    readonly byEntityKey: { [key: number]: Sequence.Entity }
}

namespace Sequence {
    export interface Entity {
        readonly entityId: string,
        readonly num: Column<number>
        // _entity_poly_seq.mon_id
        readonly compId: Column<string>
    }

    export function fromHierarchy(hierarchy: Hierarchy): Sequence {
        // const { label_comp_id } = hierarchy.residues;

        throw 'not implemented';
    }
}

export default Sequence