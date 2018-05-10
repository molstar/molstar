/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database as mmCIF } from 'mol-io/reader/cif/schema/mmcif'
import { Column } from 'mol-data/db';

export interface Entities {
    data: mmCIF['entity'],
    getEntityIndex(id: string): number
}

export interface Keys {
    // indicate whether the keys form an increasing sequence and within each chain, sequence numbers
    // are in increasing order.
    // monotonous sequences enable for example faster secondary structure assignment.
    isMonotonous: boolean,

    // assign a key to each chain index
    chainKey: Column<number>,
    // assigne a key to each chain index
    // also index to the Entities table.
    entityKey: Column<number>,

    findChainKey(entityId: string, label_asym_id: string): number
}