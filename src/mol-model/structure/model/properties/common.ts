/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mmCIF_Database as mmCIF } from 'mol-io/reader/cif/schema/mmcif'
import { EntityIndex } from '../indexing';

export interface Entities {
    data: mmCIF['entity'],
    getEntityIndex(id: string): EntityIndex
}