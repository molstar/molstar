/*
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Schema from '../../data/schema'
import * as Data from '../../data/data'

const str = Schema.Field.str()
const int = Schema.Field.int()
const float = Schema.Field.float()

const header = {
    'title': str,
    'timeInPs': float,
    'numberOfAtoms': int,
    'boxX': float,
    'boxY': float,
    'boxZ': float
}

const atoms = {
    'residueNumber': int,
    'residueName': str,
    'atomName': str,
    'atomNumber': int,
    'x': float,
    'y': float,
    'z': float,
    'vx': float,
    'vy': float,
    'vz': float
}

const schema = { header, atoms };
export default function (block: Data.Block) {
    return Schema.apply(schema, block);
}