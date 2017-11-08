/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import parseText from './cif/text/parser'
import parseBinary from './cif/binary/parser'
import { Frame } from './cif/data-model'
import { toDatabase } from './cif/schema'
import { mmCIF_Schema, mmCIF_Database } from './cif/schema/mmcif'

export default {
    parseText,
    parseBinary,
    toDatabase,
    schema: {
        mmCIF: (frame: Frame) => toDatabase<mmCIF_Schema, mmCIF_Database>(mmCIF_Schema, frame)
    }
}

export * from './cif/data-model'