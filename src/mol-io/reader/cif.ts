/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import parseText from './cif/text/parser'
import parseBinary from './cif/binary/parser'
import { Frame } from './cif/data-model'
import { toDatabase } from './cif/schema'
import { mmCIF_Schema, mmCIF_Database } from './cif/schema/mmcif'
import { CCD_Schema, CCD_Database } from './cif/schema/ccd'

export default {
    parseText,
    parseBinary,
    toDatabase,
    schema: {
        mmCIF: (frame: Frame) => toDatabase<mmCIF_Schema, mmCIF_Database>(mmCIF_Schema, frame),
        CCD: (frame: Frame) => toDatabase<CCD_Schema, CCD_Database>(CCD_Schema, frame)
    }
}

export * from './cif/data-model'