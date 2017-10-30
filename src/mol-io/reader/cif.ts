/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import parseText from './cif/text/parser'
import parseBinary from './cif/binary/parser'
import { Frame } from './cif/data-model'
import { toTypedFrame as applySchema } from './cif/schema'
import { Schema as mmCIF_Schema, Frame as mmCIF_Frame } from './cif/schema/mmcif'

export default {
    parseText,
    parseBinary,
    applySchema,
    schema: {
        mmCIF: (frame: Frame) => applySchema<typeof mmCIF_Schema, mmCIF_Frame>(mmCIF_Schema, frame)
    }
}

export * from './cif/data-model'