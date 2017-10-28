/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import parseText from './cif/text/parser'
import parseBinary from './cif/binary/parser'
import { Block } from './cif/data-model'
import { toTypedFrame as applySchema } from './cif/schema'
import mmCIF from './cif/schema/mmcif'

export default {
    parseText,
    parseBinary,
    applySchema,
    schema: {
        mmCIF: (block: Block) => applySchema(mmCIF, block)
    }
}

export * from './cif/data-model'