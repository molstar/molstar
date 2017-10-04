/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import parseText from './text/parser'
import parseBinary from './binary/parser'
import { Block } from './data-model'
import { apply as applySchema } from './schema'
import mmCIF from './schema/mmcif'

export default {
    parseText,
    parseBinary,
    applySchema,
    schema: {
        mmCIF: (block: Block) => applySchema(mmCIF, block)
    }
}

export * from './data-model'