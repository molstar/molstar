/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Selection from './query/selection'
import Query from './query/query'
import * as generators from './query/generators'
import * as modifiers from './query/modifiers'
import pred from './query/predicates'

export const Queries = {
    generators,
    modifiers,
    pred
}

export { Selection, Query }