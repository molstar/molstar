/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Selection from './query/selection'
import Query from './query/query'
import * as generators from './query/generators'
import * as props from './query/properties'

export const Queries = {
    generators,
    props
}

export { Selection, Query }