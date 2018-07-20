/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureSelection } from './query/selection'
import { StructureQuery } from './query/query'
export * from './query/context'
import * as generators from './query/queries/generators'
import * as modifiers from './query/queries/modifiers'
import pred from './query/predicates'

export const Queries = {
    generators,
    modifiers,
    pred
}

export { StructureSelection, StructureQuery }