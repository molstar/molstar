/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureSelection } from './query/selection';
import { StructureQuery } from './query/query';
export * from './query/context';
import * as generators from './query/queries/generators';
import * as modifiers from './query/queries/modifiers';
import * as filters from './query/queries/filters';
import * as combinators from './query/queries/combinators';
import * as internal from './query/queries/internal';
import pred from './query/predicates';

export const Queries = {
    generators,
    filters,
    modifiers,
    combinators,
    pred,
    internal
};

export { StructureSelection, StructureQuery };