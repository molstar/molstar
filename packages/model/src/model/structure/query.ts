/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureSelection } from './query/selection.js';
import { StructureQuery } from './query/query.js';
export * from './query/context.js';
import * as generators from './query/queries/generators.js';
import * as modifiers from './query/queries/modifiers.js';
import * as filters from './query/queries/filters.js';
import * as combinators from './query/queries/combinators.js';
import * as internal from './query/queries/internal.js';
import * as atomset from './query/queries/atom-set.js';
import { Predicates as pred } from './query/predicates.js';

export const Queries = {
    generators,
    filters,
    modifiers,
    combinators,
    pred,
    internal,
    atomset
};

export { StructureSelection, StructureQuery };
