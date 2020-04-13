/**
 * Copyright (c) 2019 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Type from '../type';
import * as Struct from './structure-query';
import { Arguments, Argument } from '../symbol';
import { symbol } from '../helpers';

const generator = {
    '@header': 'Generators',

    bundleElement: symbol(Arguments.Dictionary({
        // TODO: should we use more universal unit keys? (i.e. based on chain and "operator name")
        groupedUnits: Argument(Type.Any), // SortedArray<number>[],
        set: Argument(Type.Any), // SortedArray<UnitIndex>
        ranges: Argument(Type.Any) // SortedArray<UnitIndex>
    }), Type.Any), // returns BundleElement

    bundle: symbol(Arguments.Dictionary({
        elements: Argument(Type.Any) // BundleElement[]
    }), Struct.Types.ElementSelectionQuery, 'A selection with single structure containing represented by the bundle.'),

    // Use with caution as this is not "state saveable"
    // This query should never be used in any State Transform!
    current: symbol(Arguments.None, Struct.Types.ElementSelectionQuery, 'Current selection provided by the query context. Avoid using this in State Transforms.')
};

export default {
    '@header': 'Internal Queries',
    generator
};