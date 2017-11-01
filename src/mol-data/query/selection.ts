/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Structure from './../structure'

type Selection =
    | Structure // each atom is interpreted as a singleton structure
    | Structure[]

namespace Selection {
    export const structureCount: (sel: Selection) => number = 0 as any;
    export const union: (sel: Selection) => Structure = 0 as any;
    export const getAt: (sel: Selection, i: number) => Structure = 0 as any;

    // TODO: 'structure iterator'
    // TODO: selection builders (linear / unique)
    // TODO: spatial lookup
    // TODO: If all structures in a selection are "singletons", collapse them into a single structure
}

export default Selection