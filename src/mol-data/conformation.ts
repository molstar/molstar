/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from '../mol-base/collections/ordered-set'

interface Conformation {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,

    // Assign a secondary structure type to each residue.
    secondaryStructureType: ArrayLike<any>,
    secondaryStructureAtomOffsets: OrderedSet
}

export default Conformation
