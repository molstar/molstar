/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO: define property accessor intefaces, graphs, spatial lookups and what have you.

import * as Formats from './model/formats'
import CommonInterface from './model/interfaces/common'
import MacromoleculeInterface from './model/interfaces/common'
import OrderedSet from '../mol-base/collections/ordered-set'

interface Model {
    data: Formats.RawData,

    common: CommonInterface,
    macromolecule: MacromoleculeInterface

    // Atom offsets of the "i-th chain" stored as a closed-open range [chainOffsets[i], chainOffsets[i + 1])
    chainSegments: OrderedSet,
    // Atom offsets of the "i-th residue" stored as a closed-open range [residueOffsets[i], residueOffsets[i + 1])
    residueSegments: OrderedSet,
    // Mapping from a residue index to chain index
    residueChainIndex: ArrayLike<number>,
    // Mapping from an atom into to residue index
    atomResidueIndex: ArrayLike<number>
}

export default Model
