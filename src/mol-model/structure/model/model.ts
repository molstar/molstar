/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import UUID from 'mol-util/uuid'
import Format from './format'
import StructureSequence from './properties/sequence'
import { AtomicHierarchy, AtomicConformation } from './properties/atomic'
import { ModelSymmetry } from './properties/symmetry'
import { CoarseHierarchy, CoarseConformation } from './properties/coarse'
import { Entities } from './properties/common';
import { SecondaryStructure } from './properties/seconday-structure';

//import from_gro from './formats/gro'
import from_mmCIF from './formats/mmcif'

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    id: UUID,
    label: string,

    modelNum: number,

    sourceData: Format,

    symmetry: ModelSymmetry,
    entities: Entities,
    sequence: StructureSequence,

    atomicHierarchy: AtomicHierarchy,
    atomicConformation: AtomicConformation,

    /** Various parts of the code can "cache" custom properties here */
    properties: { readonly secondaryStructure: SecondaryStructure } & { [customName: string]: any },

    coarseHierarchy: CoarseHierarchy,
    coarseConformation: CoarseConformation
}> {

} { }

namespace Model {
    export function create(format: Format) {
        switch (format.kind) {
            //case 'gro': return from_gro(format);
            case 'mmCIF': return from_mmCIF(format);
        }
    }
}

export default Model