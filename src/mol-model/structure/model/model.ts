/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import UUID from 'mol-util/uuid'
import Format from './format'
import HierarchyProperties from './properties/hierarchy'
import ConformationProperties from './properties/conformation'
import from_mmCIF from './formats/mmcif'


/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    id: UUID,

    modelNum: number,

    sourceData: Format,

    hierarchy: HierarchyProperties,
    conformation: ConformationProperties,

    atomCount: number
}> { }

namespace Model {
    export function create(format: Format) {
        switch (format.kind) {
            case 'mmCIF': return from_mmCIF(format);
        }
    }
}

export default Model