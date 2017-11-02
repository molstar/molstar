/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import DataFormat from './data-format'
import HierarchyProperties from './properties/hierarchy'
import ConformationProperties from './properties/conformation'
import UUID from '../../mol-base/utils/uuid'

import buildMmCIF from './builders/mmcif'

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    id: UUID,

    modelNum: number,

    sourceData: DataFormat,

    hierarchy: HierarchyProperties,
    conformation: ConformationProperties,

    atomCount: number
}> { }

namespace Model {
    export function create(format: DataFormat) {
        switch (format.kind) {
            case 'mmCIF': return buildMmCIF(format);
        }
    }
}

export default Model