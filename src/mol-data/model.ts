/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Formats from './model/formats'
import HierarchyProperties from './model/properties/hierarchy'
import ConformationProperties from './model/properties/conformation'
import UUID from '../mol-base/utils/uuid'
import buildMmCIF from './model/builders/mmcif'

export interface Version {
    hierarchy: UUID,
    conformation: UUID
}

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    id: UUID,

    modelNum: number,

    sourceData: Formats.RawData,

    hierarchy: HierarchyProperties,
    conformation: ConformationProperties,

    atomCount: number
}> { }

namespace Model {
    export const ofMmCIF = buildMmCIF;
}

export default Model