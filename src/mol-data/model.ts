/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Formats from './model/formats'
import HierarchyProperties from './model/properties/hierarchy'
import ConformationProperties from './model/properties/conformation'

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    id: string,

    model_num: number,

    sourceData: Formats.RawData,

    hierarchy: HierarchyProperties,
    conformation: ConformationProperties,

    // used for diffing.
    version: {
        data: number,
        conformation: number
    },

    atomCount: number
}> { }

export default Model