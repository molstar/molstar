/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Formats from './model/formats'
//import CommonProperties from './model/properties/common'
import MacromoleculeProperties from './model/properties/macromolecule'
import Conformation from './model/properties/conformation'
import Segmentation from '../mol-base/collections/integer/segmentation'

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    id: string,

    model_num: number,

    sourceData: Formats.RawData,

    //common: CommonProperties,
    macromolecule: MacromoleculeProperties,
    conformation: Conformation,

    // used for diffing.
    version: {
        data: number,
        conformation: number
    },

    atomCount: number,
    segments: Readonly<{
        chains: Segmentation,
        residues: Segmentation,
        entities: Segmentation
    }>
}> { }

export default Model