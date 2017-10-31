/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Formats from './model/formats'
import MacromoleculeData from './model/data/macromolecule'
import ConformationData from './model/data/conformation'
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

    macromolecule: MacromoleculeData,
    conformation: ConformationData,

    // used for diffing.
    version: {
        data: number,
        conformation: number
    },

    atomCount: number,
    segments: Readonly<{
        chains: Segmentation,
        residues: Segmentation
    }>
}> { }

export default Model