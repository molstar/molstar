/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Formats from './model/formats'
import CommonInterface from './model/common'
import MacromoleculeInterface from './model/macromolecule'
import Segmentation from '../mol-base/collections/integer/segmentation'

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
interface Model extends Readonly<{
    sourceData: Formats.RawData,

    common: CommonInterface,
    macromolecule: MacromoleculeInterface

    atomCount: number,
    segments: Readonly<{
        chains: Segmentation,
        residues: Segmentation,
        entities: Segmentation
    }>
}> { }

export default Model