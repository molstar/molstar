/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Formats from './model/formats'
import CommonInterface from './model/interfaces/common'
import MacromoleculeInterface from './model/interfaces/common'
import Segmentation from '../mol-base/collections/integer/segmentation'

interface Model {
    data: Formats.RawData,

    common: CommonInterface,
    macromolecule: MacromoleculeInterface

    chains: Segmentation,
    residues: Segmentation
}

export default Model