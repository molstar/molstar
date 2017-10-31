/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SecondaryStructureType } from '../constants'
import Segmentation from '../../../mol-base/collections/integer/segmentation'

interface Positions {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>
}

interface SecondaryStructure {
    ofResidue: ArrayLike<SecondaryStructureType>,
    // atom segmentation??
    segments: Segmentation
}

interface Conformation {
    positions: Positions,
    secondaryStructure: SecondaryStructure
}

export default Conformation