/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SecondaryStructureFlag as SSF } from './secondary-structure'
import BitFlags from '../../../mol-base/utils/bit-flags'
import Segmentation from '../../../mol-base/collections/integer/segmentation'

interface Positions {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>
}

interface SecondaryStructure {
    ofResidue: ArrayLike<BitFlags<SSF>>,
    // atom segmentation
    segments: Segmentation
}

interface Conformation {
    positions: Positions,
    secondaryStructureType: SecondaryStructure
}

export default Conformation
