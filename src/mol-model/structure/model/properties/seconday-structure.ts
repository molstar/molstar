/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SecondaryStructureType } from '../types';

/** Secondary structure "indexed" by residues. */
interface SecondaryStructure {
    // assign flags to each residue
    readonly type: ArrayLike<SecondaryStructureType>,
    /** unique value for each "element". This is because single sheet is speficied by multiple records. */
    readonly key: ArrayLike<number>
}

export { SecondaryStructure }