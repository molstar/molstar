/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Secondary structure "indexed" by residues. */
export interface SecondaryStructure {
    type: number[],
    index: number[],
    flags: number[],
    /** unique value for each "element". This is because single sheet is speficied by multiple records. */
    key: number[]
}