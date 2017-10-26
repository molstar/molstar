/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

interface Conformation {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,

    // Assign a secondary structure type to each residue.
    secondaryStructureType: ArrayLike<any>
}

export default Conformation
