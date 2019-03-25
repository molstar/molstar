/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SecondaryStructureType } from '../types';

/** Secondary structure "indexed" by residues. */
interface SecondaryStructure {
    readonly type: ArrayLike<SecondaryStructureType>,

    /** index into the elements array */
    readonly key: ArrayLike<number>,
    /** indexed by key */
    readonly elements: ReadonlyArray<SecondaryStructure.Element>,
    /** string representation of DSSP annotation */
    readonly dsspString: String
}

namespace SecondaryStructure {
    export type Element = None | Bend | Turn | Helix | Sheet

    export interface None {
        kind: 'none'
    }

    export interface Bend {
        kind: 'bend',
        flags: SecondaryStructureType
    }

    export interface Turn {
        kind: 'turn',
        flags: SecondaryStructureType
    }

    export interface Helix {
        kind: 'helix',
        flags: SecondaryStructureType,
        type_id: string, // TODO: use aliased type?
        helix_class: string,
        details?: string
    }

    export interface Sheet {
        kind: 'sheet',
        flags: SecondaryStructureType,
        sheet_id: string,
        symmetry?: string
    }
}

export { SecondaryStructure }