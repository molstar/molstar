/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SecondaryStructureType } from '../types';
import { ResidueIndex } from '../indexing';

/** Secondary structure "indexed" by residues. */
interface SecondaryStructure {
    readonly type: ArrayLike<SecondaryStructureType>,
    /** index into the elements array */
    readonly key: ArrayLike<number>,
    /** indexed by key */
    readonly elements: ReadonlyArray<SecondaryStructure.Element>
    /** mapping from residue index */
    readonly getIndex: (rI: ResidueIndex) => number,
}

function SecondaryStructure(type: SecondaryStructure['type'], key: SecondaryStructure['key'], elements: SecondaryStructure['elements'], getIndex: SecondaryStructure['getIndex']) {
    return { type, key, elements, getIndex };
}

namespace SecondaryStructure {
    export type Element = None | Turn | Helix | Sheet

    export interface None {
        kind: 'none'
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

export { SecondaryStructure };