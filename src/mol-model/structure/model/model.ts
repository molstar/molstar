/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import UUID from 'mol-util/uuid'
import Format from './format'
import StructureSequence from './properties/sequence'
import { AtomicHierarchy, AtomicConformation } from './properties/atomic'
import { ModelSymmetry } from './properties/symmetry'
import { CoarseHierarchy, CoarseConformation } from './properties/coarse'
import { Entities } from './properties/common';
import { CustomProperties } from './properties/custom';
import { SecondaryStructure } from './properties/seconday-structure';

import from_mmCIF from './formats/mmcif'
import { ChemicalComponent } from './properties/chemical-component';

/**
 * Interface to the "source data" of the molecule.
 *
 * "Atoms" are integers in the range [0, atomCount).
 */
export interface Model extends Readonly<{
    id: UUID,
    label: string,

    // for IHM, corresponds to ihm_model_list.model_id
    modelNum: number,

    sourceData: Format,

    symmetry: ModelSymmetry,
    entities: Entities,
    sequence: StructureSequence,

    atomicHierarchy: AtomicHierarchy,
    atomicConformation: AtomicConformation,

    properties: {
        /** secondary structure provided by the input file */
        readonly secondaryStructure: SecondaryStructure,
        /** maps modified residue name to its parent */
        readonly modifiedResidueNameMap: Map<string, string>,
        /** maps asym id to unique serial number */
        readonly asymIdSerialMap: Map<string, number>
        /** maps residue name to `ChemicalComponent` data */
        readonly chemicalComponentMap: Map<string, ChemicalComponent>
    },

    customProperties: CustomProperties,

    /**
     * Not to be accessed directly, each custom property descriptor
     * defines property accessors that use this field to store the data.
     */
    _staticPropertyData: { [name: string]: any },
    _dynamicPropertyData: { [name: string]: any },

    coarseHierarchy: CoarseHierarchy,
    coarseConformation: CoarseConformation
}> {

} { }

export namespace Model {
    export function create(format: Format) {
        switch (format.kind) {
            // case 'gro': return from_gro(format);
            case 'mmCIF': return from_mmCIF(format);
        }
    }
}