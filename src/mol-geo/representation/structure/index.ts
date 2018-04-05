/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomGroup, AtomSet, Structure, Unit } from 'mol-model/structure';

export interface RepresentationProps {

}

export interface UnitRepresentation {
    create: (unit: Unit, atomGroup: AtomGroup, props: RepresentationProps) => boolean,
    update: (props: RepresentationProps) => boolean,
    draw: () => void
}

// export interface StructureRepresentation {
//     create: (structure: Structure, props: RepresentationProps) => boolean,
//     update: (props: RepresentationProps) => boolean
// }

export class StructureRepresentation {
    constructor () {

    }
    create (structure: Structure, props: RepresentationProps) {
        const { atoms, units } = structure;
        const unitIds = AtomSet.unitIds(atoms);
        for (let i = 0, _i = unitIds.length; i < _i; i++) {
            const unitId = unitIds[i];
            const unit = units[unitId];
            const atomGroup = AtomSet.unitGetByIndex(atoms, i);

        }

        return true
    }
    update: (props: RepresentationProps) => false
}