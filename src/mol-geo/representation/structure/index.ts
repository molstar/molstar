/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AtomGroup, AtomSet, Structure, Unit } from 'mol-model/structure';
import { RenderObject } from 'mol-gl/renderer';
import { EquivalenceClasses } from 'mol-data/util';
import { OrderedSet } from 'mol-data/int'

export interface RepresentationProps {

}

export interface UnitRepresentation {
    create: (unit: Unit, atomGroup: AtomGroup, props: RepresentationProps) => RenderObject[],
    update: (props: RepresentationProps) => boolean,
}

// export interface StructureRepresentation {
//     create: (structure: Structure, props: RepresentationProps) => boolean,
//     update: (props: RepresentationProps) => boolean
// }

export class StructureRepresentation {
    // map: uint.id -> atomGroup.hashCode[]

    constructor () {

    }
    create (structure: Structure, props: RepresentationProps) {
        const { atoms, units } = structure;
        const uniqueGroups = EquivalenceClasses<number, AtomGroup>(
            AtomGroup.hashCode,
            (a, b) => units[a.id].model.id === units[b.id].model.id && OrderedSet.areEqual(a.atoms, b.atoms));

        for (let i = 0, _i = AtomSet.unitCount(atoms); i < _i; i++) {
            const group = AtomSet.unitGetByIndex(atoms, i);
            uniqueGroups.add(group.id, group);

        }

        //uniqueGroups.groups

        return true
    }
    update (atoms: AtomSet, props: RepresentationProps) {
        // TODO check model.id, conformation.id, unit.id, atomGroup(.hashCode/.areEqual)
        return false
    }
}