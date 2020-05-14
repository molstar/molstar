/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../../mol-model/structure/model/types';
import { Unit, StructureElement, Structure } from '../../../../mol-model/structure';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { StructureGroup } from '../../units-visual';
import { LinkCylinderParams } from './link';
import { ObjectKeys } from '../../../../mol-util/type-helpers';

export const BondCylinderParams = {
    ...LinkCylinderParams,
    includeTypes: PD.MultiSelect(ObjectKeys(BondType.Names), PD.objectToOptions(BondType.Names)),
    excludeTypes: PD.MultiSelect([] as BondType.Names[], PD.objectToOptions(BondType.Names)),
};
export const DefaultBondCylinderProps = PD.getDefaultValues(BondCylinderParams);
export type BondCylinderProps = typeof DefaultBondCylinderProps

export function ignoreBondType(include: BondType.Flag, exclude: BondType.Flag, f: BondType.Flag) {
    return !BondType.is(include, f) || BondType.is(exclude, f);
}

export namespace BondIterator {
    export function fromGroup(structureGroup: StructureGroup): LocationIterator {
        const { group, structure } = structureGroup;
        const unit = group.units[0];
        const groupCount = Unit.isAtomic(unit) ? unit.bonds.edgeCount * 2 : 0;
        const instanceCount = group.units.length;
        const location = StructureElement.Location.create(structure);
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex];
            location.unit = unit;
            location.element = unit.elements[(unit as Unit.Atomic).bonds.a[groupIndex]];
            return location;
        };
        return LocationIterator(groupCount, instanceCount, getLocation);
    }

    export function fromStructure(structure: Structure): LocationIterator {
        const groupCount = structure.interUnitBonds.edgeCount;
        const instanceCount = 1;
        const location = StructureElement.Location.create(structure);
        const getLocation = (groupIndex: number) => {
            const bond = structure.interUnitBonds.edges[groupIndex];
            location.unit = bond.unitA;
            location.element = bond.unitA.elements[bond.indexA];
            return location;
        };
        return LocationIterator(groupCount, instanceCount, getLocation, true);
    }
}