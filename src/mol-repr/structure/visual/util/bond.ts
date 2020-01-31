/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../../mol-model/structure/model/types';
import { Unit, StructureElement, Structure, Bond } from '../../../../mol-model/structure';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { StructureGroup } from '../../units-visual';
import { LinkCylinderParams } from './link';

export const BondCylinderParams = {
    ...LinkCylinderParams,
    includeTypes: PD.MultiSelect(Object.keys(BondType.Names) as BondType.Names[], PD.objectToOptions(BondType.Names)),
    excludeTypes: PD.MultiSelect([] as BondType.Names[], PD.objectToOptions(BondType.Names)),
}
export const DefaultBondCylinderProps = PD.getDefaultValues(BondCylinderParams)
export type BondCylinderProps = typeof DefaultBondCylinderProps

export function ignoreBondType(include: BondType.Flag, exclude: BondType.Flag, f: BondType.Flag) {
    return !BondType.is(include, f) || BondType.is(exclude, f)
}

export namespace BondIterator {
    export function fromGroup(structureGroup: StructureGroup): LocationIterator {
        const { group } = structureGroup
        const unit = group.units[0]
        const groupCount = Unit.isAtomic(unit) ? unit.bonds.edgeCount * 2 : 0
        const instanceCount = group.units.length
        const location = StructureElement.Location.create()
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex]
            location.unit = unit
            location.element = unit.elements[(unit as Unit.Atomic).bonds.a[groupIndex]]
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation)
    }

    export function fromStructure(structure: Structure): LocationIterator {
        const groupCount = structure.interUnitBonds.edgeCount
        const instanceCount = 1
        const location = Bond.Location()
        const getLocation = (groupIndex: number) => {
            const bond = structure.interUnitBonds.edges[groupIndex]
            location.aUnit = bond.unitA
            location.aIndex = bond.indexA as StructureElement.UnitIndex
            location.bUnit = bond.unitB
            location.bIndex = bond.indexB as StructureElement.UnitIndex
            return location
        }
        return LocationIterator(groupCount, instanceCount, getLocation, true)
    }
}