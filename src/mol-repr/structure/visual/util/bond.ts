/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../../mol-model/structure/model/types';
import { Unit, StructureElement, Structure, Bond } from '../../../../mol-model/structure';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { LocationIterator } from '../../../../mol-geo/util/location-iterator';
import { LinkCylinderParams, LinkLineParams } from './link';
import { ObjectKeys } from '../../../../mol-util/type-helpers';
import { PickingId } from '../../../../mol-geo/geometry/picking';
import { EmptyLoci, Loci } from '../../../../mol-model/loci';
import { Interval, OrderedSet, SortedArray } from '../../../../mol-data/int';
import { isH, isHydrogen, StructureGroup } from './common';

export const BondParams = {
    includeTypes: PD.MultiSelect(ObjectKeys(BondType.Names), PD.objectToOptions(BondType.Names)),
    excludeTypes: PD.MultiSelect([] as BondType.Names[], PD.objectToOptions(BondType.Names)),
    ignoreHydrogens: PD.Boolean(false),
    aromaticBonds: PD.Boolean(false, { description: 'Display aromatic bonds with dashes' }),
    multipleBonds: PD.Select('symmetric', PD.arrayToOptions(['off', 'symmetric', 'offset'] as const)),
};
export const DefaultBondProps = PD.getDefaultValues(BondParams);
export type BondProps = typeof DefaultBondProps

export const BondCylinderParams = {
    ...LinkCylinderParams,
    ...BondParams,
    adjustCylinderLength: PD.Boolean(false, { description: 'Shorten cylinders to reduce overlap with spheres. Useful for for transparent bonds. Not working well with aromatic bonds.' })
};
export const DefaultBondCylinderProps = PD.getDefaultValues(BondCylinderParams);
export type BondCylinderProps = typeof DefaultBondCylinderProps

export const BondLineParams = {
    ...LinkLineParams,
    ...BondParams
};
export const DefaultBondLineProps = PD.getDefaultValues(BondLineParams);
export type BondLineProps = typeof DefaultBondLineProps

export function ignoreBondType(include: BondType.Flag, exclude: BondType.Flag, f: BondType.Flag) {
    return !BondType.is(include, f) || BondType.is(exclude, f);
}

export function makeIntraBondIgnoreTest(structure: Structure, unit: Unit.Atomic, props: BondProps): undefined | ((edgeIndex: number) => boolean) {
    const elements = unit.elements;
    const { atomicNumber } = unit.model.atomicHierarchy.derived.atom;
    const bonds = unit.bonds;
    const { a, b, edgeProps } = bonds;
    const { flags: _flags } = edgeProps;

    const { ignoreHydrogens, includeTypes, excludeTypes } = props;

    const include = BondType.fromNames(includeTypes);
    const exclude = BondType.fromNames(excludeTypes);
    const allBondTypes = BondType.isAll(include) && BondType.Flag.None === exclude;

    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) throw new Error('expected childUnit to exist if child exists');

    if (allBondTypes && !ignoreHydrogens && !child) return;

    return (edgeIndex: number) => {
        return (
            (!!childUnit && !SortedArray.has(childUnit.elements, elements[a[edgeIndex]])) ||
            (ignoreHydrogens && (isH(atomicNumber, elements[a[edgeIndex]]) || isH(atomicNumber, elements[b[edgeIndex]]))) ||
            (!allBondTypes && ignoreBondType(include, exclude, _flags[edgeIndex]))
        );
    };
}

export function makeInterBondIgnoreTest(structure: Structure, props: BondProps): undefined | ((edgeIndex: number) => boolean) {
    const bonds = structure.interUnitBonds;
    const { edges } = bonds;

    const { ignoreHydrogens, includeTypes, excludeTypes } = props;

    const include = BondType.fromNames(includeTypes);
    const exclude = BondType.fromNames(excludeTypes);
    const allBondTypes = BondType.isAll(include) && BondType.Flag.None === exclude;

    const { child } = structure;

    if (allBondTypes && !ignoreHydrogens && !child) return;

    return (edgeIndex: number) => {
        if (child) {
            const b = edges[edgeIndex];
            const childUnitA = child.unitMap.get(b.unitA);
            if (!childUnitA) return true;

            const unitA = structure.unitMap.get(b.unitA);
            const eA = unitA.elements[b.indexA];
            if (!SortedArray.has(childUnitA.elements, eA)) return true;
        }

        if (ignoreHydrogens) {
            const b = edges[edgeIndex];
            const uA = structure.unitMap.get(b.unitA);
            const uB = structure.unitMap.get(b.unitB);
            if (isHydrogen(uA, uA.elements[b.indexA]) || isHydrogen(uB, uB.elements[b.indexB])) return true;
        }

        if (!allBondTypes) {
            if (ignoreBondType(include, exclude, edges[edgeIndex].props.flag)) return true;
        }

        return false;
    };
}

export namespace BondIterator {
    export function fromGroup(structureGroup: StructureGroup): LocationIterator {
        const { group, structure } = structureGroup;
        const unit = group.units[0] as Unit.Atomic;
        const groupCount = Unit.isAtomic(unit) ? unit.bonds.edgeCount * 2 : 0;
        const instanceCount = group.units.length;
        const location = Bond.Location(structure, undefined, undefined, structure, undefined, undefined);
        const getLocation = (groupIndex: number, instanceIndex: number) => {
            const unit = group.units[instanceIndex] as Unit.Atomic;
            location.aUnit = unit;
            location.bUnit = unit;
            location.aIndex = unit.bonds.a[groupIndex];
            location.bIndex = unit.bonds.b[groupIndex];
            return location;
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation);
    }

    export function fromStructure(structure: Structure): LocationIterator {
        const groupCount = structure.interUnitBonds.edgeCount;
        const instanceCount = 1;
        const location = Bond.Location(structure, undefined, undefined, structure, undefined, undefined);
        const getLocation = (groupIndex: number) => {
            const bond = structure.interUnitBonds.edges[groupIndex];
            location.aUnit = structure.unitMap.get(bond.unitA);
            location.aIndex = bond.indexA;
            location.bUnit = structure.unitMap.get(bond.unitB);
            location.bIndex = bond.indexB;
            return location;
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
    }
}

//

export function getIntraBondLoci(pickingId: PickingId, structureGroup: StructureGroup, id: number) {
    const { objectId, instanceId, groupId } = pickingId;
    if (id === objectId) {
        const { structure, group } = structureGroup;
        const unit = group.units[instanceId];
        if (Unit.isAtomic(unit)) {
            const { target } = structure;
            const iA = unit.bonds.a[groupId];
            const iB = unit.bonds.b[groupId];
            return Bond.Loci(target, [
                Bond.Location(target, unit, iA, target, unit, iB),
                Bond.Location(target, unit, iB, target, unit, iA)
            ]);
        }
    }
    return EmptyLoci;
}

export function eachIntraBond(loci: Loci, structureGroup: StructureGroup, apply: (interval: Interval) => boolean, isMarking: boolean) {
    let changed = false;
    if (Bond.isLoci(loci)) {
        const { structure, group } = structureGroup;
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        const unit = group.units[0];
        if (!Unit.isAtomic(unit)) return false;
        const groupCount = unit.bonds.edgeCount * 2;
        for (const b of loci.bonds) {
            if (b.aUnit !== b.bUnit) continue;
            const unitIdx = group.unitIndexMap.get(b.aUnit.id);
            if (unitIdx !== undefined) {
                const idx = unit.bonds.getDirectedEdgeIndex(b.aIndex, b.bIndex);
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx))) changed = true;
                }
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        const { structure, group } = structureGroup;
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        const unit = group.units[0];
        if (!Unit.isAtomic(unit)) return false;
        const groupCount = unit.bonds.edgeCount * 2;
        for (const e of loci.elements) {
            const unitIdx = group.unitIndexMap.get(e.unit.id);
            if (unitIdx !== undefined) {
                const { offset, b } = unit.bonds;
                OrderedSet.forEach(e.indices, v => {
                    for (let t = offset[v], _t = offset[v + 1]; t < _t; t++) {
                        if (!isMarking || OrderedSet.has(e.indices, b[t])) {
                            if (apply(Interval.ofSingleton(unitIdx * groupCount + t))) changed = true;
                        }
                    }
                });
            }
        }
    }
    return changed;
}

//

export function getInterBondLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const { target } = structure;
        const b = structure.interUnitBonds.edges[groupId];
        const uA = structure.unitMap.get(b.unitA);
        const uB = structure.unitMap.get(b.unitB);
        return Bond.Loci(target, [
            Bond.Location(target, uA, b.indexA, target, uB, b.indexB),
            Bond.Location(target, uB, b.indexB, target, uA, b.indexA)
        ]);
    }
    return EmptyLoci;
}

export function eachInterBond(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean, isMarking: boolean) {
    let changed = false;
    if (Bond.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        for (const b of loci.bonds) {
            const idx = structure.interUnitBonds.getBondIndexFromLocation(b);
            if (idx !== -1) {
                if (apply(Interval.ofSingleton(idx))) changed = true;
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        if (isMarking && loci.elements.length === 1) return false; // only a single unit

        const map = new Map<number, OrderedSet<StructureElement.UnitIndex>>();
        for (const e of loci.elements) map.set(e.unit.id, e.indices);

        for (const e of loci.elements) {
            const { unit } = e;
            if (!Unit.isAtomic(unit)) continue;
            structure.interUnitBonds.getConnectedUnits(unit.id).forEach(b => {
                const otherLociIndices = map.get(b.unitB);
                if (!isMarking || otherLociIndices) {
                    OrderedSet.forEach(e.indices, v => {
                        if (!b.connectedIndices.includes(v)) return;
                        b.getEdges(v).forEach(bi => {
                            if (!isMarking || (otherLociIndices && OrderedSet.has(otherLociIndices, bi.indexB))) {
                                const idx = structure.interUnitBonds.getEdgeIndex(v, unit.id, bi.indexB, b.unitB);
                                if (apply(Interval.ofSingleton(idx))) changed = true;
                            }
                        });
                    });
                }
            });
        }
    }
    return changed;
}