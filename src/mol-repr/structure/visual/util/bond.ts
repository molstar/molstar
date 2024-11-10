/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
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
import { isHydrogen, StructureGroup } from './common';

export const BondParams = {
    includeTypes: PD.MultiSelect(ObjectKeys(BondType.Names), PD.objectToOptions(BondType.Names)),
    excludeTypes: PD.MultiSelect([] as BondType.Names[], PD.objectToOptions(BondType.Names)),
    ignoreHydrogens: PD.Boolean(false),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    aromaticBonds: PD.Boolean(true, { description: 'Display aromatic bonds with dashes' }),
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
    const bonds = unit.bonds;
    const { a, b, edgeProps } = bonds;
    const { flags: _flags } = edgeProps;

    const { ignoreHydrogens, ignoreHydrogensVariant, includeTypes, excludeTypes } = props;

    const include = BondType.fromNames(includeTypes);
    const exclude = BondType.fromNames(excludeTypes);
    const allBondTypes = BondType.isAll(include) && BondType.Flag.None === exclude;

    const { child } = structure;
    const childUnit = child?.unitMap.get(unit.id);
    if (child && !childUnit) throw new Error('expected childUnit to exist if child exists');

    if (allBondTypes && !ignoreHydrogens && !child) return;

    return (edgeIndex: number) => {
        const aI = a[edgeIndex];
        const bI = b[edgeIndex];

        if ((!!childUnit && !SortedArray.has(childUnit.elements, elements[aI]))) {
            return true;
        }

        if (!allBondTypes && ignoreBondType(include, exclude, _flags[edgeIndex])) {
            return true;
        }

        if (!ignoreHydrogens) return false;

        if (isHydrogen(structure, unit, elements[aI], ignoreHydrogensVariant) || isHydrogen(structure, unit, elements[bI], ignoreHydrogensVariant)) return true;

        return false;
    };
}

export function makeInterBondIgnoreTest(structure: Structure, props: BondProps): undefined | ((edgeIndex: number) => boolean) {
    const bonds = structure.interUnitBonds;
    const { edges } = bonds;

    const { ignoreHydrogens, ignoreHydrogensVariant, includeTypes, excludeTypes } = props;

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
            if (isHydrogen(structure, uA, uA.elements[b.indexA], ignoreHydrogensVariant) || isHydrogen(structure, uB, uB.elements[b.indexB], ignoreHydrogensVariant)) return true;
        }

        if (!allBondTypes) {
            if (ignoreBondType(include, exclude, edges[edgeIndex].props.flag)) return true;
        }

        return false;
    };
}

export function hasUnitVisibleBonds(unit: Unit.Atomic, props: { ignoreHydrogens: boolean, ignoreHydrogensVariant: 'all' | 'non-polar' | 'polar' }) {
    if (Unit.Traits.is(unit.traits, Unit.Trait.Water)) {
        return !props.ignoreHydrogens || props.ignoreHydrogensVariant === 'non-polar';
    }
    return true;
}

export function hasStructureVisibleBonds(structure: Structure, props: { ignoreHydrogens: boolean, ignoreHydrogensVariant: 'all' | 'non-polar' | 'polar' }) {
    for (const { units } of structure.unitSymmetryGroups) {
        if (Unit.isAtomic(units[0]) && hasUnitVisibleBonds(units[0], props)) return true;
    }
    return false;
}

export namespace BondIterator {
    export function fromGroup(structureGroup: StructureGroup, props?: { includeLocation2?: boolean }): LocationIterator {
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
        if (props?.includeLocation2) {
            const location2 = Bond.Location(structure, undefined, undefined, structure, undefined, undefined);
            const getLocation2 = (groupIndex: number, instanceIndex: number) => { // swapping A and B
                const unit = group.units[instanceIndex] as Unit.Atomic;
                location2.aUnit = unit;
                location2.bUnit = unit;
                location2.aIndex = unit.bonds.b[groupIndex];
                location2.bIndex = unit.bonds.a[groupIndex];
                return location2;
            };
            return LocationIterator(groupCount, instanceCount, 1, getLocation, false, () => false, getLocation2);
        }
        return LocationIterator(groupCount, instanceCount, 1, getLocation);
    }

    export function fromStructure(structure: Structure, props?: { includeLocation2?: boolean }): LocationIterator {
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
        if (props?.includeLocation2) {
            const location2 = Bond.Location(structure, undefined, undefined, structure, undefined, undefined);
            const getLocation2 = (groupIndex: number) => { // swapping A and B
                const bond = structure.interUnitBonds.edges[groupIndex];
                location2.aUnit = structure.unitMap.get(bond.unitB);
                location2.aIndex = bond.indexB;
                location2.bUnit = structure.unitMap.get(bond.unitA);
                location2.bIndex = bond.indexA;
                return location2;
            };
            return LocationIterator(groupCount, instanceCount, 1, getLocation, true, () => false, getLocation2);
        };
        return LocationIterator(groupCount, instanceCount, 1, getLocation, true);
    }

    export function fromStructureGroups(structure: Structure, props?: { includeLocation2?: boolean }): LocationIterator {
        const { bondCount, unitIndex, unitGroupIndex, unitEdgeIndex } = structure.intraUnitBondMapping;
        const groupCount = bondCount;
        const instanceCount = 1;
        const location = Bond.Location(structure, undefined, undefined, structure, undefined, undefined);
        const getLocation = (groupIndex: number) => {
            const ug = structure.unitSymmetryGroups[unitIndex[groupIndex]];
            const unit = ug.units[unitGroupIndex[groupIndex]] as Unit.Atomic;
            const edgeIndex = unitEdgeIndex[groupIndex];
            location.aUnit = unit;
            location.bUnit = unit;
            location.aIndex = unit.bonds.a[edgeIndex];
            location.bIndex = unit.bonds.b[edgeIndex];
            return location;
        };
        if (props?.includeLocation2) {
            const location2 = Bond.Location(structure, undefined, undefined, structure, undefined, undefined);
            const getLocation2 = (groupIndex: number) => { // swapping A and B
                const ug = structure.unitSymmetryGroups[unitIndex[groupIndex]];
                const unit = ug.units[unitGroupIndex[groupIndex]] as Unit.Atomic;
                const edgeIndex = unitEdgeIndex[groupIndex];
                location2.aUnit = unit;
                location2.bUnit = unit;
                location2.aIndex = unit.bonds.b[edgeIndex];
                location2.bIndex = unit.bonds.a[edgeIndex];
                return location2;
            };
            return LocationIterator(groupCount, instanceCount, 1, getLocation, true, () => false, getLocation2);
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

const __unitMap = new Map<number, OrderedSet<StructureElement.UnitIndex>>();

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

        for (const e of loci.elements) __unitMap.set(e.unit.id, e.indices);

        for (const e of loci.elements) {
            const { unit } = e;
            if (!Unit.isAtomic(unit)) continue;
            for (const b of structure.interUnitBonds.getConnectedUnits(unit.id)) {
                const otherLociIndices = __unitMap.get(b.unitB);
                if (!isMarking || otherLociIndices) {
                    OrderedSet.forEach(e.indices, v => {
                        if (!b.connectedIndices.includes(v)) return;
                        for (const bi of b.getEdges(v)) {
                            if (!isMarking || (otherLociIndices && OrderedSet.has(otherLociIndices, bi.indexB))) {
                                const idx = structure.interUnitBonds.getEdgeIndex(v, unit.id, bi.indexB, b.unitB);
                                if (apply(Interval.ofSingleton(idx))) changed = true;
                            }
                        }
                    });
                }
            }
        }

        __unitMap.clear();
    }
    return changed;
}

//

export function getStructureGroupsBondLoci(pickingId: PickingId, structure: Structure, id: number) {
    const { objectId, groupId } = pickingId;
    if (id === objectId) {
        const mapping = structure.intraUnitBondMapping;
        return getIntraBondLoci({
            objectId,
            instanceId: mapping.unitGroupIndex[groupId],
            groupId: mapping.unitEdgeIndex[groupId]
        }, {
            structure,
            group: structure.unitSymmetryGroups[mapping.unitIndex[groupId]]
        }, id);
    }
    return EmptyLoci;
}

export function eachStructureGroupsBond(loci: Loci, structure: Structure, apply: (interval: Interval) => boolean, isMarking: boolean) {
    const { unitGroupOffset } = structure.intraUnitBondMapping;

    let changed = false;
    if (Bond.isLoci(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        for (const b of loci.bonds) {
            if (b.aUnit !== b.bUnit) continue;
            const groupIdx = structure.unitSymmetryGroupsIndexMap.get(b.aUnit.id);
            const group = structure.unitSymmetryGroups[groupIdx];
            const unit = group.units[0];
            if (!Unit.isAtomic(unit)) continue;
            const o = unitGroupOffset[groupIdx];
            const groupCount = unit.bonds.edgeCount * 2;
            const unitIdx = group.unitIndexMap.get(b.aUnit.id);
            if (unitIdx !== undefined) {
                const idx = unit.bonds.getDirectedEdgeIndex(b.aIndex, b.bIndex);
                if (idx !== -1) {
                    if (apply(Interval.ofSingleton(unitIdx * groupCount + idx + o))) changed = true;
                }
            }
        }
    } else if (StructureElement.Loci.is(loci)) {
        if (!Structure.areEquivalent(loci.structure, structure)) return false;
        for (const e of loci.elements) {
            const groupIdx = structure.unitSymmetryGroupsIndexMap.get(e.unit.id);
            const group = structure.unitSymmetryGroups[groupIdx];
            const unit = group.units[0];
            if (!Unit.isAtomic(unit)) continue;
            const o = unitGroupOffset[groupIdx];
            const groupCount = unit.bonds.edgeCount * 2;
            const unitIdx = group.unitIndexMap.get(e.unit.id);
            if (unitIdx !== undefined) {
                const { offset, b } = unit.bonds;
                OrderedSet.forEach(e.indices, v => {
                    for (let t = offset[v], _t = offset[v + 1]; t < _t; t++) {
                        if (!isMarking || OrderedSet.has(e.indices, b[t])) {
                            if (apply(Interval.ofSingleton(unitIdx * groupCount + t + o))) changed = true;
                        }
                    }
                });
            }
        }
    }
    return changed;
}
