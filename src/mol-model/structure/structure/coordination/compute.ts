/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';
import { Structure } from '../structure';
import { Unit } from '../unit';
import { ElementIndex } from '../../model/indexing';
import { StructureElement } from '../element';
import { Coordination, CoordinationSite, EmptyCoordination } from './data';
import { cantorPairing } from '../../../../mol-data/util';

const _pos = Vec3();

const EmptyArray: ReadonlyArray<number> = [];

/** Minimum number of bonds for an atom to be considered a coordination site */
const MinCoordination = 4;

export function computeCoordination(structure: Structure): Coordination {
    const sites: CoordinationSite[] = [];
    const siteIndexMap = new Map<number, number[]>();

    for (let ui = 0, uil = structure.units.length; ui < uil; ++ui) {
        const unit = structure.units[ui];
        if (!Unit.isAtomic(unit)) continue;

        const { elements } = unit;

        for (let ei = 0, eil = elements.length; ei < eil; ++ei) {
            const element = elements[ei];
            const unitIndex = ei as StructureElement.UnitIndex;
            const ligandPositions: Vec3[] = [];
            const ligandUnits: Unit.Atomic[] = [];
            const ligandElements: ElementIndex[] = [];

            // Intra-unit bonds
            const bonds = unit.bonds;
            const offset = bonds.offset;
            for (let bi = offset[unitIndex]; bi < offset[unitIndex + 1]; ++bi) {
                const neighborIdx = bonds.b[bi];
                const neighborElement = elements[neighborIdx];
                unit.conformation.position(neighborElement, _pos);
                ligandPositions.push(Vec3.clone(_pos));
                ligandUnits.push(unit);
                ligandElements.push(neighborElement);
            }

            // Inter-unit bonds
            const interBondIndices = structure.interUnitBonds.getEdgeIndices(unitIndex, unit.id);
            for (let ib = 0; ib < interBondIndices.length; ++ib) {
                const edge = structure.interUnitBonds.edges[interBondIndices[ib]];
                const otherUnitId = edge.unitA !== unit.id ? edge.unitA : edge.unitB;
                const otherIndex = edge.indexA !== unitIndex ? edge.indexA : edge.indexB;
                const otherUnit = structure.unitMap.get(otherUnitId) as Unit.Atomic;
                const otherElement = otherUnit.elements[otherIndex];
                otherUnit.conformation.position(otherElement, _pos);
                ligandPositions.push(Vec3.clone(_pos));
                ligandUnits.push(otherUnit);
                ligandElements.push(otherElement);
            }

            if (ligandPositions.length >= MinCoordination) {
                const siteIndex = sites.length;
                sites.push({ unit, element, unitIndex, ligandPositions, ligandUnits, ligandElements });

                const key = coordinationKey(unit.id, element);
                if (siteIndexMap.has(key)) {
                    siteIndexMap.get(key)!.push(siteIndex);
                } else {
                    siteIndexMap.set(key, [siteIndex]);
                }
            }
        }
    }

    if (sites.length === 0) return EmptyCoordination;

    return {
        sites,
        getSiteIndices: (unit: Unit.Atomic, element: ElementIndex) => {
            return siteIndexMap.get(coordinationKey(unit.id, element)) || EmptyArray;
        },
    };
}

function coordinationKey(unitId: number, element: ElementIndex) {
    return cantorPairing(unitId, element);
}
