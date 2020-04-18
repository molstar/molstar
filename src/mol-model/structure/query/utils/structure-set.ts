/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Structure, Unit, StructureElement } from '../../structure';
import { SortedArray } from '../../../../mol-data/int';
import { StructureSubsetBuilder } from '../../structure/util/subset-builder';

export function structureUnion(source: Structure, structures: Structure[]) {
    if (structures.length === 0) return Structure.Empty;
    if (structures.length === 1) return structures[0];

    const unitMap = new Map<number, StructureElement.Set>();
    const fullUnits = new Set<number>();

    for (const { units } of structures) {
        for (let i = 0, _i = units.length; i < _i; i++) {
            const u = units[i];
            if (unitMap.has(u.id)) {
                // check if there is anything more to union in this particual unit.
                if (fullUnits.has(u.id)) continue;
                const merged = SortedArray.union(unitMap.get(u.id)!, u.elements);
                unitMap.set(u.id, merged);
                if (merged.length === source.unitMap.get(u.id).elements.length) fullUnits.add(u.id);
            } else {
                unitMap.set(u.id, u.elements);
                if (u.elements.length === source.unitMap.get(u.id).elements.length) fullUnits.add(u.id);
            }
        }
    }

    const builder = source.subsetBuilder(true);
    unitMap.forEach(buildUnion, builder);
    return builder.getStructure();
}

function buildUnion(this: StructureSubsetBuilder, elements: StructureElement.Set, id: number) {
    this.setUnit(id, elements);
}

export function structureAreEqual(sA: Structure, sB: Structure): boolean {
    if (sA === sB) return true;

    if (sA.units.length !== sB.units.length) return false;

    const aU = sA.units, bU = sB.unitMap;
    for (let i = 0, _i = aU.length; i < _i; i++) {
        const u = aU[i];
        if (!bU.has(u.id)) return false;
        const v = bU.get(u.id);
        if (!SortedArray.areEqual(u.elements, v.elements)) return false;
    }

    return true;
}

export function structureAreIntersecting(sA: Structure, sB: Structure): boolean {
    if (sA === sB) return true;

    let a, b;
    if (sA.units.length < sB.units.length) {
        a = sA; b = sB;
    } else {
        a = sB; b = sA;
    }

    const aU = a.units, bU = b.unitMap;

    for (let i = 0, _i = aU.length; i < _i; i++) {
        const u = aU[i];
        if (!bU.has(u.id)) continue;
        const v = bU.get(u.id);
        if (SortedArray.areIntersecting(u.elements, v.elements)) return true;
    }

    return false;
}

export function structureIntersect(sA: Structure, sB: Structure): Structure {
    if (sA === sB) return sA;
    if (!structureAreIntersecting(sA, sB)) return Structure.Empty;

    let a, b;
    if (sA.units.length < sB.units.length) {
        a = sA; b = sB;
    } else {
        a = sB; b = sA;
    }

    const aU = a.units, bU = b.unitMap;
    const units: Unit[] = [];

    for (let i = 0, _i = aU.length; i < _i; i++) {
        const u = aU[i];
        if (!bU.has(u.id)) continue;
        const v = bU.get(u.id);
        if (SortedArray.areIntersecting(u.elements, v.elements)) {
            const int = SortedArray.intersect(u.elements, v.elements);
            units[units.length] = u.getChild(int);
        }
    }

    return Structure.create(units, { parent: sA.parent || sB.parent });
}

export function structureSubtract(a: Structure, b: Structure): Structure {
    if (a === b) return Structure.Empty;
    if (!structureAreIntersecting(a, b)) return a;

    const aU = a.units, bU = b.unitMap;
    const units: Unit[] = [];

    for (let i = 0, _i = aU.length; i < _i; i++) {
        const u = aU[i];
        if (!bU.has(u.id)) {
            units[units.length] = u;
            continue;
        }
        const v = bU.get(u.id);
        const sub = SortedArray.subtract(u.elements, v.elements);
        if (sub.length > 0) {
            units[units.length] = u.getChild(sub);
        }
    }

    return Structure.create(units, { parent: a.parent || b.parent });
}