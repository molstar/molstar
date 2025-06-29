/**
 * Copyright (c) 2022-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SortedArray } from '../../../../mol-data/int/sorted-array';
import { sortedCantorPairing } from '../../../../mol-data/util';
import { BondType } from '../../model/types';
import { StructureElement } from '../element';
import { Unit } from '../unit';

export type UnitResonance = {
    /**
     * Lookup for triplets of atoms in delocalized bonds.
     *
     * Does not include triplets that are part of aromatic rings.
     */
    readonly delocalizedTriplets: {
        /** Return 3rd element in triplet or undefined if `a` and `b` are not part of a triplet */
        readonly getThirdElement: (a: StructureElement.UnitIndex, b: StructureElement.UnitIndex) => StructureElement.UnitIndex | undefined
        /** Return index into `triplets` or undefined if `a` is not part of any triplet */
        readonly getTripletIndices: (a: StructureElement.UnitIndex) => number[] | undefined
        readonly triplets: SortedArray<StructureElement.UnitIndex>[]
    }
}

const EmptyUnitResonance: UnitResonance = {
    delocalizedTriplets: {
        getThirdElement: () => undefined,
        getTripletIndices: () => undefined,
        triplets: []
    }
};

export function getResonance(unit: Unit.Atomic): UnitResonance {
    if (Unit.Traits.is(unit.traits, Unit.Trait.Water) || Unit.Traits.is(unit.traits, Unit.Trait.CoarseGrained)) {
        return EmptyUnitResonance;
    }
    return {
        delocalizedTriplets: getDelocalizedTriplets(unit)
    };
}

function getDelocalizedTriplets(unit: Unit.Atomic) {
    const bonds = unit.bonds;
    const { b, edgeProps, offset } = bonds;
    const { order: _order, flags: _flags } = edgeProps;
    const { elementAromaticRingIndices } = unit.rings;

    const triplets: SortedArray<StructureElement.UnitIndex>[] = [];
    const thirdElementMap = new Map<number, StructureElement.UnitIndex>();
    const indicesMap = new Map<number, number[]>();

    const add = (a: StructureElement.UnitIndex, b: StructureElement.UnitIndex, c: StructureElement.UnitIndex) => {
        const index = triplets.length;
        triplets.push(SortedArray.ofUnsortedArray([a, b, c]));
        thirdElementMap.set(sortedCantorPairing(a, b), c);
        if (indicesMap.has(a)) indicesMap.get(a)!.push(index);
        else indicesMap.set(a, [index]);
    };

    for (let i = 0 as StructureElement.UnitIndex; i < unit.elements.length; i++) {
        if (elementAromaticRingIndices.has(i)) continue;

        const count = offset[i + 1] - offset[i] + 1;
        if (count < 2) continue;

        const deloBonds: StructureElement.UnitIndex[] = [];
        for (let t = offset[i], _t = offset[i + 1]; t < _t; t++) {
            const f = _flags[t];
            if (!BondType.is(f, BondType.Flag.Aromatic)) continue;

            deloBonds.push(b[t]);
        }

        if (deloBonds.length >= 2) {
            add(i, deloBonds[0], deloBonds[1]);
            for (let j = 1, jl = deloBonds.length; j < jl; j++) {
                add(i, deloBonds[j], deloBonds[0]);
            }
        }
    }

    return {
        getThirdElement: (a: StructureElement.UnitIndex, b: StructureElement.UnitIndex) => {
            return thirdElementMap.get(sortedCantorPairing(a, b));
        },
        getTripletIndices: (a: StructureElement.UnitIndex) => {
            return indicesMap.get(a);
        },
        triplets,
    };
}
