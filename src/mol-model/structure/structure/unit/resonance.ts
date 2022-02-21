/**
 * Copyright (c) 2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

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
    delocalizedTriplets: {
        /** Return 3rd element in triplet or undefined if `a` and `b` are not part of a triplet */
        getThirdElement: (a: StructureElement.UnitIndex, b: StructureElement.UnitIndex) => StructureElement.UnitIndex | undefined
    }
}

export function getResonance(unit: Unit.Atomic) {
    return {
        delocalizedTriplets: getDelocalizedTriplets(unit)
    };
}

function getDelocalizedTriplets(unit: Unit.Atomic) {
    const bonds = unit.bonds;
    const { b, edgeProps, offset } = bonds;
    const { order: _order, flags: _flags } = edgeProps;
    const { elementAromaticRingIndices } = unit.rings;

    const map = new Map<number, StructureElement.UnitIndex>();

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
            map.set(sortedCantorPairing(i, deloBonds[0]), deloBonds[1]);
            for (let j = 1, jl = deloBonds.length; j < jl; j++) {
                map.set(sortedCantorPairing(i, deloBonds[j]), deloBonds[0]);
            }
        }
    }

    return {
        getThirdElement: (a: StructureElement.UnitIndex, b: StructureElement.UnitIndex) => {
            return map.get(sortedCantorPairing(a, b));
        }
    };
}
