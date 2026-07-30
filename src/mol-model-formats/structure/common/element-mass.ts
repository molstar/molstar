/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Himanshu Raj <himanshuraj6771@gmail.com>
 */
import { ElementAtomWeights } from '../../../mol-model/structure/model/properties/atomic/measures';
import { ElementSymbol, getElementFromAtomicNumber } from '../../../mol-model/structure/model/types';

/**
 * Array of element symbols and their corresponding atomic masses, sorted by mass.
 */
const ElementMassesByMass: [ElementSymbol, number][] = [];
for (const key in ElementAtomWeights) {
    const mass = ElementAtomWeights[Number(key)];
    if (mass !== undefined) {
        ElementMassesByMass.push([getElementFromAtomicNumber(Number(key)), mass]);
    }
}
ElementMassesByMass.sort((a, b) => a[1] - b[1]);

/**
 * Resolve the closest matching element symbol for a given atomic mass.
 *
 * @param mass atomic mass to resolve, in amu
 * @param tolerance max allowed percentage difference between the given mass and the closest element mass
 * @returns the closest element symbol within `tolerance`, or `undefined` if none is close enough
 */
export function getElementSymbolFromMass(mass: number, tolerance = 5.0): ElementSymbol | undefined {
    if (!Number.isFinite(mass) || mass <= 0) return undefined;

    let minDiff = Infinity;
    let closestSymbol: ElementSymbol | undefined;

    for (const [symbol, elementMass] of ElementMassesByMass) {
        const diff = (Math.abs(elementMass - mass) / elementMass) * 100;
        if (diff < minDiff) {
            minDiff = diff;
            closestSymbol = symbol;
        } else if (diff > minDiff) {
            // Array is sorted by mass, so the difference can only increase from here.
            break;
        }
    }

    return minDiff <= tolerance ? closestSymbol : undefined;
}