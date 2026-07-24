/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Himanshu Raj <himanshuraj6771@gmail.com>
 */
import { AtomicNumbers, ElementAtomWeights } from '../../../mol-model/structure/model/properties/atomic/measures';

// Reverse map: atomic number -> canonical element symbol.
// AtomicNumbers includes isotope aliases (D, T -> 1), skip those so we
// resolve to the canonical symbol ('H') instead of an alias.
const SymbolByAtomicNumber: { [n: number]: string } = {};
for (const symbol in AtomicNumbers) {
    if (symbol === 'D' || symbol === 'T') continue;
    const num = AtomicNumbers[symbol];
    if (num === undefined) continue;
    if (SymbolByAtomicNumber[num] === undefined) {
        SymbolByAtomicNumber[num] = symbol;
    }
}

const ElementMassesByMass: [string, number][] = [];
for (const key in ElementAtomWeights) {
    const atomicNumber = Number(key);
    const mass = ElementAtomWeights[atomicNumber];
    const symbol = SymbolByAtomicNumber[atomicNumber];
    if (mass === undefined || symbol === undefined) continue;
    ElementMassesByMass.push([symbol, mass]);
}
// Sort ascending by mass so getElementSymbolFromMass can early-break during lookup.
ElementMassesByMass.sort((a, b) => a[1] - b[1]);

/**
 * Resolve the closest matching element symbol for a given atomic mass.
 *
 * @param mass atomic mass to resolve, in amu
 * @param tolerance max allowed |mass - standard atomic weight| for a match
 * @returns the closest element symbol within `tolerance`, or `undefined` if none is close enough
 */
export function getElementSymbolFromMass(mass: number, tolerance = 5.0): string | undefined {
    if (!Number.isFinite(mass) || mass <= 0) return undefined;

    let minDiff = Infinity;
    let closestSymbol: string | undefined;

    for (const [symbol, elementMass] of ElementMassesByMass) {
        const diff = Math.abs(elementMass - mass);
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