/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

// Step 9c per-bond-type reference lengths (Angstrom), from Sayle "PDB: Cruft to Content". A bond
// of the given order is accepted when its length is at or below the reference. Pairs absent here
// (e.g. C=N) are never assigned a localized multiple bond by 9c.

import { sortedCantorPairing } from '../../../../mol-data/util';
import { Vec3 } from '../../../../mol-math/linear-algebra';
import { radToDeg } from '../../../../mol-math/misc';
import { Elements } from '../../../../mol-model/structure/model/properties/atomic/types';
import { ElementSymbol } from '../../../../mol-model/structure/model/types';
import { UnitIndex } from '../../../../mol-model/structure/structure/element/element';
import { getElementIdx } from '../../../../mol-model/structure/structure/unit/bonds/common';
import { isDebugMode } from '../../../../mol-util/debug';
import { atomId, eachBondedAtom, typeSymbol } from '../util';
import { isBondFromMainAltLoc, State } from './common';

// Max deviation from 0°/180° of the substituent dihedral across a candidate double bond for it
// to count as planar (`isBondPlanar`); same tolerance as the 6-ring torsion gate. A real double
// stays nearly coplanar even when distorted (CHR C1=C12 ~10.6°, 8RQ CBD=CBE ~3.75°), while a
// twisted single across two sp2 centres (8RQ CBC-CBD ~27°) or a pyramidal sp3 centre (~±60°) is
// rejected.
const BondPlanarMaxDeg = 15;

// from Labute (2005)
const __UnambiguousSingleBondThresholds: { [e: number]: number | undefined } = {
    // C-C    // C-N    // N-N     // C-O     // N-O     // O-O     // C-Si    // N-Si    // C-P     // O-Si    // N-P     // C-S     // O-P     // N-S     // O-S     // Si-Si   // Si-P    // P-P     // Si-S    // P-S     // S-S     // C-Se    // N-Se    // O-Se    // Si-Se    // P-Se     // S-Se     // Se-Se
    84: 1.54, 98: 1.47, 112: 1.45, 113: 1.43, 128: 1.43, 144: 1.47, 224: 1.86, 245: 1.75, 246: 1.85, 267: 1.63, 268: 1.68, 269: 1.75, 291: 1.57, 292: 1.76, 316: 1.57, 420: 2.36, 450: 2.26, 480: 2.26, 481: 2.15, 512: 2.07, 544: 2.05, 854: 1.97, 895: 1.85, 937: 1.97, 1210: 2.42, 1259: 2.27, 1309: 2.19, 2380: 2.34,
};

/**
 * A bond longer than this is too long to be a double/aromatic bond. Returns -1 when the pair is not listed.
 */
export function getUnambiguousSingleBondThreshold(i: number, j: number) {
    if (i < 0 || j < 0) return -1;
    const r = __UnambiguousSingleBondThresholds[sortedCantorPairing(i, j)];
    if (r === void 0) return -1;
    return r;
}

const carbonylBondMaxLength = 1.28; // C=O
export const carbonylBondMaxLengthSq = carbonylBondMaxLength * carbonylBondMaxLength;

/**
 * Reference maximum length (Zhang 2012) for a multiple bond of the given order between the two
 *  elements, or undefined for an unlisted pair.
 */
function multipleBondMaxLength(a: ElementSymbol, b: ElementSymbol, order: number): number | undefined {
    if (order === 3) {
        if (a === Elements.C && b === Elements.C) return 1.25;
        if (a === Elements.C && b === Elements.N) return 1.22;
        return undefined;
    }
    if (order === 2) {
        if (a === Elements.C) {
            if (b === Elements.C) return 1.38;
            if (b === Elements.O) return carbonylBondMaxLength;
            if (b === Elements.S) return 1.70;
            if (b === Elements.N) return 1.29;
            return undefined;
        }
        if (a === Elements.N) {
            if (b === Elements.N) return 1.32;
            if (b === Elements.O) return 1.24;
            return undefined;
        }
    }
    return undefined;
}

/**
 * Squared maximum length for a multiple bond of the given order between two elements.
 * `tolerance` is optionally added to the reference length.
 */
export function multipleBondMaxSq(elA: ElementSymbol, elB: ElementSymbol, order: number, tolerance = 0): number | undefined {
    const a = elA < elB ? elA : elB;
    const b = elA < elB ? elB : elA;
    const len = multipleBondMaxLength(a, b, order);
    if (len === undefined) return undefined;
    const max = len + tolerance;
    return max * max;
}

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpVecC = Vec3();
const tmpVecD = Vec3();

/** True when the bond u-v is longer than the Labute upper single-bond-length reference for its element
 *  pair (`getDoubleBondMaxLength`) — too long to be a double/aromatic bond. Unlisted pairs (-1) never
 *  exceed. */
export function exceedsMaxDoubleLength(state: State, u: UnitIndex, v: UnitIndex) {
    const { unit } = state;
    const t = getUnambiguousSingleBondThreshold(getElementIdx(typeSymbol(unit,u)), getElementIdx(typeSymbol(unit,v)));
    if (t < 0) return false;
    unit.conformation.position(unit.elements[u], tmpVecA);
    unit.conformation.position(unit.elements[v], tmpVecB);
    // t -= 0.05; // Labute (2005) 0.05 Å tolerance for a bond to be considered "too long" for a double. Disabled for 4GR1/RGS where C7'-O11 is 1.38008 vs threshold 1.43 - 0.05!
    isDebugMode && console.log('Bond ', atomId(unit, u), '-', atomId(unit, v), ' length ', Math.sqrt(Vec3.squaredDistance(tmpVecA, tmpVecB)), ' vs threshold ', t);
    return Vec3.squaredDistance(tmpVecA, tmpVecB) > t * t;
}

/**
 * Labute's conservative planarity test for a candidate π bond: is the bond coplanar with *any* pair of
 * its substituents? Scans every heavy-neighbour pair (n_u≠v, n_v≠u) and returns true as soon as one
 * n_u-u-v-n_v dihedral lies within `BondPlanarMaxDeg` of 0°/180° (the "smallest dihedral" form). A
 * terminal end (no heavy substituent) can't disprove planarity, so returns true. Only when every pair is
 * twisted beyond the threshold is the bond non-planar (→ single).
 */
export function hasPlanarDihedral(state: State, u: UnitIndex, v: UnitIndex) {
    const { structure, unit } = state;
    let planarFound = false;

    unit.conformation.position(unit.elements[u], tmpVecB);
    unit.conformation.position(unit.elements[v], tmpVecC);

    eachBondedAtom(structure, unit, u, (unitB, otherIndex) => {
        if (planarFound) return;
        if (unitB.id === unit.id && otherIndex === v) return;
        if (!isBondFromMainAltLoc(state, unit, u, unitB, otherIndex)) return;
        unitB.conformation.position(unitB.elements[otherIndex], tmpVecA);
        eachBondedAtom(structure, unit, v, (unitD, otherIndex2) => {
            if (planarFound) return;
            if (unitD.id === unit.id && otherIndex2 === u) return;
            if (!isBondFromMainAltLoc(state, unit, v, unitD, otherIndex2)) return;
            unitD.conformation.position(unitD.elements[otherIndex2], tmpVecD);
            const deg = Math.abs(radToDeg(Vec3.dihedralAngle(tmpVecA, tmpVecB, tmpVecC, tmpVecD)));
            if (deg <= BondPlanarMaxDeg || deg >= 180 - BondPlanarMaxDeg) {
                planarFound = true;
                isDebugMode && console.log('Bond ', atomId(unit, u), '-', atomId(unit, v), ' has planar dihedral ', deg.toFixed(2), '° with substituents ', atomId(unitB, otherIndex), '-', atomId(unitD, otherIndex2));
            }
        });
    });
    return planarFound;
}