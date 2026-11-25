/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { sortedCantorPairing } from '@molstar/core/data/util';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { radToDeg } from '@molstar/core/math/misc';
import { isDebugMode } from '@molstar/core/util/debug';
import { Elements } from '@molstar/model/model/structure/model/properties/atomic/types';
import type { ElementSymbol } from '@molstar/model/model/structure/model/types';
import type { UnitIndex } from '@molstar/model/model/structure/structure/element/util';
import { getElementIdx } from '@molstar/model/model/structure/structure/unit/bonds/common';
import { atomId, typeSymbol } from '@molstar/model/props/computed/chemistry/util';
import type { State } from './common.js';

// Max deviation from 0°/180° of the substituent dihedral across a candidate double bond for it
// to count as planar (`isBondPlanar`); same tolerance as the 6-ring torsion gate. A real double
// stays nearly coplanar even when distorted (CHR C1=C12 ~10.6°, 8RQ CBD=CBE ~3.75°), while a
// twisted single across two sp2 centres (8RQ CBC-CBD ~27°) or a pyramidal sp3 centre (~±60°) is
// rejected.
const BondPlanarMaxDeg = 15;

// from Labute (2005)
const __UnambiguousSingleBondThresholds: { [e: number]: number | undefined } = {
  84: 1.54, // C-C
  98: 1.47, // C-N
  112: 1.45, // N-N
  113: 1.43, // C-O
  128: 1.43, // N-O
  144: 1.47, // O-O
  224: 1.86, // C-Si
  245: 1.75, // N-Si
  246: 1.85, // C-P
  267: 1.63, // O-Si
  268: 1.68, // N-P
  269: 1.75, // C-S
  291: 1.57, // O-P
  292: 1.76, // N-S
  316: 1.57, // O-S
  420: 2.36, // Si-Si
  450: 2.26, // Si-P
  480: 2.26, // P-P
  481: 2.15, // Si-S
  512: 2.07, // P-S
  544: 2.05, // S-S
  854: 1.97, // C-Se
  895: 1.85, // N-Se
  937: 1.97, // O-Se
  1210: 2.42, // Si-Se
  1259: 2.27, // P-Se
  1309: 2.19, // S-Se
  2380: 2.34, // Se-Se
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
 * Reference maximum length (Sayle 2001, Zhang 2012) for a multiple bond of the given order between the two
 *  elements, or undefined for an unlisted pair: an unlisted pair gets no localized multiple bond in step 9c.
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
      if (b === Elements.S) return 1.7;
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
export function multipleBondMaxSq(
  elA: ElementSymbol,
  elB: ElementSymbol,
  order: number,
  tolerance = 0,
): number | undefined {
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
 *  pair (`getUnambiguousSingleBondThreshold`) — too long to be a double/aromatic bond. Unlisted pairs (-1) never
 *  exceed. */
export function exceedsMaxDoubleLength(state: State, u: UnitIndex, v: UnitIndex) {
  const { unit } = state;
  const t = getUnambiguousSingleBondThreshold(getElementIdx(typeSymbol(unit, u)), getElementIdx(typeSymbol(unit, v)));
  if (t < 0) return false;
  const { elements, conformation } = unit;
  conformation.invariantPosition(elements[u], tmpVecA);
  conformation.invariantPosition(elements[v], tmpVecB);
  // t -= 0.05; // Labute (2005) 0.05 Å tolerance for a bond to be considered "too long" for a double. Disabled for 4GR1/RGS where C7'-O11 is 1.38008 vs threshold 1.43 - 0.05!
  isDebugMode &&
    console.log(
      'Bond ',
      atomId(unit, u),
      '-',
      atomId(unit, v),
      ' length ',
      Math.sqrt(Vec3.squaredDistance(tmpVecA, tmpVecB)),
      ' vs threshold ',
      t,
    );
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
  const { unit, unitIndices, neighbourPos, neighbourLocal } = state;
  // Both endpoints are residue atoms; use their residue-local index to read the neighbour store.
  const localU = unitIndices.indexOf(u);
  const localV = unitIndices.indexOf(v);
  if (localU === -1 || localV === -1) return false;

  const { elements, conformation } = unit;
  conformation.invariantPosition(elements[u], tmpVecB);
  conformation.invariantPosition(elements[v], tmpVecC);

  const npU = neighbourPos[localU],
    nlU = neighbourLocal[localU];
  const npV = neighbourPos[localV],
    nlV = neighbourLocal[localV];

  // Substituent coordinates are already collected (altloc-filtered) on both sides; iterate the
  // stored positions, skipping the bond partner (v on u's side, u on v's side).
  for (let a = 0, al = nlU.length; a < al; a++) {
    if (nlU[a] === localV) continue;
    Vec3.set(tmpVecA, npU[3 * a], npU[3 * a + 1], npU[3 * a + 2]);
    for (let b = 0, bl = nlV.length; b < bl; b++) {
      if (nlV[b] === localU) continue;
      Vec3.set(tmpVecD, npV[3 * b], npV[3 * b + 1], npV[3 * b + 2]);
      const deg = Math.abs(radToDeg(Vec3.dihedralAngle(tmpVecA, tmpVecB, tmpVecC, tmpVecD)));
      if (deg <= BondPlanarMaxDeg || deg >= 180 - BondPlanarMaxDeg) {
        isDebugMode &&
          console.log('Bond ', atomId(unit, u), '-', atomId(unit, v), ' has planar dihedral ', deg.toFixed(2), '°');
        return true;
      }
    }
  }
  return false;
}
