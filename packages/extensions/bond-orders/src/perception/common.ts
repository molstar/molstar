/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { isDebugMode } from '@molstar/core/util/debug';
import type { NumberArray } from '@molstar/core/util/type-helpers';
import type { Unit } from '@molstar/model/model/structure';
import { Elements } from '@molstar/model/model/structure/model/properties/atomic/types';
import { BondType, type ElementSymbol } from '@molstar/model/model/structure/model/types';
import type { UnitIndex } from '@molstar/model/model/structure/structure/element/util';
import { getElementIdx, isHydrogen } from '@molstar/model/model/structure/structure/unit/bonds/common';
import type { IntraUnitBonds } from '@molstar/model/model/structure/structure/unit/bonds/data';
import { altLoc, formalCharge, typeSymbol } from '@molstar/model/props/computed/chemistry/util';
import type { InterUnitNeighbours } from './inter-unit-neighbours.js';

/** The rings of the graph being perceived, from `computeRings(unit, bonds)` and `createIndex(rings, [])`. */
export interface UnitRingsInfo {
  readonly all: ReadonlyArray<ArrayLike<UnitIndex>>;
  readonly elementRingIndices: ReadonlyMap<UnitIndex, ArrayLike<number>>;
}

export interface State {
  unit: Unit.Atomic;
  /** the graph being perceived: a copy of the built-in bonds, written in place */
  bonds: IntraUnitBonds;
  rings: UnitRingsInfo;
  /** number of atoms in the residue conformation (unique altloc) */
  n: number;
  /** UnitIndex per residue atom */
  unitIndices: UnitIndex[];
  /** element symbol per residue atom */
  el: ElementSymbol[];
  /** explicit hydrogen count per residue atom */
  hCount: Int8Array;
  /** formal charge per residue atom */
  charge: Int8Array;
  /** total covalent degree per residue atom, bonds to other units included */
  degree: Int8Array;
  /** intra-unit heavy neighbours (residue based index) per residue atom */
  heavyNeighbours: number[][];
  /** intra-unit all neighbours, explicit Hs included (residue based index) per residue atom */
  allNeighbours: number[][];
  /** Flat xyz of every neighbour of each residue atom — including neighbours in other units */
  neighbourPos: number[][];
  /** Residue-local index of each stored neighbour from `neighbourPos` or -1 when the neighbour is outside this residue */
  neighbourLocal: number[][];
  /** Element of each stored neighbour, parallel to `neighbourPos` / `neighbourLocal` */
  neighbourEl: ElementSymbol[][];
  /**
   * Σ orders of the bonds from each residue atom to atoms outside the residue (another residue of the unit or another
   * unit): fixed during perception, which never writes them (_spec §1.11)
   */
  outsideValence: Int8Array;
  /** 1 if one of the bonds from the residue atom to atoms outside the residue is multiple */
  outsideMultiple: Uint8Array;
  /** geometry (from bond angles) per local atom, as `AtomGeometry` */
  geometry: Int8Array;
  /** remaining open valence per local atom */
  open: Int8Array;
  /** 1 if the atom belongs to a ring within the residue */
  inRing: Uint8Array;
  /** 1 if the atom is a pyramidal nitrogen */
  pyramidalNitrogen: Uint8Array;
  /** `Ambiguity` bit-flags per local atom (geometry that can't be decided) */
  ambiguous: Uint8Array;
  assignedBonds: Uint8Array;
  /** rings of size 5 or 6 that are fully contained in the residue, as local UnitIndex arrays */
  localRings5or6: number[][];
  /**
   * fused ring systems or isolated ring systems (ring components)
   * `atoms`: all atoms in the ring system.
   * `rings`: list of atoms per ring in the system.
   */
  localRingSystems: { atoms: number[]; rings: number[][] }[];
  /**
   * Identifier of the first alternate location found in the residue, default to `''`.
   * This is used to skip bonds to atoms in other alternate locations.
   */
  firstAltLoc: string;
}

const tmpVecA = Vec3();
const tmpVecB = Vec3();

/**
 * Collect the residue's neighbour lists: intra-unit ones from the graph being perceived (`bonds`), neighbours in other
 * units from `interUnit`; neighbours in another altloc than the residue's reference one are skipped. Fills `degree`,
 * `hCount`, the residue-local heavy/all lists and the neighbour-coordinate store (positions and residue-local index,
 * `-1` for a neighbour outside the residue).
 */
export function State(
  unit: Unit.Atomic,
  bonds: IntraUnitBonds,
  rings: UnitRingsInfo,
  interUnit: InterUnitNeighbours,
  start: UnitIndex,
  end: UnitIndex,
): State {
  const { edgeProps, offset, b } = bonds;
  const { order, flags } = edgeProps;
  const { elements, conformation } = unit;
  const { label_alt_id, type_symbol } = unit.model.atomicHierarchy.atoms;

  const assignedBonds = new Uint8Array(order.length);
  const { unitIndices, referenceAltLoc } = filterReferenceAltLoc(unit, start, end);
  const unitIndexToLocalIndex = new Map<UnitIndex, number>(unitIndices.map((u, i) => [u, i]));
  const n = unitIndices.length;
  const el = new Array<ElementSymbol>(n);
  const hCount = new Int8Array(n);
  const charge = new Int8Array(n);
  const degree = new Int8Array(n);
  const inRing = new Uint8Array(n);
  const heavyNeighbours: number[][] = [];
  const allNeighbours: number[][] = [];
  const neighbourPos: number[][] = [];
  const neighbourLocal: number[][] = [];
  const neighbourEl: ElementSymbol[][] = [];
  const outsideValence = new Int8Array(n);
  const outsideMultiple = new Uint8Array(n);

  for (let i = 0; i < n; i++) {
    const u = unitIndices[i];
    el[i] = typeSymbol(unit, u);
    const c = formalCharge(unit, u);
    charge[i] = Number.isFinite(c) ? c : 0;
    const heavyList: number[] = [];
    const allList: number[] = [];
    const posList: number[] = [];
    const localList: number[] = [];
    const elList: ElementSymbol[] = [];

    // inter-unit neighbours
    const partners = interUnit.get(elements[u]);
    if (partners) {
      for (const { element: eI, order: bondOrder } of partners) {
        // with `mergeBySymmetry`, a unit holds several chains of the model structure: a partner inside it is already
        // in the unit's bonds
        if (SortedArray.has(elements, eI)) continue;

        const eAlt = label_alt_id.value(eI);
        if (eAlt && eAlt !== referenceAltLoc) continue;

        degree[i]++;
        conformation.invariantPosition(eI, tmpVecA);
        posList.push(tmpVecA[0], tmpVecA[1], tmpVecA[2]);
        localList.push(-1);
        elList.push(type_symbol.value(eI));
        outsideValence[i] += bondOrder;
        if (bondOrder > 1) outsideMultiple[i] = 1;
      }
    }

    // intra-unit neighbours, from the graph being perceived
    for (let t = offset[u]; t < offset[u + 1]; t++) {
      if (!BondType.isCovalent(flags[t])) continue;
      const v = b[t];
      const vAlt = altLoc(unit, v);
      if (vAlt && vAlt !== referenceAltLoc) continue; // skip bonds to atoms in other altlocs
      degree[i]++;
      conformation.invariantPosition(elements[v], tmpVecA);
      posList.push(tmpVecA[0], tmpVecA[1], tmpVecA[2]);
      const vEl = typeSymbol(unit, v);
      elList.push(vEl);
      const j = unitIndexToLocalIndex.get(v);
      if (j === void 0) {
        // not in this residue: never written during perception, its order is fixed (_spec §1.11)
        localList.push(-1);
        outsideValence[i] += order[t];
        if (order[t] > 1) outsideMultiple[i] = 1;
        continue;
      }
      localList.push(j);
      if (isHydrogenElement(vEl)) hCount[i]++;
      else heavyList.push(j);
      allList.push(j);
    }

    heavyNeighbours.push(heavyList);
    allNeighbours.push(allList);
    neighbourPos.push(posList);
    neighbourLocal.push(localList);
    neighbourEl.push(elList);

    for (let t = offset[u]; t < offset[u + 1]; t++) {
      if (!isPerceivable(flags[t], order[t])) assignedBonds[t] = 1;
    }
  }

  // Only the rings fully inside the residue (_spec §1.11): a ring through an atom outside it (another residue, or
  // another altloc: atoms of different altlocs are never bonded, so each altloc has its own ring) sets no `inRing`, is
  // not a 5/6 ring and is in no ring system.
  const { all: allRings, elementRingIndices } = rings;
  const seenRings = new Set<number>();
  /** the complete rings met from the residue's atoms (every atom in the residue), as local atoms */
  const residueRings: number[][] = [];
  const localRings5or6: number[][] = [];
  for (let i = 0; i < n; i++) {
    const u = unitIndices[i];
    const ringsForAtom = elementRingIndices.get(u);
    if (!ringsForAtom) continue;
    for (let r = 0; r < ringsForAtom.length; r++) {
      const ri = ringsForAtom[r];
      if (seenRings.has(ri)) continue;
      seenRings.add(ri);
      const ring = allRings[ri];
      let local: number[] | undefined = new Array<number>(ring.length);
      for (let k = 0; k < ring.length; k++) {
        const localIndex = unitIndexToLocalIndex.get(ring[k]);
        if (localIndex === void 0) {
          // not in this residue (altloc or different residue)
          local = void 0;
          break;
        }
        local[k] = localIndex;
      }
      if (!local) continue;
      residueRings.push(local);
      for (const k of local) inRing[k] = 1;
      if (ring.length < 5 || ring.length > 6) continue;
      const seq = ringTraversalSort(heavyNeighbours, local);
      if (seq) localRings5or6.push(seq);
    }
  }

  // Ring systems collected for aromaticity perception:
  // the residue's complete rings (`residueRings`) grouped by shared atoms, each connected group a system.
  // This differs from Mol* ringComponents which may merge rings across residues and altlocs. When these kind of rings
  // are filtered out, the remaining rings are grouped by shared atoms or left alone (single ring systems).
  // 3 and 4-membered rings are ignored for aromaticity unless part of a larger system with larger rings
  const ringGraph = new IntAdjacencyGraph.UniqueEdgeBuilder(residueRings.length);
  const ringOfAtom = new Int32Array(n).fill(-1); // first ring "seen" for each atom. Used to link rings when an atom is already assigned to a ring
  for (let r = 0; r < residueRings.length; r++) {
    for (const k of residueRings[r]) {
      if (ringOfAtom[k] >= 0) ringGraph.addEdge(ringOfAtom[k], r);
      else ringOfAtom[k] = r;
    }
  }
  const { componentCount, componentIndex } = IntAdjacencyGraph.connectedComponents(ringGraph.getGraph());
  const localRingSystems: State['localRingSystems'] = [];
  for (let c = 0; c < componentCount; c++) {
    const atomsSet = new Set<number>();
    const systemRings: number[][] = [];
    let hasRingGreaterThan4 = false;
    for (let r = 0; r < residueRings.length; r++) {
      if (componentIndex[r] !== c) continue;
      if (residueRings[r].length > 4) hasRingGreaterThan4 = true;
      for (const k of residueRings[r]) atomsSet.add(k);
      systemRings.push(residueRings[r]);
    }
    if (hasRingGreaterThan4) localRingSystems.push({ atoms: [...atomsSet], rings: systemRings });
  }

  return {
    unit,
    bonds,
    rings,
    n,
    unitIndices,
    firstAltLoc: referenceAltLoc,
    el,
    hCount,
    charge,
    degree,
    heavyNeighbours,
    allNeighbours,
    neighbourPos,
    neighbourLocal,
    neighbourEl,
    outsideValence,
    outsideMultiple,
    inRing,
    geometry: new Int8Array(n),
    open: new Int8Array(n),
    ambiguous: new Uint8Array(n),
    pyramidalNitrogen: new Uint8Array(n),
    assignedBonds,
    localRings5or6,
    localRingSystems,
  };
}

/**
 * order ring atoms indices by cyclic walk
 */
function ringTraversalSort(neighbours: number[][], ring: number[]): number[] | undefined {
  const len = ring.length;

  let current = ring[0],
    prev = -1;
  const seq = [current];
  for (let i = 1; i < len; i++) {
    for (const nb of neighbours[current]) {
      if (!ring.includes(nb) || prev === nb) continue;
      seq.push(nb);
      prev = current;
      current = nb;
      break;
    }
  }
  if (seq.length !== len) return undefined; // shoud never happen
  return seq;
}

function filterReferenceAltLoc(unit: Unit.Atomic, start: UnitIndex, end: UnitIndex) {
  const { label_alt_id } = unit.model.atomicHierarchy.atoms;
  const { elements } = unit;
  const unitIndices: UnitIndex[] = [];
  let referenceAltLoc = '';

  for (let i = start; i < end; ++i) {
    const eI = elements[i];
    const altId = label_alt_id.value(eI);
    if (altId) {
      if (!referenceAltLoc) referenceAltLoc = altId;
      if (altId !== referenceAltLoc) {
        continue;
      }
    }
    unitIndices.push(i);
  }

  return { unitIndices, referenceAltLoc };
}

/**
 * Cases where an atom's sp2/sp3 hybridization can't be decided from coordinates alone.
 */
export const enum Ambiguity {
  None = 0,
  /** degree-2 atom in a non-planar 5-ring: the ~108° angle can't tell sp2 from sp3 */
  Sp2OrSp3InRing5 = 1,
  /** degree-2 atom whose single angle sits just below the sp2 threshold */
  Sp2OrSp3Borderline = 2,
}

/** An order nobody stated: `Computed` on distance bonds, `OrderUnknown` on struct_conn/CONECT bonds without one. */
const UnknownOrder = BondType.Flag.Computed | BondType.Flag.OrderUnknown;

/**
 * A bond is eligible for order perception if it is a single covalent bond whose order nobody stated (`UnknownOrder`).
 * Bonds with an explicit order (`chem_comp_bond`, struct_conn `doub`/`sing`, mol/mol2/sdf) carry neither flag and are
 * left untouched.
 */
export function isPerceivable(flags: number, order: number) {
  return order === 1 && (flags & UnknownOrder) !== 0 && BondType.isCovalent(flags);
}

/**
 * Whether a bond inside the residue (both atoms from `start` to `end`, unit indices) is perceivable. Bonds to other
 * residues do not count: an unknown peptide or glycosidic bond (found by distance) does not make a residue with stated
 * orders perceived (_spec §1.9).
 */
export function hasPerceivableBond(bonds: IntraUnitBonds, start: UnitIndex, end: UnitIndex) {
  const { offset, b, edgeProps } = bonds;
  const { order, flags } = edgeProps;
  for (let t = offset[start], tl = offset[end]; t < tl; t++) {
    const v = b[t];
    if (v >= start && v < end && isPerceivable(flags[t], order[t])) return true;
  }
  return false;
}

/**
 * Position of residue atom `i`. Positions are always the model's (`invariantPosition`), never through the unit's
 * operator: the result is one of the model and the elements (_spec §1.3).
 */
export function pos(state: State, i: number, out: Vec3) {
  const { elements, conformation } = state.unit;
  return conformation.invariantPosition(elements[state.unitIndices[i]], out);
}

export function distSq(state: State, i: number, j: number) {
  pos(state, i, tmpVecA);
  pos(state, j, tmpVecB);
  return Vec3.squaredDistance(tmpVecA, tmpVecB);
}

export function isBondFromMainAltLoc(state: State, u: UnitIndex, v: UnitIndex) {
  const { firstAltLoc, unit } = state;
  if (!firstAltLoc) return true;
  const uAlt = altLoc(unit, u);
  if (uAlt && uAlt !== firstAltLoc) return false;
  const vAlt = altLoc(unit, v);
  if (vAlt && vAlt !== firstAltLoc) return false;
  return true;
}

/** Set the order (and optionally OR-in flags) of the undirected bond u-v in both directed slots. */
export function setBond(
  bonds: IntraUnitBonds,
  u: UnitIndex,
  v: UnitIndex,
  order: number,
  addFlags: number,
  assignedBonds?: Uint8Array,
) {
  const ord = bonds.edgeProps.order as NumberArray;
  const flg = bonds.edgeProps.flags as NumberArray;
  let eidx = bonds.getDirectedEdgeIndex(u, v);
  if (eidx !== -1) {
    ord[eidx] = order;
    if (addFlags) flg[eidx] |= addFlags;
    if (assignedBonds) assignedBonds[eidx] = 1;
  }
  eidx = bonds.getDirectedEdgeIndex(v, u);
  if (eidx !== -1) {
    ord[eidx] = order;
    if (addFlags) flg[eidx] |= addFlags;
    if (assignedBonds) assignedBonds[eidx] = 1;
  }
  isDebugMode && console.log('Set bond ', u, '-', v, ' to order ', order, ' flags ', addFlags);
}

/** Mark the bond u-v as assigned; without `v`, every intra-unit covalent bond of `u`. */
export function setAssignedBond(state: State, u: UnitIndex, v?: UnitIndex) {
  const { bonds } = state;
  if (v === void 0) {
    const { offset, b, edgeProps } = bonds;
    for (let t = offset[u]; t < offset[u + 1]; t++) {
      if (BondType.isCovalent(edgeProps.flags[t])) setAssignedBond(state, u, b[t] as UnitIndex);
    }
    return;
  }
  let eidx = bonds.getDirectedEdgeIndex(u, v);
  if (eidx !== -1) state.assignedBonds[eidx] = 1;
  eidx = bonds.getDirectedEdgeIndex(v, u);
  if (eidx !== -1) state.assignedBonds[eidx] = 1;
}

export function isAssignedBond(state: State, u: UnitIndex, v: UnitIndex) {
  const eidx = state.bonds.getEdgeIndex(u, v);
  if (eidx === -1) return true; // meaning we don't touch this bond
  return state.assignedBonds[eidx] === 1;
}

/** Whether atom `i` already has an incident multiple (double/triple) bond, in the residue or leaving it. */
export function hasMultipleBond(state: State, i: number) {
  if (state.outsideMultiple[i]) return true;
  const { unitIndices, heavyNeighbours } = state;
  const u = unitIndices[i];
  for (const j of heavyNeighbours[i]) {
    if (getOrder(state.bonds, u, unitIndices[j]) > 1) return true;
  }
  return false;
}

export function getOrder(bonds: IntraUnitBonds, u: UnitIndex, v: UnitIndex) {
  const eidx = bonds.getEdgeIndex(u, v);
  return eidx !== -1 ? bonds.edgeProps.order[eidx] : 0;
}

export function getFlags(bonds: IntraUnitBonds, u: UnitIndex, v: UnitIndex) {
  const eidx = bonds.getEdgeIndex(u, v);
  return eidx !== -1 ? bonds.edgeProps.flags[eidx] : BondType.Flag.None;
}

export function isHydrogenElement(el: ElementSymbol) {
  const elI = getElementIdx(el);
  return isHydrogen(elI);
}

/** Standard (Daylight-like) heavy valence used for open-valence bookkeeping. -1 = unknown. */
function defaultValence(el: ElementSymbol): number {
  switch (el) {
    case Elements.C:
    case Elements.SI:
      return 4;
    case Elements.N:
    case Elements.P:
    case Elements.B:
      return 3;
    case Elements.O:
    case Elements.S:
    case Elements.SE:
      return 2;
    case Elements.F:
    case Elements.CL:
    case Elements.BR:
    case Elements.I:
    case Elements.H:
    case Elements.D:
    case Elements.T:
      return 1;
    default:
      return -1;
  }
}

export function computeOpenValence(state: State) {
  const { n, bonds, unitIndices } = state;
  for (let i = 0; i < n; i++) {
    let val = defaultValence(state.el[i]);
    if (val < 0) {
      state.open[i] = 0;
      continue;
    }
    val += state.charge[i];
    // explicit H, the fixed bonds leaving the residue, and the current orders to heavy neighbours in it
    let used = state.hCount[i] + state.outsideValence[i];
    const u = unitIndices[i];
    for (const j of state.heavyNeighbours[i]) {
      used += getOrder(bonds, u, unitIndices[j]);
    }
    state.open[i] = Math.max(0, val - used);
  }
}
