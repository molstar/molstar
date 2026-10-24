/**
 * Copyright (c) 2017-2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType, MoleculeType } from '../../../model/types.js';
import { Structure } from '../../structure.js';
import { Unit } from '../../unit.js';
import {
  getElementIdx,
  getElementThreshold,
  isHydrogen,
  type BondComputationProps,
  MetalsSet,
  DefaultBondComputationProps,
  getPairingThreshold,
  MaxPairingThreshold,
} from './common.js';
import { InterUnitBonds, type InterUnitEdgeProps } from './data.js';
import { SortedArray } from '@molstar/core/data/int';
import { Vec3, Mat4 } from '@molstar/core/math/linear-algebra';
import { Result } from '@molstar/core/math/geometry/lookup3d/common';
import type { StructureElement } from '../../element.js';
import type { ElementIndex } from '../../../model/indexing.js';
import { getInterBondOrderFromTable } from '../../../model/properties/atomic/bonds.js';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { InterUnitGraph } from '@molstar/core/math/graph/inter-unit-graph';
import { StructConn } from '@molstar/model/formats/structure/property/bonds/struct_conn';
import { equalEps } from '@molstar/core/math/linear-algebra/3d/common';
import { Model } from '../../../model.js';
import { cantorPairing, invertCantorPairing, sortedCantorPairing } from '@molstar/core/data/util';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3distance = Vec3.distance;
const v3set = Vec3.set;
const v3squaredDistance = Vec3.squaredDistance;
const v3transformMat4 = Vec3.transformMat4;

const tmpDistVecA = Vec3();
const tmpDistVecB = Vec3();
function getDistance(unitA: Unit.Atomic, indexA: ElementIndex, unitB: Unit.Atomic, indexB: ElementIndex) {
  unitA.conformation.position(indexA, tmpDistVecA);
  unitB.conformation.position(indexB, tmpDistVecB);
  return v3distance(tmpDistVecA, tmpDistVecB);
}

const _imageTransform = Mat4();
const _imageA = Vec3();
const _bToA = Mat4();
const _bCenterInA = Vec3();
const _candidateResult = Result.create<StructureElement.UnitIndex>();
let _candidates = new Int32Array(1024);

/** All of `unitA`'s elements, in order. */
function allCandidates(unitA: Unit.Atomic) {
  const n = unitA.elements.length;
  if (_candidates.length < n) _candidates = new Int32Array(Math.max(n, 2 * _candidates.length));
  for (let i = 0; i < n; i++) _candidates[i] = i;
  return _candidates.subarray(0, n);
}

/**
 * The elements of `unitA` (sorted UnitIndices) that may lie within `radius` of `unitB`'s boundary sphere: fetched from
 * `unitA`'s grid around the sphere centre mapped into `unitA`'s frame, instead of transforming every element of `unitA`
 * (05-model-structure P4). The radius is padded so that the exact test in `unitB`'s frame still decides. Every element
 * is a candidate, without a grid query (or building `unitA`'s grid), when `unitA`'s sphere lies inside that radius,
 * mostly ions and waters next to a chain, and for operators that do not keep distances (mmCIF allows them with a
 * warning).
 */
function getCandidates(unitA: Unit.Atomic, unitB: Unit.Atomic, radius: number) {
  Mat4.mul(_bToA, unitA.conformation.operator.inverse, unitB.conformation.operator.matrix);
  if (!Mat4.isRigid(_bToA, 1e-4)) return allCandidates(unitA);
  Vec3.transformMat4(_bCenterInA, unitB.boundary.sphere.center, _bToA);
  const sphereA = unitA.boundary.sphere;
  if (Vec3.distance(_bCenterInA, sphereA.center) + sphereA.radius <= radius) return allCandidates(unitA);
  const { indices, count } = unitA.lookup3d.find(
    _bCenterInA[0],
    _bCenterInA[1],
    _bCenterInA[2],
    radius * (1 + 1e-3) + 1e-3,
    _candidateResult,
  );
  if (_candidates.length < count) _candidates = new Int32Array(Math.max(count, 2 * _candidates.length));
  for (let i = 0; i < count; i++) _candidates[i] = indices[i];
  // ascending, as the scan was: keeps the edge insertion order
  return _candidates.subarray(0, count).sort();
}

function findPairBonds(
  unitA: Unit.Atomic,
  unitB: Unit.Atomic,
  props: BondComputationProps,
  builder: InterUnitGraph.Builder<number, StructureElement.UnitIndex, InterUnitEdgeProps>,
) {
  const { maxRadius } = props;
  const distanceRadius = Math.min(maxRadius, MaxPairingThreshold);

  const { elements: atomsA, residueIndex: residueIndexA } = unitA;
  const { x: xA, y: yA, z: zA } = unitA.model.atomicConformation;
  const { elements: atomsB, residueIndex: residueIndexB } = unitB;

  const {
    type_symbol: type_symbolA,
    label_alt_id: label_alt_idA,
    label_atom_id: label_atom_idA,
    label_comp_id: label_comp_idA,
  } = unitA.model.atomicHierarchy.atoms;
  const {
    type_symbol: type_symbolB,
    label_alt_id: label_alt_idB,
    label_atom_id: label_atom_idB,
    label_comp_id: label_comp_idB,
  } = unitB.model.atomicHierarchy.atoms;
  const { auth_seq_id: auth_seq_idA } = unitA.model.atomicHierarchy.residues;
  const { auth_seq_id: auth_seq_idB } = unitB.model.atomicHierarchy.residues;
  const { occupancy: occupancyA } = unitA.model.atomicConformation;
  const { occupancy: occupancyB } = unitB.model.atomicConformation;
  const hasOccupancy = occupancyA.isDefined && occupancyB.isDefined;

  const structConn = unitA.model === unitB.model && StructConn.Provider.get(unitA.model);
  const indexPairs = !props.forceCompute && unitA.model === unitB.model && IndexPairBonds.Provider.get(unitA.model);

  const { atomSourceIndex: sourceIndex } = unitA.model.atomicHierarchy;
  const { invertedIndex } = indexPairs ? Model.getInvertedAtomSourceIndex(unitB.model) : { invertedIndex: void 0 };

  const structConnExhaustive = unitA.model === unitB.model && StructConn.isExhaustive(unitA.model);

  // the lookup queries need to happen in the "unitB space".
  // that means _imageA = inverseOperB(operA(aI))
  const imageTransform = Mat4.mul(
    _imageTransform,
    unitB.conformation.operator.inverse,
    unitA.conformation.operator.matrix,
  );
  const isNotIdentity = !Mat4.isIdentity(imageTransform);

  const { center: bCenter, radius: bRadius } = unitB.boundary.sphere;
  const testDistanceSq = (bRadius + maxRadius) * (bRadius + maxRadius);

  builder.startUnitPair(unitA.id, unitB.id);
  const opKeyA = unitA.conformation.operator.key;
  const opKeyB = unitB.conformation.operator.key;

  const candidates = getCandidates(unitA, unitB, bRadius + maxRadius);
  for (let ci = 0, cl = candidates.length; ci < cl; ci++) {
    const _aI = candidates[ci] as StructureElement.UnitIndex;
    const aI = atomsA[_aI];
    v3set(_imageA, xA[aI], yA[aI], zA[aI]);
    if (isNotIdentity) v3transformMat4(_imageA, _imageA, imageTransform);
    if (v3squaredDistance(_imageA, bCenter) > testDistanceSq) continue;

    if (!props.forceCompute && indexPairs) {
      const { maxDistance } = indexPairs;
      const {
        offset,
        b,
        edgeProps: { order, distance, flag, key, operatorA, operatorB },
      } = indexPairs.bonds;

      const srcA = sourceIndex.value(aI);
      const aeI = getElementIdx(type_symbolA.value(aI));
      for (let i = offset[srcA], il = offset[srcA + 1]; i < il; ++i) {
        const bI = invertedIndex![b[i]];

        const _bI = SortedArray.indexOf(unitB.elements, bI) as StructureElement.UnitIndex;
        if (_bI < 0) continue;

        const opA = operatorA[i];
        const opB = operatorB[i];
        if (opA >= 0 && opB >= 0) {
          if (opA === opB) continue;
          if (opA !== opKeyA || opB !== opKeyB) continue;
        }

        const beI = getElementIdx(type_symbolA.value(bI));

        const d = distance[i];
        const dist = getDistance(unitA, aI, unitB, bI);

        let add = false;
        if (d >= 0) {
          add = equalEps(dist, d, 0.3);
        } else if (maxDistance >= 0) {
          add = dist < maxDistance;
        } else {
          const pairingThreshold = getPairingThreshold(aeI, beI, getElementThreshold(aeI), getElementThreshold(beI));
          add = dist < pairingThreshold;

          if (isHydrogen(aeI) && isHydrogen(beI)) {
            // TODO handle molecular hydrogen
            add = false;
          }
        }

        if (add) {
          builder.add(_aI, _bI, { order: order[i], flag: flag[i], key: key[i] });
        }
      }
      continue; // assume `indexPairs` supplies all bonds
    }

    const structConnEntries = props.forceCompute ? void 0 : structConn && structConn.byAtomIndex.get(aI);
    if (structConnEntries && structConnEntries.length) {
      let added = false;
      for (const se of structConnEntries) {
        const { partnerA, partnerB } = se;
        const p = partnerA.atomIndex === aI ? partnerB : partnerA;
        const _bI = SortedArray.indexOf(unitB.elements, p.atomIndex) as StructureElement.UnitIndex;
        if (_bI < 0) continue;

        // check if the bond is within MAX_RADIUS for this pair of units
        if (getDistance(unitA, aI, unitB, p.atomIndex) > maxRadius) continue;

        builder.add(_aI, _bI, { order: se.order, flag: se.flags, key: se.rowIndex });
        added = true;
      }
      // assume, for an atom, that if any inter unit bond is given
      // all are given and thus we don't need to compute any other
      if (added) continue;
    }
    if (structConnExhaustive || props.noDistance) continue;

    const occA = occupancyA.value(aI);

    const { lookup3d } = unitB;
    // only distance bonds from here: none is longer than the largest pairing threshold (perf-audit-2 05/P16)
    const { indices, count, squaredDistances } = lookup3d.find(_imageA[0], _imageA[1], _imageA[2], distanceRadius);
    if (count === 0) continue;

    const aeI = getElementIdx(type_symbolA.value(aI));
    const isHa = isHydrogen(aeI);
    const thresholdA = getElementThreshold(aeI);
    const altA = label_alt_idA.value(aI);
    const metalA = MetalsSet.has(aeI);
    const atomIdA = label_atom_idA.value(aI);
    const compIdA = label_comp_idA.value(residueIndexA[aI]);

    for (let ni = 0; ni < count; ni++) {
      const _bI = indices[ni] as StructureElement.UnitIndex;
      const bI = atomsB[_bI];

      const altB = label_alt_idB.value(bI);
      if (altA && altB && altA !== altB) continue;

      // Do not include bonds between images of the same residue with partial occupancy.
      // TODO: is this condition good enough?
      // - It works for cases like 3WQJ (label_asym_id: I) which have partial occupancy.
      // - Does NOT work for cases like 1RB8 (DC 7) with full occupancy.
      if (hasOccupancy && occupancyB.value(bI) < 1 && occA < 1) {
        if (auth_seq_idA.value(residueIndexA[aI]) === auth_seq_idB.value(residueIndexB[bI])) {
          continue;
        }
      }

      if (structConn && unitA.model === unitB.model) {
        const residuePair = sortedCantorPairing(residueIndexA[aI], residueIndexB[bI]);
        // Do not add bonds for residue pairs that have a structConn entry
        if (structConn.residueCantorPairs.has(residuePair)) continue;
      }

      const beI = getElementIdx(type_symbolB.value(bI)!);

      const isHb = isHydrogen(beI);
      if (isHa && isHb) continue;

      const isMetal = (metalA || MetalsSet.has(beI)) && !(isHa || isHb);

      const dist = Math.sqrt(squaredDistances[ni]);
      if (dist === 0) continue;

      const pairingThreshold = getPairingThreshold(aeI, beI, thresholdA, getElementThreshold(beI));
      if (dist <= pairingThreshold) {
        const atomIdB = label_atom_idB.value(bI);
        const compIdB = label_comp_idB.value(residueIndexB[bI]);
        builder.add(_aI, _bI, {
          order: getInterBondOrderFromTable(compIdA, compIdB, atomIdA, atomIdB),
          flag: (isMetal ? BondType.Flag.MetallicCoordination : BondType.Flag.Covalent) | BondType.Flag.Computed,
          key: -1,
        });
      }
    }
  }

  builder.finishUnitPair();
}

function canAddFromIndexPairBonds(structure: Structure) {
  for (const m of structure.models) {
    const indexPairs = IndexPairBonds.Provider.get(m);
    if (!indexPairs?.hasOperators) return false;
  }
  for (const u of structure.units) {
    if (u.conformation.operator.key === -1) return false;
  }
  return true;
}

function addIndexPairBonds(
  structure: Structure,
  builder: InterUnitGraph.Builder<number, StructureElement.UnitIndex, InterUnitEdgeProps>,
) {
  const opUnits = new Map<number, Set<Unit>>();
  for (const u of structure.units) {
    const { key } = u.conformation.operator;
    if (opUnits.has(key)) opUnits.get(key)!.add(u);
    else opUnits.set(key, new Set([u]));
  }

  for (const m of structure.models) {
    const indexPairs = IndexPairBonds.Provider.get(m)!;
    const { a, b } = indexPairs.bonds;
    const { order, flag, key, operatorA, operatorB } = indexPairs.bonds.edgeProps;
    const { invertedIndex } = Model.getInvertedAtomSourceIndex(m);

    const atomsToUnits = new Map<ElementIndex, Set<Unit>>();
    for (const u of structure.units) {
      if (u.model !== m) continue;

      for (let i = 0, il = u.elements.length; i < il; ++i) {
        const aI = u.elements[i];
        if (atomsToUnits.has(aI)) atomsToUnits.get(aI)!.add(u);
        else atomsToUnits.set(aI, new Set([u]));
      }
    }

    const pairs = new Map<number, Set<number>>();
    for (let i = 0, il = operatorA.length; i < il; ++i) {
      let unitsA: Set<Unit> | undefined;
      let unitsB: Set<Unit> | undefined;

      if (operatorA[i] === operatorB[i]) {
        unitsA = atomsToUnits.get(invertedIndex[a[i]]);
        unitsB = atomsToUnits.get(invertedIndex[b[i]]);
      } else {
        unitsA = opUnits.get(operatorA[i]);
        unitsB = opUnits.get(operatorB[i]);
      }
      if (!unitsA || !unitsB) continue;

      for (const uA of unitsA) {
        if (operatorA[i] !== uA.conformation.operator.key) continue;

        for (const uB of unitsB) {
          if (operatorB[i] !== uB.conformation.operator.key) continue;
          if (uA === uB || !Unit.isAtomic(uA) || !Unit.isAtomic(uB)) continue;
          if (uA.id > uB.id) continue;

          const h = cantorPairing(uA.id, uB.id);
          if (pairs.has(h)) pairs.get(h)!.add(i);
          else pairs.set(h, new Set([i]));
        }
      }
    }

    const unitIds: [number, number] = [-1, -1];
    pairs.forEach((indices, h) => {
      const [unitIdA, unitIdB] = invertCantorPairing(unitIds, h);
      const uA = structure.unitMap.get(unitIdA);
      const uB = structure.unitMap.get(unitIdB);
      builder.startUnitPair(unitIdA, unitIdB);
      indices.forEach((i) => {
        const aI = invertedIndex[a[i]];
        const _aI = SortedArray.indexOf(uA.elements, aI) as StructureElement.UnitIndex;
        if (_aI < 0) return;

        const bI = invertedIndex[b[i]];
        const _bI = SortedArray.indexOf(uB.elements, bI) as StructureElement.UnitIndex;
        if (_bI < 0) return;

        builder.add(_aI, _bI, { order: order[i], flag: flag[i], key: key[i] });
      });
      builder.finishUnitPair();
    });
  }
}

export interface InterBondComputationProps extends BondComputationProps {
  validUnit: (unit: Unit) => boolean;
  validUnitPair: (structure: Structure, unitA: Unit, unitB: Unit) => boolean;
  ignoreWater: boolean;
  ignoreIon: boolean;
}

const DefaultInterBondComputationProps = {
  ...DefaultBondComputationProps,
  ignoreWater: true,
  ignoreIon: true,
};

function findBonds(structure: Structure, props: InterBondComputationProps) {
  const builder = new InterUnitGraph.Builder<number, StructureElement.UnitIndex, InterUnitEdgeProps>();
  const hasIndexPairBonds = structure.models.some((m) => IndexPairBonds.Provider.get(m));
  const hasExhaustiveStructConn = structure.models.some((m) => StructConn.isExhaustive(m));

  if (props.noCompute || (structure.isCoarseGrained && !hasIndexPairBonds && !hasExhaustiveStructConn)) {
    return new InterUnitBonds(builder.getMap());
  }

  if (!props.forceCompute && canAddFromIndexPairBonds(structure)) {
    addIndexPairBonds(structure, builder);
    return new InterUnitBonds(builder.getMap());
  }

  Structure.eachUnitPair(
    structure,
    (unitA: Unit, unitB: Unit) => {
      findPairBonds(unitA as Unit.Atomic, unitB as Unit.Atomic, props, builder);
    },
    {
      maxRadius: props.maxRadius,
      validUnit: (unit: Unit) => props.validUnit(unit),
      validUnitPair: (unitA: Unit, unitB: Unit) => props.validUnitPair(structure, unitA, unitB),
    },
  );

  return new InterUnitBonds(builder.getMap());
}

function computeInterUnitBonds(structure: Structure, props?: Partial<InterBondComputationProps>): InterUnitBonds {
  const p = { ...DefaultInterBondComputationProps, ...props };
  return findBonds(structure, {
    ...p,
    validUnit: (props && props.validUnit) || ((u) => Unit.isAtomic(u)),
    validUnitPair:
      (props && props.validUnitPair) ||
      ((s, a, b) => {
        const isValidPair = Structure.validUnitPair(s, a, b);
        if (!isValidPair) return false;

        const mtA = a.model.atomicHierarchy.derived.residue.moleculeType;
        const mtB = b.model.atomicHierarchy.derived.residue.moleculeType;
        const notWater =
          (!Unit.isAtomic(a) || mtA[a.residueIndex[a.elements[0]]] !== MoleculeType.Water) &&
          (!Unit.isAtomic(b) || mtB[b.residueIndex[b.elements[0]]] !== MoleculeType.Water);

        const notIonA = !Unit.isAtomic(a) || mtA[a.residueIndex[a.elements[0]]] !== MoleculeType.Ion;
        const notIonB = !Unit.isAtomic(b) || mtB[b.residueIndex[b.elements[0]]] !== MoleculeType.Ion;
        const notIon = notIonA && notIonB;

        const check = (notWater || !p.ignoreWater) && (notIon || !p.ignoreIon);
        if (!check) {
          // In case both units have a struct conn record, ignore other criteria
          return hasCommonStructConnRecord(a, b);
        }
        return true;
      }),
  });
}

/**
 * The struct_conn partners of a unit's atoms (sorted), per struct_conn object and unit elements: most units
 * have none, and a unit pairs with many units (its neighbours, symmetry copies, the children that reuse it), so the
 * rows are read once per unit, not per pair (perf-audit-2 05/P17)
 */
const structConnPartners = new WeakMap<StructConn, WeakMap<SortedArray<ElementIndex>, ArrayLike<ElementIndex>>>();
const NoPartners: ArrayLike<ElementIndex> = [];

function getStructConnPartners(structConn: StructConn, elements: SortedArray<ElementIndex>) {
  let byElements = structConnPartners.get(structConn);
  if (!byElements) structConnPartners.set(structConn, (byElements = new WeakMap()));
  let partners = byElements.get(elements);
  if (partners) return partners;

  const found: number[] = [];
  for (let i = 0, _i = elements.length; i < _i; i++) {
    const aI = elements[i];
    const entries = structConn.byAtomIndex.get(aI);
    if (!entries?.length) continue;
    for (const e of entries) {
      found.push(e.partnerA.atomIndex === aI ? e.partnerB.atomIndex : e.partnerA.atomIndex);
    }
  }
  partners = found.length === 0 ? NoPartners : SortedArray.ofUnsortedArray<ElementIndex>(found as ElementIndex[]);
  byElements.set(elements, partners);
  return partners;
}

/** Whether a struct_conn record links an atom of `unitA` and one of `unitB`: the smaller unit's partners in the other */
function hasCommonStructConnRecord(unitA: Unit, unitB: Unit) {
  if (unitA.model !== unitB.model || !Unit.isAtomic(unitA) || !Unit.isAtomic(unitB)) return false;
  const structConn = StructConn.Provider.get(unitA.model);
  if (!structConn) return false;

  const smaller = unitA.elements.length < unitB.elements.length ? unitA : unitB;
  const { elements } = smaller === unitA ? unitB : unitA;
  const partners = getStructConnPartners(structConn, smaller.elements);
  for (let i = 0, _i = partners.length; i < _i; i++) {
    if (SortedArray.has(elements, partners[i])) return true;
  }
  return false;
}

export { computeInterUnitBonds };
