/**
 * Copyright (c) 2017-2026 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Paul Pillot <paul.pillot@tandemai.com>
 */

import { BondType } from '../../../model/types.js';
import { type IntraUnitBondProps, IntraUnitBonds } from './data.js';
import type { Unit } from '../../unit.js';
import { IntAdjacencyGraph } from '@molstar/core/math/graph';
import {
  type BondComputationProps,
  getElementIdx,
  MetalsSet,
  getElementThreshold,
  isHydrogen,
  DefaultBondComputationProps,
  getPairingThreshold,
} from './common.js';
import { SortedArray } from '@molstar/core/data/int';
import { getIntraBondOrderFromTable } from '../../../model/properties/atomic/bonds.js';
import type { StructureElement } from '../../element.js';
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { ComponentBond } from '@molstar/model/formats/structure/property/bonds/chem_comp';
import { StructConn } from '@molstar/model/formats/structure/property/bonds/struct_conn';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { ElementIndex } from '../../../model/indexing.js';
import { equalEps } from '@molstar/core/math/linear-algebra/3d/common';
import { Model } from '../../../model/model.js';
import { sortedCantorPairing } from '@molstar/core/data/util';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3distance = Vec3.distance;

const CoarseGrainedBondMaxRadius = 6;
const CoarseGrainedIntraResidueBondMaxDistance = 5.5;
const CoarseGrainedInterResidueBondMaxDistance = 3.9;

function getGraph(
  atomA: StructureElement.UnitIndex[],
  atomB: StructureElement.UnitIndex[],
  _order: number[],
  _flags: number[],
  _key: number[],
  atomCount: number,
  props?: IntraUnitBondProps,
): IntraUnitBonds {
  const builder = new IntAdjacencyGraph.EdgeBuilder(atomCount, atomA, atomB);
  const flags = new Uint16Array(builder.slotCount);
  const order = new Int8Array(builder.slotCount);
  const key = new Uint32Array(builder.slotCount);
  for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
    builder.addNextEdge();
    builder.assignProperty(flags, _flags[i]);
    builder.assignProperty(order, _order[i]);
    builder.assignProperty(key, _key[i]);
  }

  return builder.createGraph({ flags, order, key }, props);
}

const tmpDistVecA = Vec3();
const tmpDistVecB = Vec3();
function getDistance(unit: Unit.Atomic, indexA: ElementIndex, indexB: ElementIndex) {
  unit.conformation.position(indexA, tmpDistVecA);
  unit.conformation.position(indexB, tmpDistVecB);
  return v3distance(tmpDistVecA, tmpDistVecB);
}

const __structConnAdded = new Set<StructureElement.UnitIndex>();

function findIndexPairBonds(unit: Unit.Atomic) {
  const indexPairs = IndexPairBonds.Provider.get(unit.model)!;
  const { elements: atoms } = unit;
  const { type_symbol } = unit.model.atomicHierarchy.atoms;
  const atomCount = unit.elements.length;
  const { maxDistance } = indexPairs;
  const {
    offset,
    b,
    edgeProps: { order, distance, flag, key, operatorA, operatorB },
  } = indexPairs.bonds;

  const { atomSourceIndex: sourceIndex } = unit.model.atomicHierarchy;
  const { invertedIndex } = Model.getInvertedAtomSourceIndex(unit.model);

  const atomA: StructureElement.UnitIndex[] = [];
  const atomB: StructureElement.UnitIndex[] = [];
  const flags: number[] = [];
  const orders: number[] = [];
  const keys: number[] = [];

  const opKey = unit.conformation.operator.key;

  for (let _aI = 0 as StructureElement.UnitIndex; _aI < atomCount; _aI++) {
    const aI = atoms[_aI];
    const aeI = getElementIdx(type_symbol.value(aI));
    const isHa = isHydrogen(aeI);

    const srcA = sourceIndex.value(aI);

    for (let i = offset[srcA], il = offset[srcA + 1]; i < il; ++i) {
      const bI = invertedIndex[b[i]];
      if (aI >= bI) continue;

      const _bI = SortedArray.indexOf(unit.elements, bI) as StructureElement.UnitIndex;
      if (_bI < 0) continue;

      const opA = operatorA[i];
      const opB = operatorB[i];
      if ((opA >= 0 && opA !== opKey) || (opB >= 0 && opB !== opKey)) continue;

      const beI = getElementIdx(type_symbol.value(bI));

      const d = distance[i];
      const dist = getDistance(unit, aI, bI);

      let add = false;
      if (d >= 0) {
        add = equalEps(dist, d, 0.3);
      } else if (maxDistance >= 0) {
        add = dist < maxDistance;
      } else {
        const pairingThreshold = getPairingThreshold(aeI, beI, getElementThreshold(aeI), getElementThreshold(beI));
        add = dist < pairingThreshold;

        if (isHa && isHydrogen(beI)) {
          // TODO handle molecular hydrogen
          add = false;
        }
      }

      if (add) {
        atomA[atomA.length] = _aI;
        atomB[atomB.length] = _bI;
        orders[orders.length] = order[i];
        flags[flags.length] = flag[i];
        keys[keys.length] = key[i];
      }
    }
  }

  return getGraph(atomA, atomB, orders, flags, keys, atomCount, {
    canRemap: false,
    cacheable: indexPairs.cacheable,
  });
}

/**
 * Above every distance threshold (`getPairingThreshold` is at most (2.88 + 2.88) / 1.95 ≈ 2.95 Å): distance bonds of
 * atomic models need no larger neighbour query. Dictionary bonds keep `maxRadius`, through the residue's atoms
 * (perf-audit 05/P7).
 */
const DistanceBondMaxRadius = 3.1;

/** Above this many atoms a residue's dictionary pairs come from the grid: walking the residue per atom is quadratic */
const LargeResidueAtomCount = 256;

const structConnResidueCache = new WeakMap<StructConn, Set<number>>();
/** Residues of the struct_conn entries: no pair of residues without one is in `residueCantorPairs` */
function getStructConnResidues(structConn: StructConn) {
  let residues = structConnResidueCache.get(structConn);
  if (!residues) {
    residues = new Set();
    for (const e of structConn.entries) {
      residues.add(e.partnerA.residueIndex);
      residues.add(e.partnerB.residueIndex);
    }
    structConnResidueCache.set(structConn, residues);
  }
  return residues;
}

function findBonds(unit: Unit.Atomic, props: BondComputationProps): IntraUnitBonds {
  const isCoarseGrained = Model.isCoarseGrained(unit.model);
  const maxRadius = isCoarseGrained ? CoarseGrainedBondMaxRadius : props.maxRadius;
  // the grid's inclusion test is `distSq <= r * r`; dictionary pairs found through the residue use the same
  const maxRadiusSq = maxRadius * maxRadius;
  const queryRadius = isCoarseGrained ? maxRadius : Math.min(maxRadius, DistanceBondMaxRadius);

  const { x, y, z } = unit.model.atomicConformation;
  const atomCount = unit.elements.length;
  const { elements: atoms, residueIndex, chainIndex } = unit;
  const { type_symbol, label_atom_id, label_alt_id, label_comp_id } = unit.model.atomicHierarchy.atoms;
  const { label_seq_id } = unit.model.atomicHierarchy.residues;
  const { traceElementIndex } = unit.model.atomicHierarchy.derived.residue;
  const { offsets: residueOffsets } = unit.model.atomicHierarchy.residueAtomSegments;
  const { index } = unit.model.atomicHierarchy;
  const { byEntityKey } = unit.model.sequence;
  const query3d = unit.lookup3d;

  const structConn = StructConn.Provider.get(unit.model);
  const component = ComponentBond.Provider.get(unit.model);

  const structConnExhaustive = StructConn.isExhaustive(unit.model);

  const atomA: StructureElement.UnitIndex[] = [];
  const atomB: StructureElement.UnitIndex[] = [];
  const flags: number[] = [];
  const order: number[] = [];
  const key: number[] = [];

  let lastResidue = -1;
  let componentEntry: ComponentBond.Entry | undefined = void 0;

  let isWatery = true,
    isDictionaryBased = true,
    isSequenced = true;

  const structConnAdded = __structConnAdded;
  const hasStructConnEntries = !!structConn?.residueCantorPairs.size;
  const structConnResidues = hasStructConnEntries ? getStructConnResidues(structConn) : undefined;

  // element index per atom, looked up once instead of per neighbour (`getElementIdx` is in [-1, 118])
  const elementIdx = new Int8Array(atomCount);
  for (let _aI = 0; _aI < atomCount; _aI++) {
    const elem = type_symbol.value(atoms[_aI]);
    if (isWatery && elem !== 'H' && elem !== 'O') isWatery = false;
    elementIdx[_aI] = getElementIdx(elem);
  }

  /** A bond from distance (`dist` > 0) between atoms not both hydrogen */
  const addDistanceBond = (
    _aI: StructureElement.UnitIndex,
    _bI: StructureElement.UnitIndex,
    aI: ElementIndex,
    bI: ElementIndex,
    raI: number,
    rbI: number,
    dist: number,
    thresholdA: number,
    isMetal: boolean,
    compId: string,
    atomIdA: string,
    seqIdA: number,
  ) => {
    let flag = false;
    if (isCoarseGrained) {
      if (raI === rbI) {
        // intra residue bonds
        flag = dist <= CoarseGrainedIntraResidueBondMaxDistance;
      } else {
        // inter residue "backbone" bonds
        flag =
          dist <= CoarseGrainedInterResidueBondMaxDistance &&
          traceElementIndex[raI] === aI &&
          traceElementIndex[rbI] === bI;
      }
    } else {
      const beI = elementIdx[_bI];
      const pairingThreshold = getPairingThreshold(elementIdx[_aI], beI, thresholdA, getElementThreshold(beI));
      flag = dist <= pairingThreshold;
    }

    if (flag) {
      atomA[atomA.length] = _aI;
      atomB[atomB.length] = _bI;
      order[order.length] = getIntraBondOrderFromTable(compId, atomIdA, label_atom_id.value(bI));
      flags[flags.length] =
        (isMetal ? BondType.Flag.MetallicCoordination : BondType.Flag.Covalent) | BondType.Flag.Computed;
      key[key.length] = -1;

      const seqIdB = label_seq_id.value(rbI);

      if (seqIdA === seqIdB) isDictionaryBased = false;
      if (Math.abs(seqIdA - seqIdB) > 1) isSequenced = false;
    }
  };

  for (let _aI = 0 as StructureElement.UnitIndex; _aI < atomCount; _aI++) {
    const aI = atoms[_aI];

    const structConnEntries = props.forceCompute ? void 0 : structConn && structConn.byAtomIndex.get(aI);
    let hasStructConn = false;
    if (structConnEntries) {
      for (const se of structConnEntries) {
        const { partnerA, partnerB } = se;
        // symmetry must be the same for intra-unit bonds
        if (partnerA.symmetry !== partnerB.symmetry) continue;

        const p = partnerA.atomIndex === aI ? partnerB : partnerA;
        const _bI = SortedArray.indexOf(unit.elements, p.atomIndex) as StructureElement.UnitIndex;
        if (_bI < 0 || atoms[_bI] < aI) continue;

        atomA[atomA.length] = _aI;
        atomB[atomB.length] = _bI;
        flags[flags.length] = se.flags;
        order[order.length] = se.order;
        key[key.length] = se.rowIndex;

        if (!hasStructConn) structConnAdded.clear();
        hasStructConn = true;
        structConnAdded.add(_bI);
      }
    }
    if (structConnExhaustive) continue;

    const raI = residueIndex[aI];
    const seqIdA = label_seq_id.value(raI);
    const compId = label_comp_id.value(aI);

    if (!props.forceCompute && raI !== lastResidue) {
      if (!!component && component.entries.has(compId)) {
        const entitySeq = byEntityKey[index.getEntityFromChain(chainIndex[aI])];
        if (entitySeq && entitySeq.sequence.microHet.has(seqIdA)) {
          // compute for sequence positions with micro-heterogeneity
          componentEntry = void 0;
        } else {
          componentEntry = component.entries.get(compId)!;
        }
      } else {
        componentEntry = void 0;
      }
    }
    lastResidue = raI;

    const aeI = elementIdx[_aI];
    const isHa = isHydrogen(aeI);
    const atomIdA = label_atom_id.value(aI);
    const componentPairs = componentEntry ? componentEntry.get(atomIdA, isHa) : void 0;

    const thresholdA = getElementThreshold(aeI);
    const altA = label_alt_id.value(aI);
    const metalA = MetalsSet.has(aeI);
    // only residues with a struct_conn entry can be in a pair that has one
    const checkResiduePairs = hasStructConnEntries && structConnResidues!.has(raI);

    // "component dictionary" bonds, within `maxRadius`: the later atoms of the residue (it is contiguous in the unit);
    // in a residue of more than `LargeResidueAtomCount` atoms, those the grid finds within `maxRadius`
    if (componentPairs) {
      const skipResidue = checkResiduePairs && structConn!.residueCantorPairs.has(sortedCantorPairing(raI, raI));
      const xA = x[aI],
        yA = y[aI],
        zA = z[aI];
      const near =
        residueOffsets[raI + 1] - residueOffsets[raI] > LargeResidueAtomCount
          ? query3d.find(xA, yA, zA, maxRadius)
          : undefined;
      const candidateCount = near ? near.count : atomCount - _aI - 1;
      for (let k = 0; k < candidateCount; k++) {
        let _bI: StructureElement.UnitIndex, bI: ElementIndex, distSq: number;
        if (near) {
          _bI = near.indices[k];
          bI = atoms[_bI];
          if (bI <= aI || residueIndex[bI] !== raI) continue;
          distSq = near.squaredDistances[k];
        } else {
          _bI = (_aI + 1 + k) as StructureElement.UnitIndex;
          bI = atoms[_bI];
          if (residueIndex[bI] !== raI) break;
          const dx = x[bI] - xA;
          const dy = y[bI] - yA;
          const dz = z[bI] - zA;
          distSq = dx * dx + dy * dy + dz * dz;
          // as the grid's test: NaN coordinates are never within
          if (!(distSq <= maxRadiusSq)) continue;
        }
        if (hasStructConn && structConnAdded.has(_bI)) continue;

        const altB = label_alt_id.value(bI);
        if (altA && altB && altA !== altB) continue;

        // Do not add bonds for residues that have a structConn entry
        if (skipResidue) continue;

        const beI = elementIdx[_bI];
        const isHb = isHydrogen(beI);
        if (isHa && isHb) continue;

        const isMetal = (metalA || MetalsSet.has(beI)) && !(isHa || isHb);

        const e = componentPairs.get(label_atom_id.value(bI), isHb);
        if (e) {
          atomA[atomA.length] = _aI;
          atomB[atomB.length] = _bI;
          order[order.length] = e.order;
          let flag = e.flags;
          if (isMetal) {
            flag &= ~BondType.Flag.Covalent;
            flag |= BondType.Flag.MetallicCoordination;
          }
          flags[flags.length] = flag;
          key[key.length] = e.key;
          continue;
        }
        if (!isHa && !isHb) continue; // If Hydrogen not in the component dictionary, fallback to distance-based.
        if (props.noDistance) continue;

        const dist = Math.sqrt(distSq);
        if (dist === 0) continue;
        addDistanceBond(_aI, _bI, aI, bI, raI, raI, dist, thresholdA, isMetal, compId, atomIdA, seqIdA);
      }
    }

    if (props.noDistance) continue;

    const { indices, count, squaredDistances } = query3d.find(x[aI], y[aI], z[aI], queryRadius);

    for (let ni = 0; ni < count; ni++) {
      const _bI = indices[ni];
      if (hasStructConn && structConnAdded.has(_bI)) continue;

      const bI = atoms[_bI];
      if (bI <= aI) continue;

      const rbI = residueIndex[bI];
      // the residue's own atoms went through the dictionary above
      if (componentPairs && raI === rbI) continue;

      const altB = label_alt_id.value(bI);
      if (altA && altB && altA !== altB) continue;

      // Do not add bonds for residues that have a structConn entry
      if (checkResiduePairs) {
        const residuePair = sortedCantorPairing(raI, rbI);
        if (structConn!.residueCantorPairs.has(residuePair)) continue;
      }

      const beI = elementIdx[_bI];

      const isHb = isHydrogen(beI);
      if (isHa && isHb) continue;

      const isMetal = (metalA || MetalsSet.has(beI)) && !(isHa || isHb);

      const dist = Math.sqrt(squaredDistances[ni]);
      if (dist === 0) continue;
      addDistanceBond(_aI, _bI, aI, bI, raI, rbI, dist, thresholdA, isMetal, compId, atomIdA, seqIdA);
    }
  }

  return getGraph(atomA, atomB, order, flags, key, atomCount, {
    canRemap: isWatery || (isDictionaryBased && isSequenced),
  });
}

function canGetFromIndexPairBonds(unit: Unit.Atomic) {
  if (unit.conformation.operator.key === -1) return false;

  const indexPairs = IndexPairBonds.Provider.get(unit.model);
  return !!indexPairs?.hasOperators;
}

function getIndexPairBonds(unit: Unit.Atomic) {
  const indexPairs = IndexPairBonds.Provider.get(unit.model)!;
  const bonds = indexPairs.bySameOperator.get(unit.conformation.operator.key);
  if (!bonds) return IntraUnitBonds.Empty;

  const {
    a,
    b,
    edgeProps: { key, flag, order },
  } = indexPairs.bonds;
  const { invertedIndex } = Model.getInvertedAtomSourceIndex(unit.model);
  const { elements } = unit;

  const atomA: StructureElement.UnitIndex[] = [];
  const atomB: StructureElement.UnitIndex[] = [];
  const flags: number[] = [];
  const orders: number[] = [];
  const keys: number[] = [];

  for (let j = 0, jl = bonds.length; j < jl; ++j) {
    const i = bonds[j];
    if (a[i] >= b[i]) continue;

    const aI = invertedIndex[a[i]];
    const _aI = SortedArray.indexOf(elements, aI) as StructureElement.UnitIndex;
    if (_aI < 0) continue;

    const bI = invertedIndex[b[i]];
    const _bI = SortedArray.indexOf(elements, bI) as StructureElement.UnitIndex;
    if (_bI < 0) continue;

    atomA[atomA.length] = _aI;
    atomB[atomB.length] = _bI;
    flags[flags.length] = flag[i];
    orders[orders.length] = order[i];
    keys[keys.length] = key[i];
  }

  return getGraph(atomA, atomB, orders, flags, keys, elements.length, {
    canRemap: false,
    cacheable: indexPairs.cacheable,
  });
}

function computeIntraUnitBonds(unit: Unit.Atomic, props?: Partial<BondComputationProps>) {
  const p = { ...DefaultBondComputationProps, ...props };
  if (p.noCompute) return IntraUnitBonds.Empty;
  if (unit.elements.length <= 1) return IntraUnitBonds.Empty;

  if (!p.forceCompute && IndexPairBonds.Provider.get(unit.model)) {
    return canGetFromIndexPairBonds(unit) ? getIndexPairBonds(unit) : findIndexPairBonds(unit);
  } else {
    return findBonds(unit, p);
  }
}

export { computeIntraUnitBonds };
