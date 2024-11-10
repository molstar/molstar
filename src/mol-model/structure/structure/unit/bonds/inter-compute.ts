/**
 * Copyright (c) 2017-2024 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType, MoleculeType } from '../../../model/types';
import { Structure } from '../../structure';
import { Unit } from '../../unit';
import { getElementIdx, getElementThreshold, isHydrogen, BondComputationProps, MetalsSet, DefaultBondComputationProps, getPairingThreshold } from './common';
import { InterUnitBonds, InterUnitEdgeProps } from './data';
import { SortedArray } from '../../../../../mol-data/int';
import { Vec3, Mat4 } from '../../../../../mol-math/linear-algebra';
import { StructureElement } from '../../element';
import { ElementIndex } from '../../../model/indexing';
import { getInterBondOrderFromTable } from '../../../model/properties/atomic/bonds';
import { IndexPairBonds } from '../../../../../mol-model-formats/structure/property/bonds/index-pair';
import { InterUnitGraph } from '../../../../../mol-math/graph/inter-unit-graph';
import { StructConn } from '../../../../../mol-model-formats/structure/property/bonds/struct_conn';
import { equalEps } from '../../../../../mol-math/linear-algebra/3d/common';
import { Model } from '../../../model';
import { cantorPairing, invertCantorPairing } from '../../../../../mol-data/util';

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

function findPairBonds(unitA: Unit.Atomic, unitB: Unit.Atomic, props: BondComputationProps, builder: InterUnitGraph.Builder<number, StructureElement.UnitIndex, InterUnitEdgeProps>) {
    const { maxRadius } = props;

    const { elements: atomsA, residueIndex: residueIndexA } = unitA;
    const { x: xA, y: yA, z: zA } = unitA.model.atomicConformation;
    const { elements: atomsB, residueIndex: residueIndexB } = unitB;
    const atomCount = unitA.elements.length;

    const { type_symbol: type_symbolA, label_alt_id: label_alt_idA, label_atom_id: label_atom_idA, label_comp_id: label_comp_idA } = unitA.model.atomicHierarchy.atoms;
    const { type_symbol: type_symbolB, label_alt_id: label_alt_idB, label_atom_id: label_atom_idB, label_comp_id: label_comp_idB } = unitB.model.atomicHierarchy.atoms;
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
    const imageTransform = Mat4.mul(_imageTransform, unitB.conformation.operator.inverse, unitA.conformation.operator.matrix);
    const isNotIdentity = !Mat4.isIdentity(imageTransform);

    const { center: bCenter, radius: bRadius } = unitB.boundary.sphere;
    const testDistanceSq = (bRadius + maxRadius) * (bRadius + maxRadius);

    builder.startUnitPair(unitA.id, unitB.id);
    const opKeyA = unitA.conformation.operator.key;
    const opKeyB = unitB.conformation.operator.key;

    for (let _aI = 0 as StructureElement.UnitIndex; _aI < atomCount; _aI++) {
        const aI = atomsA[_aI];
        v3set(_imageA, xA[aI], yA[aI], zA[aI]);
        if (isNotIdentity) v3transformMat4(_imageA, _imageA, imageTransform);
        if (v3squaredDistance(_imageA, bCenter) > testDistanceSq) continue;

        if (!props.forceCompute && indexPairs) {
            const { maxDistance } = indexPairs;
            const { offset, b, edgeProps: { order, distance, flag, key, operatorA, operatorB } } = indexPairs.bonds;

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
                    const pairingThreshold = getPairingThreshold(
                        aeI, beI, getElementThreshold(aeI), getElementThreshold(beI)
                    );
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
        if (structConnExhaustive) continue;

        const occA = occupancyA.value(aI);

        const { lookup3d } = unitB;
        const { indices, count, squaredDistances } = lookup3d.find(_imageA[0], _imageA[1], _imageA[2], maxRadius);
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
                    key: -1
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

function addIndexPairBonds(structure: Structure, builder: InterUnitGraph.Builder<number, StructureElement.UnitIndex, InterUnitEdgeProps>) {
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

        const pairs = new Map<number, Set<number>>();
        for (let i = 0, il = operatorA.length; i < il; ++i) {
            const unitsA = opUnits.get(operatorA[i]);
            const unitsB = opUnits.get(operatorB[i]);
            if (!unitsA || !unitsB) continue;

            for (const uA of unitsA) {
                for (const uB of unitsB) {
                    if (uA === uB || !Unit.isAtomic(uA) || !Unit.isAtomic(uB)) continue;
                    if (uA.id > uB.id) continue;

                    const h = cantorPairing(uA.id, uB.id);
                    if (pairs.has(h)) pairs.get(h)!.add(i);
                    else pairs.set(h, new Set([i]));
                }
            }
        }

        const { invertedIndex } = Model.getInvertedAtomSourceIndex(m);
        const unitIds: [number, number] = [-1, -1];
        pairs.forEach((indices, h) => {
            const [unitIdA, unitIdB] = invertCantorPairing(unitIds, h);
            const uA = structure.unitMap.get(unitIdA);
            const uB = structure.unitMap.get(unitIdB);
            builder.startUnitPair(unitIdA, unitIdB);
            indices.forEach(i => {
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
    validUnit: (unit: Unit) => boolean
    validUnitPair: (structure: Structure, unitA: Unit, unitB: Unit) => boolean
    ignoreWater: boolean
    ignoreIon: boolean
}

const DefaultInterBondComputationProps = {
    ...DefaultBondComputationProps,
    ignoreWater: true,
    ignoreIon: true,
};

function findBonds(structure: Structure, props: InterBondComputationProps) {
    const builder = new InterUnitGraph.Builder<number, StructureElement.UnitIndex, InterUnitEdgeProps>();
    const hasIndexPairBonds = structure.models.some(m => IndexPairBonds.Provider.get(m));
    const hasExhaustiveStructConn = structure.models.some(m => StructConn.isExhaustive(m));

    if (props.noCompute || (structure.isCoarseGrained && !hasIndexPairBonds && !hasExhaustiveStructConn)) {
        return new InterUnitBonds(builder.getMap());
    }

    if (!props.forceCompute && canAddFromIndexPairBonds(structure)) {
        addIndexPairBonds(structure, builder);
        return new InterUnitBonds(builder.getMap());
    }

    Structure.eachUnitPair(structure, (unitA: Unit, unitB: Unit) => {
        findPairBonds(unitA as Unit.Atomic, unitB as Unit.Atomic, props, builder);
    }, {
        maxRadius: props.maxRadius,
        validUnit: (unit: Unit) => props.validUnit(unit),
        validUnitPair: (unitA: Unit, unitB: Unit) => props.validUnitPair(structure, unitA, unitB)
    });

    return new InterUnitBonds(builder.getMap());
}

function computeInterUnitBonds(structure: Structure, props?: Partial<InterBondComputationProps>): InterUnitBonds {
    const p = { ...DefaultInterBondComputationProps, ...props };
    return findBonds(structure, {
        ...p,
        validUnit: (props && props.validUnit) || (u => Unit.isAtomic(u)),
        validUnitPair: (props && props.validUnitPair) || ((s, a, b) => {
            const isValidPair = Structure.validUnitPair(s, a, b);
            if (!isValidPair) return false;

            const mtA = a.model.atomicHierarchy.derived.residue.moleculeType;
            const mtB = b.model.atomicHierarchy.derived.residue.moleculeType;
            const notWater = (
                (!Unit.isAtomic(a) || mtA[a.residueIndex[a.elements[0]]] !== MoleculeType.Water) &&
                (!Unit.isAtomic(b) || mtB[b.residueIndex[b.elements[0]]] !== MoleculeType.Water)
            );

            const notIonA = (!Unit.isAtomic(a) || mtA[a.residueIndex[a.elements[0]]] !== MoleculeType.Ion);
            const notIonB = (!Unit.isAtomic(b) || mtB[b.residueIndex[b.elements[0]]] !== MoleculeType.Ion);
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

function hasCommonStructConnRecord(unitA: Unit, unitB: Unit) {
    if (unitA.model !== unitB.model || !Unit.isAtomic(unitA) || !Unit.isAtomic(unitB)) return false;
    const structConn = StructConn.Provider.get(unitA.model);
    if (!structConn) return false;

    const smaller = unitA.elements.length < unitB.elements.length ? unitA : unitB;
    const bigger = unitA.elements.length >= unitB.elements.length ? unitA : unitB;

    const { elements: xs } = smaller;
    const { elements: ys } = bigger;
    const { indexOf } = SortedArray;

    for (let i = 0, _i = xs.length; i < _i; i++) {
        const aI = xs[i];
        const entries = structConn.byAtomIndex.get(aI);
        if (!entries?.length) continue;

        for (const e of entries) {
            const bI = e.partnerA.atomIndex === aI ? e.partnerB.atomIndex : e.partnerA.atomIndex;
            if (indexOf(ys, bI) >= 0) return true;
        }
    }
    return false;
}

export { computeInterUnitBonds };
