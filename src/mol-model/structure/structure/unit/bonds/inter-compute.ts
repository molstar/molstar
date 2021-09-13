/**
 * Copyright (c) 2017-2021 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType, MoleculeType } from '../../../model/types';
import { Structure } from '../../structure';
import { Unit } from '../../unit';
import { getElementIdx, getElementPairThreshold, getElementThreshold, isHydrogen, BondComputationProps, MetalsSet, DefaultBondComputationProps } from './common';
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

const tmpDistVecA = Vec3();
const tmpDistVecB = Vec3();
function getDistance(unitA: Unit.Atomic, indexA: ElementIndex, unitB: Unit.Atomic, indexB: ElementIndex) {
    unitA.conformation.position(indexA, tmpDistVecA);
    unitB.conformation.position(indexB, tmpDistVecB);
    return Vec3.distance(tmpDistVecA, tmpDistVecB);
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

    for (let _aI = 0 as StructureElement.UnitIndex; _aI < atomCount; _aI++) {
        const aI = atomsA[_aI];
        Vec3.set(_imageA, xA[aI], yA[aI], zA[aI]);
        if (isNotIdentity) Vec3.transformMat4(_imageA, _imageA, imageTransform);
        if (Vec3.squaredDistance(_imageA, bCenter) > testDistanceSq) continue;

        if (!props.forceCompute && indexPairs) {
            const { maxDistance } = indexPairs;
            const { offset, b, edgeProps: { order, distance, flag } } = indexPairs.bonds;

            const srcA = sourceIndex.value(aI);
            for (let i = offset[srcA], il = offset[srcA + 1]; i < il; ++i) {
                const bI = invertedIndex![b[i]];

                const _bI = SortedArray.indexOf(unitB.elements, bI) as StructureElement.UnitIndex;
                if (_bI < 0) continue;
                if (type_symbolA.value(aI) === 'H' && type_symbolB.value(bI) === 'H') continue;

                const d = distance[i];
                const dist = getDistance(unitA, aI, unitB, bI);
                if ((d !== -1 && equalEps(dist, d, 0.5)) || dist < maxDistance) {
                    builder.add(_aI, _bI, { order: order[i], flag: flag[i] });
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

                builder.add(_aI, _bI, { order: se.order, flag: se.flags });
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
                if (auth_seq_idA.value(aI) === auth_seq_idB.value(bI)) {
                    continue;
                }
            }

            const beI = getElementIdx(type_symbolB.value(bI)!);

            const isHb = isHydrogen(beI);
            if (isHa && isHb) continue;

            const isMetal = (metalA || MetalsSet.has(beI)) && !(isHa || isHb);

            const dist = Math.sqrt(squaredDistances[ni]);
            if (dist === 0) continue;

            const thresholdAB = getElementPairThreshold(aeI, beI);
            const pairingThreshold = thresholdAB > 0
                ? thresholdAB
                : beI < 0
                    ? thresholdA
                    : (thresholdA + getElementThreshold(beI)) / 1.95; // not sure if avg or min but max is too big

            if (dist <= pairingThreshold) {
                const atomIdB = label_atom_idB.value(bI);
                const compIdB = label_comp_idB.value(residueIndexB[bI]);
                builder.add(_aI, _bI, {
                    order: getInterBondOrderFromTable(compIdA, compIdB, atomIdA, atomIdB),
                    flag: (isMetal ? BondType.Flag.MetallicCoordination : BondType.Flag.Covalent) | BondType.Flag.Computed
                });
            }
        }
    }

    builder.finishUnitPair();
}

export interface InterBondComputationProps extends BondComputationProps {
    validUnitPair: (structure: Structure, unitA: Unit, unitB: Unit) => boolean
    ignoreWater: boolean
}

const DefaultInterBondComputationProps = {
    ...DefaultBondComputationProps,
    ignoreWater: true
};

function findBonds(structure: Structure, props: InterBondComputationProps) {
    const builder = new InterUnitGraph.Builder<number, StructureElement.UnitIndex, InterUnitEdgeProps>();
    const hasIndexPairBonds = structure.models.some(m => IndexPairBonds.Provider.get(m));

    if (props.noCompute || (structure.isCoarseGrained && !hasIndexPairBonds)) {
        // TODO add function that only adds bonds defined in structConn and avoids using
        //      structure.lookup and unit.lookup (expensive for large structure and not
        //      needed for archival files or files with an MD topology)
        return new InterUnitBonds(builder.getMap());
    }

    Structure.eachUnitPair(structure, (unitA: Unit, unitB: Unit) => {
        findPairBonds(unitA as Unit.Atomic, unitB as Unit.Atomic, props, builder);
    }, {
        maxRadius: props.maxRadius,
        validUnit: (unit: Unit) => Unit.isAtomic(unit),
        validUnitPair: (unitA: Unit, unitB: Unit) => props.validUnitPair(structure, unitA, unitB)
    });

    return new InterUnitBonds(builder.getMap());
}

function computeInterUnitBonds(structure: Structure, props?: Partial<InterBondComputationProps>): InterUnitBonds {
    const p = { ...DefaultInterBondComputationProps, ...props };
    return findBonds(structure, {
        ...p,
        validUnitPair: (props && props.validUnitPair) || ((s, a, b) => {
            const mtA = a.model.atomicHierarchy.derived.residue.moleculeType;
            const mtB = b.model.atomicHierarchy.derived.residue.moleculeType;
            const notWater = (
                (!Unit.isAtomic(a) || mtA[a.residueIndex[a.elements[0]]] !== MoleculeType.Water) &&
                (!Unit.isAtomic(b) || mtB[b.residueIndex[b.elements[0]]] !== MoleculeType.Water)
            );
            return Structure.validUnitPair(s, a, b) && (notWater || !p.ignoreWater);
        }),
    });
}

export { computeInterUnitBonds };
