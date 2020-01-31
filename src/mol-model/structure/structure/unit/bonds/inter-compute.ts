/**
 * Copyright (c) 2017-2019 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../model/types';
import Structure from '../../structure';
import Unit from '../../unit';
import { getElementIdx, getElementPairThreshold, getElementThreshold, isHydrogen, BondComputationProps, MetalsSet, DefaultBondComputationProps } from './common';
import { InterUnitBonds } from './data';
import { UniqueArray } from '../../../../../mol-data/generic';
import { SortedArray } from '../../../../../mol-data/int';
import { Vec3, Mat4 } from '../../../../../mol-math/linear-algebra';
import StructureElement from '../../element';
import { StructConn } from '../../../../../mol-model-formats/structure/mmcif/bonds';
import { ElementIndex } from '../../../model/indexing';
import { getInterBondOrderFromTable } from '../../../model/properties/atomic/bonds';
import { IndexPairBonds } from '../../../../../mol-model-formats/structure/mmcif/bonds/index-pair';

const MAX_RADIUS = 4;

function addMapEntry<A, B>(map: Map<A, B[]>, a: A, b: B) {
    if (map.has(a)) map.get(a)!.push(b);
    else map.set(a, [b]);
}

interface PairState {
    mapAB: Map<number, InterUnitBonds.BondInfo[]>,
    mapBA: Map<number, InterUnitBonds.BondInfo[]>,
    bondedA: UniqueArray<StructureElement.UnitIndex, StructureElement.UnitIndex>,
    bondedB: UniqueArray<StructureElement.UnitIndex, StructureElement.UnitIndex>
}

function addBond(indexA: StructureElement.UnitIndex, indexB: StructureElement.UnitIndex, order: number, flag: BondType.Flag, state: PairState) {
    addMapEntry(state.mapAB, indexA, { indexB, props: { order, flag } });
    addMapEntry(state.mapBA, indexB, { indexB: indexA, props: { order, flag } });

    UniqueArray.add(state.bondedA, indexA, indexA);
    UniqueArray.add(state.bondedB, indexB, indexB);
}

const tmpDistVecA = Vec3()
const tmpDistVecB = Vec3()
function getDistance(unitA: Unit.Atomic, indexA: ElementIndex, unitB: Unit.Atomic, indexB: ElementIndex) {
    unitA.conformation.position(indexA, tmpDistVecA)
    unitB.conformation.position(indexB, tmpDistVecB)
    return Vec3.distance(tmpDistVecA, tmpDistVecB)
}

const _imageTransform = Mat4.zero();

function findPairBonds(unitA: Unit.Atomic, unitB: Unit.Atomic, props: BondComputationProps, map: Map<number, InterUnitBonds.UnitPairBonds[]>) {
    const state: PairState = { mapAB: new Map(), mapBA: new Map(), bondedA: UniqueArray.create(), bondedB: UniqueArray.create() };
    let bondCount = 0;

    const { elements: atomsA, residueIndex: residueIndexA } = unitA;
    const { x: xA, y: yA, z: zA } = unitA.model.atomicConformation;
    const { elements: atomsB, residueIndex: residueIndexB } = unitB;
    const atomCount = unitA.elements.length;

    // const { type_symbol: type_symbolA, label_alt_id: label_alt_idA } = unitA.model.atomicHierarchy.atoms;
    // const { type_symbol: type_symbolB, label_alt_id: label_alt_idB } = unitB.model.atomicHierarchy.atoms;
    const { type_symbol: type_symbolA, label_alt_id: label_alt_idA, label_atom_id: label_atom_idA } = unitA.model.atomicHierarchy.atoms;
    const { type_symbol: type_symbolB, label_alt_id: label_alt_idB, label_atom_id: label_atom_idB } = unitB.model.atomicHierarchy.atoms;
    const { label_comp_id: label_comp_idA, auth_seq_id: auth_seq_idA } = unitA.model.atomicHierarchy.residues;
    const { label_comp_id: label_comp_idB, auth_seq_id: auth_seq_idB } = unitB.model.atomicHierarchy.residues;
    const { occupancy: occupancyA } = unitA.model.atomicConformation;
    const { occupancy: occupancyB } = unitB.model.atomicConformation;

    const { lookup3d } = unitB;
    const structConn = unitA.model === unitB.model && unitA.model.sourceData.kind === 'mmCIF' ? StructConn.get(unitA.model) : void 0;
    const indexPairs = unitA.model === unitB.model ? IndexPairBonds.get(unitA.model) : void 0;

    // the lookup queries need to happen in the "unitB space".
    // that means imageA = inverseOperB(operA(aI))
    const imageTransform = Mat4.mul(_imageTransform, unitB.conformation.operator.inverse, unitA.conformation.operator.matrix);
    const isNotIdentity = !Mat4.isIdentity(imageTransform);
    const imageA = Vec3.zero();

    const { center: bCenter, radius: bRadius } = lookup3d.boundary.sphere;
    const testDistanceSq = (bRadius + MAX_RADIUS) * (bRadius + MAX_RADIUS);

    for (let _aI = 0 as StructureElement.UnitIndex; _aI < atomCount; _aI++) {
        const aI = atomsA[_aI];
        Vec3.set(imageA, xA[aI], yA[aI], zA[aI]);
        if (isNotIdentity) Vec3.transformMat4(imageA, imageA, imageTransform);
        if (Vec3.squaredDistance(imageA, bCenter) > testDistanceSq) continue;

        if (!props.forceCompute && indexPairs) {
            for (let i = indexPairs.offset[aI], il = indexPairs.offset[aI + 1]; i < il; ++i) {
                const _bI = SortedArray.indexOf(unitA.elements, indexPairs.b[i]) as StructureElement.UnitIndex;
                if (_bI < 0) continue;
                addBond(_aI, _bI, indexPairs.edgeProps.order[i], BondType.Flag.Covalent, state);
                bondCount++;
            }
            continue // assume `indexPairs` supplies all bonds
        }

        const structConnEntries = props.forceCompute ? void 0 : structConn && structConn.getAtomEntries(aI);
        if (structConnEntries && structConnEntries.length) {
            let added = false;
            for (const se of structConnEntries) {
                for (const p of se.partners) {
                    const _bI = SortedArray.indexOf(unitB.elements, p.atomIndex) as StructureElement.UnitIndex;
                    if (_bI < 0) continue;
                    // check if the bond is within MAX_RADIUS for this pair of units
                    if (getDistance(unitA, aI, unitB, p.atomIndex) > MAX_RADIUS) continue;
                    addBond(_aI, _bI, se.order, se.flags, state);
                    bondCount++;
                    added = true;
                }
            }
            // assume, for an atom, that if any inter unit bond is given
            // all are given and thus we don't need to compute any other
            if (added) continue;
        }

        const { indices, count, squaredDistances } = lookup3d.find(imageA[0], imageA[1], imageA[2], MAX_RADIUS);
        if (count === 0) continue;

        const aeI = getElementIdx(type_symbolA.value(aI));
        const isHa = isHydrogen(aeI);
        const thresholdA = getElementThreshold(aeI);
        const altA = label_alt_idA.value(aI);
        const metalA = MetalsSet.has(aeI);
        const atomIdA = label_atom_idA.value(aI);
        const compIdA = label_comp_idA.value(residueIndexA[aI]);
        const occA = occupancyA.value(aI);

        for (let ni = 0; ni < count; ni++) {
            const _bI = indices[ni] as StructureElement.UnitIndex;
            const bI = atomsB[_bI];

            const altB = label_alt_idB.value(bI);
            if (altA && altB && altA !== altB) continue;

            // Do not include bonds between images of the same residue.
            // TODO: is this condition good enough?
            if (occupancyB.value(bI) < 1 && occA < 1)  {
                if (auth_seq_idA.value(aI) === auth_seq_idB.value(bI)) {
                    continue;
                }
            }

            const beI = getElementIdx(type_symbolB.value(bI)!);
            const isMetal = metalA || MetalsSet.has(beI);

            const isHb = isHydrogen(beI);
            if (isHa && isHb) continue;

            const dist = Math.sqrt(squaredDistances[ni]);
            if (dist === 0) continue;

            if (isHa || isHb) {
                if (dist < props.maxCovalentHydrogenBondingLength) {
                    addBond(
                        _aI, _bI,
                        1, // covalent bonds involving a hydrogen are always of order 1
                        BondType.Flag.Covalent | BondType.Flag.Computed,
                        state
                    );
                    bondCount++;
                }
                continue;
            }

            const thresholdAB = getElementPairThreshold(aeI, beI);
            const pairingThreshold = thresholdAB > 0
                ? thresholdAB
                : beI < 0 ? thresholdA : Math.max(thresholdA, getElementThreshold(beI));

            if (dist <= pairingThreshold) {
                const atomIdB = label_atom_idB.value(bI);
                const compIdB = label_comp_idB.value(residueIndexB[bI]);
                addBond(
                    _aI, _bI,
                    getInterBondOrderFromTable(compIdA, compIdB, atomIdA, atomIdB),
                    (isMetal ? BondType.Flag.MetallicCoordination : BondType.Flag.Covalent) | BondType.Flag.Computed,
                    state
                );
                bondCount++;
            }
        }
    }

    if (bondCount) {
        addMapEntry(map, unitA.id, new InterUnitBonds.UnitPairBonds(unitA, unitB, bondCount, state.bondedA.array, state.mapAB));
        addMapEntry(map, unitB.id, new InterUnitBonds.UnitPairBonds(unitB, unitA, bondCount, state.bondedB.array, state.mapBA));
    }

    return bondCount;
}

export interface InterBondComputationProps extends BondComputationProps {
    validUnitPair: (structure: Structure, unitA: Unit, unitB: Unit) => boolean
}

function findBonds(structure: Structure, props: InterBondComputationProps) {
    const map = new Map<number, InterUnitBonds.UnitPairBonds[]>();

    if (props.noCompute) {
        // TODO add function that only adds bonds defined in structConn and avoids using
        //      structure.lookup and unit.lookup (expensive for large structure and not
        //      needed for archival files or files with an MD topology)
        return new InterUnitBonds(map);
    }

    Structure.eachUnitPair(structure, (unitA: Unit, unitB: Unit) => {
        findPairBonds(unitA as Unit.Atomic, unitB as Unit.Atomic, props, map)
    }, {
        maxRadius: MAX_RADIUS,
        validUnit: (unit: Unit) => Unit.isAtomic(unit),
        validUnitPair: (unitA: Unit, unitB: Unit) => props.validUnitPair(structure, unitA, unitB)
    })

    return new InterUnitBonds(map);
}

function computeInterUnitBonds(structure: Structure, props?: Partial<InterBondComputationProps>): InterUnitBonds {
    return findBonds(structure, {
        ...DefaultBondComputationProps,
        validUnitPair: (props && props.validUnitPair) || Structure.validUnitPair,
    });
}

export { computeInterUnitBonds };
