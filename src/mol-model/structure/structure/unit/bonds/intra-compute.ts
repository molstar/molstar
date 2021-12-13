/**
 * Copyright (c) 2017-2021 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../model/types';
import { IntraUnitBonds } from './data';
import { Unit } from '../../unit';
import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import { BondComputationProps, getElementIdx, MetalsSet, getElementThreshold, isHydrogen, getElementPairThreshold, DefaultBondComputationProps } from './common';
import { SortedArray } from '../../../../../mol-data/int';
import { getIntraBondOrderFromTable } from '../../../model/properties/atomic/bonds';
import { StructureElement } from '../../element';
import { IndexPairBonds } from '../../../../../mol-model-formats/structure/property/bonds/index-pair';
import { ComponentBond } from '../../../../../mol-model-formats/structure/property/bonds/chem_comp';
import { StructConn } from '../../../../../mol-model-formats/structure/property/bonds/struct_conn';
import { Vec3 } from '../../../../../mol-math/linear-algebra';
import { ElementIndex } from '../../../model/indexing';
import { equalEps } from '../../../../../mol-math/linear-algebra/3d/common';
import { Model } from '../../../model/model';

function getGraph(atomA: StructureElement.UnitIndex[], atomB: StructureElement.UnitIndex[], _order: number[], _flags: number[], atomCount: number, canRemap: boolean): IntraUnitBonds {
    const builder = new IntAdjacencyGraph.EdgeBuilder(atomCount, atomA, atomB);
    const flags = new Uint16Array(builder.slotCount);
    const order = new Int8Array(builder.slotCount);
    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
        builder.addNextEdge();
        builder.assignProperty(flags, _flags[i]);
        builder.assignProperty(order, _order[i]);
    }

    return builder.createGraph({ flags, order }, { canRemap });
}

const tmpDistVecA = Vec3();
const tmpDistVecB = Vec3();
function getDistance(unit: Unit.Atomic, indexA: ElementIndex, indexB: ElementIndex) {
    unit.conformation.position(indexA, tmpDistVecA);
    unit.conformation.position(indexB, tmpDistVecB);
    return Vec3.distance(tmpDistVecA, tmpDistVecB);
}

const __structConnAdded = new Set<StructureElement.UnitIndex>();

function findIndexPairBonds(unit: Unit.Atomic) {
    const indexPairs = IndexPairBonds.Provider.get(unit.model)!;
    const { elements: atoms } = unit;
    const { type_symbol } = unit.model.atomicHierarchy.atoms;
    const atomCount = unit.elements.length;
    const { maxDistance } = indexPairs;
    const { offset, b, edgeProps: { order, distance, flag } } = indexPairs.bonds;

    const { atomSourceIndex: sourceIndex } = unit.model.atomicHierarchy;
    const { invertedIndex } = Model.getInvertedAtomSourceIndex(unit.model);

    const atomA: StructureElement.UnitIndex[] = [];
    const atomB: StructureElement.UnitIndex[] = [];
    const flags: number[] = [];
    const orders: number[] = [];

    for (let _aI = 0 as StructureElement.UnitIndex; _aI < atomCount; _aI++) {
        const aI = atoms[_aI];
        const isHa = type_symbol.value(aI) === 'H';

        const srcA = sourceIndex.value(aI);

        for (let i = offset[srcA], il = offset[srcA + 1]; i < il; ++i) {
            const bI = invertedIndex[b[i]];
            if (aI >= bI) continue;

            const _bI = SortedArray.indexOf(unit.elements, bI) as StructureElement.UnitIndex;
            if (_bI < 0) continue;
            if (isHa && type_symbol.value(bI) === 'H') continue;

            const d = distance[i];
            const dist = getDistance(unit, aI, bI);
            if ((d !== -1 && equalEps(dist, d, 0.5)) || dist < maxDistance) {
                atomA[atomA.length] = _aI;
                atomB[atomB.length] = _bI;
                orders[orders.length] = order[i];
                flags[flags.length] = flag[i];
            }
        }
    }

    return getGraph(atomA, atomB, orders, flags, atomCount, false);
}

function findBonds(unit: Unit.Atomic, props: BondComputationProps): IntraUnitBonds {
    const { maxRadius } = props;

    const { x, y, z } = unit.model.atomicConformation;
    const atomCount = unit.elements.length;
    const { elements: atoms, residueIndex, chainIndex } = unit;
    const { type_symbol, label_atom_id, label_alt_id, label_comp_id } = unit.model.atomicHierarchy.atoms;
    const { label_seq_id } = unit.model.atomicHierarchy.residues;
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

    let lastResidue = -1;
    let componentMap: Map<string, Map<string, { flags: number, order: number }>> | undefined = void 0;

    let isWatery = true, isDictionaryBased = true, isSequenced = true;

    const structConnAdded = __structConnAdded;

    for (let _aI = 0 as StructureElement.UnitIndex; _aI < atomCount; _aI++) {
        const aI = atoms[_aI];

        const elemA = type_symbol.value(aI);
        if (isWatery && (elemA !== 'H' || elemA !== 'O')) isWatery = false;

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
                    componentMap = void 0;
                } else {
                    componentMap = component.entries.get(compId)!.map;
                }
            } else {
                componentMap = void 0;
            }
        }
        lastResidue = raI;

        const aeI = getElementIdx(elemA);
        const atomIdA = label_atom_id.value(aI);
        const componentPairs = componentMap ? componentMap.get(atomIdA) : void 0;

        const { indices, count, squaredDistances } = query3d.find(x[aI], y[aI], z[aI], maxRadius);
        const isHa = isHydrogen(aeI);
        const thresholdA = getElementThreshold(aeI);
        const altA = label_alt_id.value(aI);
        const metalA = MetalsSet.has(aeI);

        for (let ni = 0; ni < count; ni++) {
            const _bI = indices[ni];
            if (hasStructConn && structConnAdded.has(_bI)) continue;

            const bI = atoms[_bI];
            if (bI <= aI) continue;

            const altB = label_alt_id.value(bI);
            if (altA && altB && altA !== altB) continue;

            const beI = getElementIdx(type_symbol.value(bI)!);

            const isHb = isHydrogen(beI);
            if (isHa && isHb) continue;

            const isMetal = (metalA || MetalsSet.has(beI)) && !(isHa || isHb);

            const rbI = residueIndex[bI];
            // handle "component dictionary" bonds.
            if (raI === rbI && componentPairs) {
                const e = componentPairs.get(label_atom_id.value(bI)!);
                if (e) {
                    atomA[atomA.length] = _aI;
                    atomB[atomB.length] = _bI;
                    order[order.length] = e.order;
                    let flag = e.flags;
                    if (isMetal) {
                        if (flag | BondType.Flag.Covalent) flag ^= BondType.Flag.Covalent;
                        flag |= BondType.Flag.MetallicCoordination;
                    }
                    flags[flags.length] = flag;
                }
                continue;
            }

            const dist = Math.sqrt(squaredDistances[ni]);
            if (dist === 0) continue;

            const thresholdAB = getElementPairThreshold(aeI, beI);
            const pairingThreshold = thresholdAB > 0
                ? thresholdAB
                : beI < 0
                    ? thresholdA
                    : (thresholdA + getElementThreshold(beI)) / 1.95; // not sure if avg or min but max is too big

            if (dist <= pairingThreshold) {
                atomA[atomA.length] = _aI;
                atomB[atomB.length] = _bI;
                order[order.length] = getIntraBondOrderFromTable(compId, atomIdA, label_atom_id.value(bI));
                flags[flags.length] = (isMetal ? BondType.Flag.MetallicCoordination : BondType.Flag.Covalent) | BondType.Flag.Computed;

                const seqIdB = label_seq_id.value(rbI);

                if (seqIdA === seqIdB) isDictionaryBased = false;
                if (Math.abs(seqIdA - seqIdB) > 1) isSequenced = false;
            }
        }
    }

    const canRemap = isWatery || (isDictionaryBased && isSequenced);
    return getGraph(atomA, atomB, order, flags, atomCount, canRemap);
}

function computeIntraUnitBonds(unit: Unit.Atomic, props?: Partial<BondComputationProps>) {
    const p = { ...DefaultBondComputationProps, ...props };
    if (p.noCompute || Model.isCoarseGrained(unit.model)) {
        // TODO add function that only adds bonds defined in structConn of chemCompBond
        //      and avoid using unit.lookup
        return IntraUnitBonds.Empty;
    }

    if (!p.forceCompute && IndexPairBonds.Provider.get(unit.model)) {
        return findIndexPairBonds(unit);
    } else {
        return findBonds(unit, p);
    }
}

export { computeIntraUnitBonds };