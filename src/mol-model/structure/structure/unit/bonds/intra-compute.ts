/**
 * Copyright (c) 2017-2020 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../model/types';
import { IntraUnitBonds } from './data';
import Unit from '../../unit';
import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import { BondComputationProps, getElementIdx, MetalsSet, getElementThreshold, isHydrogen, getElementPairThreshold, DefaultBondComputationProps } from './common';
import { SortedArray } from '../../../../../mol-data/int';
import { getIntraBondOrderFromTable } from '../../../model/properties/atomic/bonds';
import StructureElement from '../../element';
import { IndexPairBonds } from '../../../../../mol-model-formats/structure/property/bonds/index-pair';
import { ComponentBond } from '../../../../../mol-model-formats/structure/property/bonds/comp';
import { StructConn } from '../../../../../mol-model-formats/structure/property/bonds/struct_conn';

function getGraph(atomA: StructureElement.UnitIndex[], atomB: StructureElement.UnitIndex[], _order: number[], _flags: number[], atomCount: number): IntraUnitBonds {
    const builder = new IntAdjacencyGraph.EdgeBuilder(atomCount, atomA, atomB);
    const flags = new Uint16Array(builder.slotCount);
    const order = new Int8Array(builder.slotCount);
    for (let i = 0, _i = builder.edgeCount; i < _i; i++) {
        builder.addNextEdge();
        builder.assignProperty(flags, _flags[i]);
        builder.assignProperty(order, _order[i]);
    }

    return builder.createGraph({ flags, order });
}

const __structConnAdded = new Set<StructureElement.UnitIndex>();

function _computeBonds(unit: Unit.Atomic, props: BondComputationProps): IntraUnitBonds {
    const MAX_RADIUS = 4;

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
    const indexPairs = IndexPairBonds.Provider.get(unit.model);

    const atomA: StructureElement.UnitIndex[] = [];
    const atomB: StructureElement.UnitIndex[] = [];
    const flags: number[] = [];
    const order: number[] = [];

    let lastResidue = -1;
    let componentMap: Map<string, Map<string, { flags: number, order: number }>> | undefined = void 0;

    const structConnAdded = __structConnAdded;

    for (let _aI = 0 as StructureElement.UnitIndex; _aI < atomCount; _aI++) {
        const aI =  atoms[_aI];

        if (!props.forceCompute && indexPairs) {
            const { edgeProps } = indexPairs;
            for (let i = indexPairs.offset[aI], il = indexPairs.offset[aI + 1]; i < il; ++i) {
                const _bI = SortedArray.indexOf(unit.elements, indexPairs.b[i]) as StructureElement.UnitIndex;
                if (_bI < 0) continue;
                if (edgeProps.symmetryA[i] !== edgeProps.symmetryB[i]) continue;
                atomA[atomA.length] = _aI;
                atomB[atomB.length] = _bI;
                order[order.length] = edgeProps.order[i];
                flags[flags.length] = BondType.Flag.Covalent;
            }
            continue; // assume `indexPairs` supplies all bonds
        }

        const structConnEntries = props.forceCompute ? void 0 : structConn && structConn.byAtomIndex.get(aI);
        let hasStructConn = false;
        if (structConnEntries) {
            for (const se of structConnEntries) {
                const { partnerA, partnerB } = se;
                // symmetry must be the same for intra-unit bonds
                if (partnerA.symmetry !== partnerB.symmetry) continue;

                const p = partnerA.atomIndex === aI ? partnerB : partnerA;
                const _bI = SortedArray.indexOf(unit.elements, p.atomIndex) as StructureElement.UnitIndex;
                if (_bI < 0) continue;

                atomA[atomA.length] = _aI;
                atomB[atomB.length] = _bI;
                flags[flags.length] = se.flags;
                order[order.length] = se.order;

                if (!hasStructConn) structConnAdded.clear();
                hasStructConn = true;
                structConnAdded.add(_bI);
            }
        }

        const raI = residueIndex[aI];
        const compId = label_comp_id.value(aI);

        if (!props.forceCompute && raI !== lastResidue) {
            if (!!component && component.entries.has(compId)) {
                const entitySeq = byEntityKey[index.getEntityFromChain(chainIndex[aI])];
                if (entitySeq && entitySeq.sequence.microHet.has(label_seq_id.value(raI))) {
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

        const aeI = getElementIdx(type_symbol.value(aI));
        const atomIdA = label_atom_id.value(aI);
        const componentPairs = componentMap ? componentMap.get(atomIdA) : void 0;

        const { indices, count, squaredDistances } = query3d.find(x[aI], y[aI], z[aI], MAX_RADIUS);
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
                    : (thresholdA + getElementThreshold(beI)) / 2; // not sure if avg or min but max is too big

            if (dist <= pairingThreshold) {
                atomA[atomA.length] = _aI;
                atomB[atomB.length] = _bI;
                order[order.length] = getIntraBondOrderFromTable(compId, atomIdA, label_atom_id.value(bI));
                flags[flags.length] = (isMetal ? BondType.Flag.MetallicCoordination : BondType.Flag.Covalent) | BondType.Flag.Computed;
            }
        }
    }

    return getGraph(atomA, atomB, order, flags, atomCount);
}

function computeIntraUnitBonds(unit: Unit.Atomic, props?: Partial<BondComputationProps>) {
    const p = { ...DefaultBondComputationProps, ...props };
    if (p.noCompute) {
        // TODO add function that only adds bonds defined in structConn of chemCompBond
        //      and avoid using unit.lookup
        return IntraUnitBonds.Empty;
    }
    return _computeBonds(unit, p);
}

export { computeIntraUnitBonds };