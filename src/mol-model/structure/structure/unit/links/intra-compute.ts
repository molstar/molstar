/**
 * Copyright (c) 2017-2019 Mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LinkType } from '../../../model/types'
import { IntraUnitLinks } from './data'
import Unit from '../../unit'
import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import { LinkComputationProps, getElementIdx, MetalsSet, getElementThreshold, isHydrogen, getElementPairThreshold, DefaultLinkComputationProps } from './common';
import { SortedArray } from '../../../../../mol-data/int';
import { StructConn, ComponentBond } from '../../../../../mol-model-formats/structure/mmcif/bonds';
import { getIntraBondOrderFromTable } from '../../../../../mol-model/structure/model/properties/atomic/bonds';

function getGraph(atomA: number[], atomB: number[], _order: number[], _flags: number[], atomCount: number): IntraUnitLinks {
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

function _computeBonds(unit: Unit.Atomic, props: LinkComputationProps): IntraUnitLinks {
    const MAX_RADIUS = 4;

    const { x, y, z } = unit.model.atomicConformation;
    const atomCount = unit.elements.length;
    const { elements: atoms, residueIndex } = unit;
    const { type_symbol, label_atom_id, label_alt_id } = unit.model.atomicHierarchy.atoms;
    const { label_comp_id } = unit.model.atomicHierarchy.residues;
    const query3d = unit.lookup3d;

    const structConn = unit.model.sourceData.kind === 'mmCIF' ? StructConn.get(unit.model) : void 0;
    const component = unit.model.sourceData.kind === 'mmCIF' ? ComponentBond.get(unit.model) : void 0;

    const atomA: number[] = [];
    const atomB: number[] = [];
    const flags: number[] = [];
    const order: number[] = [];

    let lastResidue = -1;
    let componentMap: Map<string, Map<string, { flags: number, order: number }>> | undefined = void 0;

    for (let _aI = 0; _aI < atomCount; _aI++) {
        const aI =  atoms[_aI];
        const raI = residueIndex[aI];
        const compId = label_comp_id.value(raI);

        if (!props.forceCompute && raI !== lastResidue) {
            if (!!component && component.entries.has(compId)) {
                componentMap = component.entries.get(compId)!.map;
            } else {
                componentMap = void 0;
            }
        }
        lastResidue = raI;

        const aeI = getElementIdx(type_symbol.value(aI));
        const atomIdA = label_atom_id.value(aI)
        const componentPairs = componentMap ? componentMap.get(atomIdA) : void 0;

        const { indices, count, squaredDistances } = query3d.find(x[aI], y[aI], z[aI], MAX_RADIUS);
        const isHa = isHydrogen(aeI);
        const thresholdA = getElementThreshold(aeI);
        const altA = label_alt_id.value(aI);
        const metalA = MetalsSet.has(aeI);
        const structConnEntries = props.forceCompute ? void 0 : structConn && structConn.getAtomEntries(aI);

        if (structConnEntries) {
            for (const se of structConnEntries) {
                if (se.distance > MAX_RADIUS) continue;

                for (const p of se.partners) {
                    const _bI = SortedArray.indexOf(unit.elements, p.atomIndex);
                    if (_bI < 0) continue;
                    atomA[atomA.length] = _aI;
                    atomB[atomB.length] = _bI;
                    flags[flags.length] = se.flags;
                    order[order.length] = se.order;
                }
            }
        }

        for (let ni = 0; ni < count; ni++) {
            const _bI = indices[ni];
            const bI = atoms[_bI];
            if (bI <= aI) continue;

            const altB = label_alt_id.value(bI);
            if (altA && altB && altA !== altB) continue;

            const beI = getElementIdx(type_symbol.value(bI)!);
            const isMetal = metalA || MetalsSet.has(beI);

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
                        if (flag | LinkType.Flag.Covalent) flag ^= LinkType.Flag.Covalent;
                        flag |= LinkType.Flag.MetallicCoordination;
                    }
                    flags[flags.length] = flag;
                }
                continue;
            }

            const isHb = isHydrogen(beI);
            if (isHa && isHb) continue;

            const dist = Math.sqrt(squaredDistances[ni]);
            if (dist === 0) continue;

            if (isHa || isHb) {
                if (dist < props.maxCovalentHydrogenBondingLength) {
                    atomA[atomA.length] = _aI;
                    atomB[atomB.length] = _bI;
                    order[order.length] = 1; // covalent bonds involving hydrogen are always of order 1
                    flags[flags.length] = LinkType.Flag.Covalent | LinkType.Flag.Computed;
                }
                continue;
            }

            const thresholdAB = getElementPairThreshold(aeI, beI);
            const pairingThreshold = thresholdAB > 0
                ? thresholdAB
                : beI < 0 ? thresholdA : Math.max(thresholdA, getElementThreshold(beI));

            if (dist <= pairingThreshold) {
                atomA[atomA.length] = _aI;
                atomB[atomB.length] = _bI;
                order[order.length] = getIntraBondOrderFromTable(compId, atomIdA, label_atom_id.value(bI));
                flags[flags.length] = (isMetal ? LinkType.Flag.MetallicCoordination : LinkType.Flag.Covalent) | LinkType.Flag.Computed;
            }
        }
    }

    return getGraph(atomA, atomB, order, flags, atomCount);
}

function computeIntraUnitBonds(unit: Unit.Atomic, props?: Partial<LinkComputationProps>) {
    return _computeBonds(unit, { ...DefaultLinkComputationProps, ...props });
}

export { computeIntraUnitBonds }