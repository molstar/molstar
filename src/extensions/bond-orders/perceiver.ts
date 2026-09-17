/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Paul Pillot <paul.pillot@tandemai.com>
 *
 * Placeholder intra-residue bond-order overlay. The real algorithm is ported
 * into `perceiveIntra` later; the signature (structure, unit, bonds, rings, mode) does not change.
 */

import { IntAdjacencyGraph } from '../../mol-math/graph';
import { Structure, StructureElement, Unit } from '../../mol-model/structure';
import { hasIntraBondOrderFromTable } from '../../mol-model/structure/model/properties/atomic/bonds';
import { BondType, WaterNames } from '../../mol-model/structure/model/types';
import { IntraUnitBonds } from '../../mol-model/structure/structure/unit/bonds/data';
import type { UnitRings } from '../../mol-model/structure/structure/unit/rings';

export type BondOrdersMode = 'model' | 'force';

/**
 * Clone `bonds` edge props and overlay placeholder orders on intra-residue covalent edges.
 * Topology (`offset` / `a` / `b`) is unchanged. `structure` and `rings` are unused by the
 * placeholder but are required by the real algorithm.
 */
export function perceiveIntra(structure: Structure, unit: Unit.Atomic, bonds: IntraUnitBonds, rings: UnitRings, mode: BondOrdersMode): IntraUnitBonds {
    void structure;
    void rings;

    const srcOrder = bonds.edgeProps.order;
    const srcFlags = bonds.edgeProps.flags;
    const order = new Int8Array(srcOrder.length);
    const flags = new Uint16Array(srcFlags.length);
    for (let i = 0, il = srcOrder.length; i < il; i++) {
        order[i] = srcOrder[i];
        flags[i] = srcFlags[i];
    }

    const { offset, b } = bonds;
    const { elements, residueIndex } = unit;
    const { type_symbol, label_comp_id } = unit.model.atomicHierarchy.atoms;

    for (let u = 0, ul = elements.length; u < ul; u++) {
        const eU = elements[u];
        const rU = residueIndex[eU];
        const compId = label_comp_id.value(eU);
        if (WaterNames.has(compId) || hasIntraBondOrderFromTable(compId)) continue;

        for (let i = offset[u], il = offset[u + 1]; i < il; i++) {
            const v = b[i];
            if (u >= v) continue;
            if (residueIndex[elements[v]] !== rU) continue;
            if (!BondType.isCovalent(flags[i])) continue;

            if (mode === 'model') {
                if (order[i] !== 1 || !BondType.is(flags[i], BondType.Flag.Computed)) continue;
            }

            if (type_symbol.value(eU) === 'C' && type_symbol.value(elements[v]) === 'C') {
                order[i] = 2;
                const iRev = bonds.getDirectedEdgeIndex(v as StructureElement.UnitIndex, u as StructureElement.UnitIndex);
                if (iRev >= 0) order[iRev] = 2;
            }
        }
    }

    return IntAdjacencyGraph.create(bonds.offset, bonds.a, bonds.b, bonds.edgeCount, {
        ...bonds.edgeProps,
        order,
        flags
    }, bonds.props) satisfies IntraUnitBonds;
}
