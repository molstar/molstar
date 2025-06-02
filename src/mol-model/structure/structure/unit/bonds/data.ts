/**
 * Copyright (c) 2017-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../model/types';
import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import { Unit } from '../../unit';
import { StructureElement } from '../../element';
import { Bond } from '../bonds';
import { InterUnitGraph } from '../../../../../mol-math/graph/inter-unit-graph';

export type IntraUnitBondProps = {
    /** Can remap even with `dynamicBonds` on, e.g., for water molecules */
    readonly canRemap?: boolean
    /** Can be cached in `ElementSetIntraBondCache` */
    readonly cacheable?: boolean
}

export type IntraUnitBonds = IntAdjacencyGraph<StructureElement.UnitIndex, {
    readonly order: ArrayLike<number>,
    readonly flags: ArrayLike<BondType.Flag>
    readonly key: ArrayLike<number>,
}, IntraUnitBondProps>

export namespace IntraUnitBonds {
    export const Empty: IntraUnitBonds = IntAdjacencyGraph.create([], [], [], 0, { flags: [], order: [], key: [] });
}

export type InterUnitEdgeProps = { readonly order: number, readonly flag: BondType.Flag, readonly key: number }

export class InterUnitBonds extends InterUnitGraph<number, StructureElement.UnitIndex, InterUnitEdgeProps> {
    /** Get inter-unit bond given a bond-location */
    getBondFromLocation(l: Bond.Location) {
        return Unit.isAtomic(l.aUnit) && Unit.isAtomic(l.bUnit) ? this.getEdge(l.aIndex, l.aUnit.id, l.bIndex, l.bUnit.id) : undefined;
    }

    /** Get inter-unit bond index given a bond-location */
    getBondIndexFromLocation(l: Bond.Location) {
        return Unit.isAtomic(l.aUnit) && Unit.isAtomic(l.bUnit) ? this.getEdgeIndex(l.aIndex, l.aUnit.id, l.bIndex, l.bUnit.id) : -1;
    }
}

export namespace InterUnitBonds {
    export class UnitPairBonds extends InterUnitGraph.UnitPairEdges<number, StructureElement.UnitIndex, InterUnitEdgeProps> {}
    export type BondInfo = InterUnitGraph.EdgeInfo<StructureElement.UnitIndex, InterUnitEdgeProps>
}