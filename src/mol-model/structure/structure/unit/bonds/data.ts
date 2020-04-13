/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../model/types';
import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import Unit from '../../unit';
import StructureElement from '../../element';
import { Bond } from '../bonds';
import { InterUnitGraph } from '../../../../../mol-math/graph/inter-unit-graph';

type IntraUnitBonds = IntAdjacencyGraph<StructureElement.UnitIndex, { readonly order: ArrayLike<number>, readonly flags: ArrayLike<BondType.Flag> }>

namespace IntraUnitBonds {
    export const Empty: IntraUnitBonds = IntAdjacencyGraph.create([], [], [], 0, { flags: [], order: [] });
}

type InterUnitEdgeProps = { readonly order: number, readonly flag: BondType.Flag }

class InterUnitBonds extends InterUnitGraph<Unit.Atomic, StructureElement.UnitIndex, InterUnitEdgeProps> {
    /** Get inter-unit bond given a bond-location */
    getBondFromLocation(l: Bond.Location) {
        return Unit.isAtomic(l.aUnit) && Unit.isAtomic(l.bUnit) ? this.getEdge(l.aIndex, l.aUnit, l.bIndex, l.bUnit) : undefined;
    }

    /** Get inter-unit bond index given a bond-location */
    getBondIndexFromLocation(l: Bond.Location) {
        return Unit.isAtomic(l.aUnit) && Unit.isAtomic(l.bUnit) ? this.getEdgeIndex(l.aIndex, l.aUnit, l.bIndex, l.bUnit) : -1;
    }
}

namespace InterUnitBonds {
    export class UnitPairBonds extends InterUnitGraph.UnitPairEdges<Unit.Atomic, StructureElement.UnitIndex, InterUnitEdgeProps> {}
    export type BondInfo = InterUnitGraph.EdgeInfo<StructureElement.UnitIndex, InterUnitEdgeProps>
}

export { IntraUnitBonds, InterUnitBonds, InterUnitEdgeProps };