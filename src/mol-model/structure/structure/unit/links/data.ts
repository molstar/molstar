/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LinkType } from '../../../model/types'
import { IntAdjacencyGraph } from '../../../../../mol-math/graph';
import Unit from '../../unit';
import StructureElement from '../../element';
import { Link } from '../links';
import { InterUnitGraph } from '../../../../../mol-math/graph/inter-unit-graph';

type IntraUnitLinks = IntAdjacencyGraph<{ readonly order: ArrayLike<number>, readonly flags: ArrayLike<LinkType.Flag> }>

namespace IntraUnitLinks {
    export const Empty: IntraUnitLinks = IntAdjacencyGraph.create([], [], [], 0, { flags: [], order: [] });
}

type InterUnitEdgeProps = { readonly order: number, readonly flag: LinkType.Flag }

class InterUnitBonds extends InterUnitGraph<Unit.Atomic, StructureElement.UnitIndex, InterUnitEdgeProps> {
    /** Get inter-unit bond given a link-location */
    getBondFromLocation(l: Link.Location) {
        return Unit.isAtomic(l.aUnit) && Unit.isAtomic(l.bUnit) ? this.getEdge(l.aIndex, l.aUnit, l.bIndex, l.bUnit) : undefined;
    }

    /** Get inter-unit bond index given a link-location */
    getBondIndexFromLocation(l: Link.Location) {
        return Unit.isAtomic(l.aUnit) && Unit.isAtomic(l.bUnit) ? this.getEdgeIndex(l.aIndex, l.aUnit, l.bIndex, l.bUnit) : -1;
    }
}

namespace InterUnitBonds {
    export class UnitPairBonds extends InterUnitGraph.UnitPairEdges<Unit.Atomic, StructureElement.UnitIndex, InterUnitEdgeProps> {}
    export type BondInfo = InterUnitGraph.EdgeInfo<StructureElement.UnitIndex, InterUnitEdgeProps>
}

export { IntraUnitLinks, InterUnitBonds }