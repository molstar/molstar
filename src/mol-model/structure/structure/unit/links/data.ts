/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LinkType } from '../../../model/types'
import { IntAdjacencyGraph } from 'mol-math/graph';
import Unit from '../../unit';

type IntraUnitLinks = IntAdjacencyGraph<{ readonly order: ArrayLike<number>, readonly flags: ArrayLike<LinkType.Flag> }>

namespace IntraUnitLinks {
    export const Empty: IntraUnitLinks = IntAdjacencyGraph.create([], [], [], 0, { flags: [], order: [] });
}

class InterUnitBonds {
    getLinkedUnits(unit: Unit): ReadonlyArray<InterUnitBonds.UnitPairBonds> {
        if (!this.map.has(unit.id)) return emptyArray;
        return this.map.get(unit.id)!;
    }

    constructor(private map: Map<number, InterUnitBonds.UnitPairBonds[]>) {
    }
}

namespace InterUnitBonds {
    export class UnitPairBonds {
        hasBonds(indexA: number) {
            return this.linkMap.has(indexA);
        }

        getBonds(indexA: number): ReadonlyArray<InterUnitBonds.BondInfo> {
            if (!this.linkMap.has(indexA)) return emptyArray;
            return this.linkMap.get(indexA)!;
        }

        get areUnitsOrdered() {
            return this.unitA.id < this.unitB.id;
        }

        constructor(public unitA: Unit.Atomic, public unitB: Unit.Atomic,
            public bondCount: number, public linkedElementIndices: ReadonlyArray<number>,
            private linkMap: Map<number, BondInfo[]>) {
        }
    }

    export interface BondInfo {
        /** indexInto */
        readonly indexB: number,
        readonly order: number,
        readonly flag: LinkType.Flag
    }
}

const emptyArray: any[] = [];

export { IntraUnitLinks, InterUnitBonds }