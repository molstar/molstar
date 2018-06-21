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
    readonly bondCount: number
    readonly bonds: ReadonlyArray<InterUnitBonds.Bond>
    private readonly bondKeyIndex: Map<string, number>

    getLinkedUnits(unit: Unit): ReadonlyArray<InterUnitBonds.UnitPairBonds> {
        if (!this.map.has(unit.id)) return emptyArray;
        return this.map.get(unit.id)!;
    }

    /** Index into this.bonds */
    getBondIndex(indexA: number, unitA: Unit, indexB: number, unitB: Unit): number {
        const key = InterUnitBonds.getBondKey(indexA, unitA, indexB, unitB)
        const index = this.bondKeyIndex.get(key)
        return index !== undefined ? index : -1
    }

    getBond(indexA: number, unitA: Unit, indexB: number, unitB: Unit): InterUnitBonds.Bond | undefined {
        const index = this.getBondIndex(indexA, unitA, indexB, unitB)
        return index !== -1 ? this.bonds[index] : undefined
    }

    constructor(private map: Map<number, InterUnitBonds.UnitPairBonds[]>) {
        let count = 0
        const bonds: (InterUnitBonds.Bond)[] = []
        const bondKeyIndex = new Map<string, number>()
        this.map.forEach(pairBondsArray => {
            pairBondsArray.forEach(pairBonds => {
                count += pairBonds.bondCount
                pairBonds.linkedElementIndices.forEach(indexA => {
                    pairBonds.getBonds(indexA).forEach(bondInfo => {
                        const { unitA, unitB } = pairBonds
                        const key = InterUnitBonds.getBondKey(indexA, unitA, bondInfo.indexB, unitB)
                        bondKeyIndex.set(key, bonds.length)
                        bonds.push({ ...bondInfo, indexA, unitA, unitB })
                    })
                })
            })
        })
        this.bondCount = count
        this.bonds = bonds
        this.bondKeyIndex = bondKeyIndex
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

    export interface Bond {
        readonly unitA: Unit.Atomic,
        readonly unitB: Unit.Atomic,
        readonly indexA: number,
        readonly indexB: number,
        readonly order: number,
        readonly flag: LinkType.Flag
    }

    export function getBondKey(indexA: number, unitA: Unit, indexB: number, unitB: Unit) {
        return `${indexA}|${unitA.id}|${indexB}|${unitB.id}`
    }
}

const emptyArray: any[] = [];

export { IntraUnitLinks, InterUnitBonds }