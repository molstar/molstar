/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { LinkType } from '../../../model/types'
import { IntAdjacencyGraph } from 'mol-math/graph';
import Unit from '../../unit';
import StructureElement from '../../element';
import { Link } from '../links';

type IntraUnitLinks = IntAdjacencyGraph<{ readonly order: ArrayLike<number>, readonly flags: ArrayLike<LinkType.Flag> }>

namespace IntraUnitLinks {
    export const Empty: IntraUnitLinks = IntAdjacencyGraph.create([], [], [], 0, { flags: [], order: [] });
}

class InterUnitBonds {
    /** Number of inter-unit bonds */
    readonly bondCount: number
    /** Array of inter-unit bonds */
    readonly bonds: ReadonlyArray<InterUnitBonds.Bond>
    private readonly bondKeyIndex: Map<string, number>
    private readonly elementKeyIndex: Map<string, number[]>

    /** Get an array of unit-pair-bonds that are linked to the given unit */
    getLinkedUnits(unit: Unit): ReadonlyArray<InterUnitBonds.UnitPairBonds> {
        if (!this.map.has(unit.id)) return emptyArray;
        return this.map.get(unit.id)!;
    }

    /** Index into this.bonds */
    getBondIndex(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit): number {
        const bondKey = InterUnitBonds.getBondKey(indexA, unitA, indexB, unitB)
        const index = this.bondKeyIndex.get(bondKey)
        return index !== undefined ? index : -1
    }

    /** Get inter-unit bond given a pair of indices and units */
    getBond(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit): InterUnitBonds.Bond | undefined {
        const index = this.getBondIndex(indexA, unitA, indexB, unitB)
        return index !== -1 ? this.bonds[index] : undefined
    }

    /** Get inter-unit bond given a link-location */
    getBondFromLocation(l: Link.Location) {
        return this.getBond(l.aIndex, l.aUnit, l.bIndex, l.bUnit);
    }

    /** Indices into this.bonds */
    getBondIndices(index: StructureElement.UnitIndex, unit: Unit): ReadonlyArray<number> {
        const elementKey = InterUnitBonds.getElementKey(index, unit)
        const indices = this.elementKeyIndex.get(elementKey)
        return indices !== undefined ? indices : []
    }

    constructor(private map: Map<number, InterUnitBonds.UnitPairBonds[]>) {
        let count = 0
        const bonds: (InterUnitBonds.Bond)[] = []
        const bondKeyIndex = new Map<string, number>()
        const elementKeyIndex = new Map<string, number[]>()

        this.map.forEach(pairBondsArray => {
            pairBondsArray.forEach(pairBonds => {
                count += pairBonds.bondCount
                pairBonds.linkedElementIndices.forEach(indexA => {
                    pairBonds.getBonds(indexA).forEach(bondInfo => {
                        const { unitA, unitB } = pairBonds
                        
                        const bondKey = InterUnitBonds.getBondKey(indexA, unitA, bondInfo.indexB, unitB)
                        bondKeyIndex.set(bondKey, bonds.length)
                        
                        const elementKey = InterUnitBonds.getElementKey(indexA, unitA)
                        const e = elementKeyIndex.get(elementKey)
                        if (e === undefined) elementKeyIndex.set(elementKey, [bonds.length])
                        else e.push(bonds.length)

                        bonds.push({ ...bondInfo, indexA, unitA, unitB })
                    })
                })
            })
        })

        this.bondCount = count
        this.bonds = bonds
        this.bondKeyIndex = bondKeyIndex
        this.elementKeyIndex = elementKeyIndex
    }
}

namespace InterUnitBonds {
    export class UnitPairBonds {
        hasBonds(indexA: StructureElement.UnitIndex) {
            return this.linkMap.has(indexA);
        }

        getBonds(indexA: StructureElement.UnitIndex): ReadonlyArray<InterUnitBonds.BondInfo> {
            if (!this.linkMap.has(indexA)) return emptyArray;
            return this.linkMap.get(indexA)!;
        }

        get areUnitsOrdered() {
            return this.unitA.id < this.unitB.id;
        }

        constructor(public unitA: Unit.Atomic, public unitB: Unit.Atomic,
            public bondCount: number, public linkedElementIndices: ReadonlyArray<StructureElement.UnitIndex>,
            private linkMap: Map<number, BondInfo[]>) {
        }
    }

    export interface BondInfo {
        /** indexInto */
        readonly indexB: StructureElement.UnitIndex,
        readonly order: number,
        readonly flag: LinkType.Flag
    }

    export interface Bond {
        readonly unitA: Unit.Atomic,
        readonly unitB: Unit.Atomic,
        readonly indexA: StructureElement.UnitIndex,
        readonly indexB: StructureElement.UnitIndex,
        readonly order: number,
        readonly flag: LinkType.Flag
    }

    export function getBondKey(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit) {
        return `${indexA}|${unitA.id}|${indexB}|${unitB.id}`
    }

    export function getElementKey(index: StructureElement.UnitIndex, unit: Unit) {
        return `${index}|${unit.id}`
    }
}

const emptyArray: any[] = [];

export { IntraUnitLinks, InterUnitBonds }