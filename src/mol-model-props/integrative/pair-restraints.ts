/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureElement, Unit } from '../../mol-model/structure';

const emptyArray: number[] = [];

export interface PairRestraint {
    readonly unitA: Unit,
    readonly unitB: Unit,
    readonly indexA: StructureElement.UnitIndex,
    readonly indexB: StructureElement.UnitIndex,
}

function getPairKey(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit) {
    return `${indexA}|${unitA.id}|${indexB}|${unitB.id}`;
}

export class PairRestraints<T extends PairRestraint> {
    readonly count: number
    private readonly pairKeyIndices: Map<string, number[]>

    /** Indices into this.pairs */
    getPairIndices(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit): ReadonlyArray<number> {
        const key = getPairKey(indexA, unitA, indexB, unitB);
        return this.pairKeyIndices.get(key) || emptyArray;
    }

    getPairs(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit): T[] {
        const indices = this.getPairIndices(indexA, unitA, indexB, unitB);
        return indices.map(idx => this.pairs[idx]);
    }

    constructor(public pairs: ReadonlyArray<T>) {
        const pairKeyIndices = new Map<string, number[]>();
        this.pairs.forEach((p, i) => {
            const key = getPairKey(p.indexA, p.unitA, p.indexB, p.unitB);
            const indices = pairKeyIndices.get(key);
            if (indices) indices.push(i);
            else pairKeyIndices.set(key, [i]);
        });

        this.count = pairs.length;
        this.pairKeyIndices = pairKeyIndices;
    }
}