/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Unit from '../../unit';

const emptyArray: number[] = []

class CrossLinkRestraints {
    readonly count: number
    private readonly pairKeyIndices: Map<string, number[]>

    /** Indices into this.pairs */
    getPairIndices(indexA: number, unitA: Unit, indexB: number, unitB: Unit): ReadonlyArray<number> {
        const key = CrossLinkRestraints.getPairKey(indexA, unitA, indexB, unitB)
        const indices = this.pairKeyIndices.get(key)
        return indices !== undefined ? indices : emptyArray
    }

    getPairs(indexA: number, unitA: Unit, indexB: number, unitB: Unit): CrossLinkRestraints.Pair[] | undefined {
        const indices = this.getPairIndices(indexA, unitA, indexB, unitB)
        return indices.length ? indices.map(idx => this.pairs[idx]) : undefined
    }

    constructor(public pairs: ReadonlyArray<CrossLinkRestraints.Pair>) {
        const pairKeyIndices = new Map<string, number[]>()
        this.pairs.forEach((p, i) => {
            const key = CrossLinkRestraints.getPairKey(p.indexA, p.unitA, p.indexB, p.unitB)
            const indices = pairKeyIndices.get(key)
            if (indices) indices.push(i)
            else pairKeyIndices.set(key, [i])
        })

        this.count = pairs.length
        this.pairKeyIndices = pairKeyIndices
    }
}

namespace CrossLinkRestraints {
    export interface Pair {
        readonly unitA: Unit,
        readonly unitB: Unit,
        readonly indexA: number,
        readonly indexB: number,

        readonly restraintType: 'harmonic' | 'upper bound' | 'lower bound',
        readonly distanceThreshold: number,
        readonly psi: number,
        readonly sigma1: number,
        readonly sigma2: number,
    }

    export function getPairKey(indexA: number, unitA: Unit, indexB: number, unitB: Unit) {
        return `${indexA}|${unitA.id}|${indexB}|${unitB.id}`
    }
}

export { CrossLinkRestraints }