/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import Unit from '../../unit';
import { StructureElement } from '../../../structure';

const emptyArray: number[] = []

interface PairRestraint {
    readonly unitA: Unit,
    readonly unitB: Unit,
    readonly indexA: StructureElement.UnitIndex,
    readonly indexB: StructureElement.UnitIndex,
}

function getPairKey(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit) {
    return `${indexA}|${unitA.id}|${indexB}|${unitB.id}`
}

export class PairRestraints<T extends PairRestraint> {
    readonly count: number
    private readonly pairKeyIndices: Map<string, number[]>

    /** Indices into this.pairs */
    getPairIndices(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit): ReadonlyArray<number> {
        const key = getPairKey(indexA, unitA, indexB, unitB)
        const indices = this.pairKeyIndices.get(key)
        return indices !== undefined ? indices : emptyArray
    }

    getPairs(indexA: StructureElement.UnitIndex, unitA: Unit, indexB: StructureElement.UnitIndex, unitB: Unit): T[] | undefined {
        const indices = this.getPairIndices(indexA, unitA, indexB, unitB)
        return indices.length ? indices.map(idx => this.pairs[idx]) : undefined
    }

    constructor(public pairs: ReadonlyArray<T>) {
        const pairKeyIndices = new Map<string, number[]>()
        this.pairs.forEach((p, i) => {
            const key = getPairKey(p.indexA, p.unitA, p.indexB, p.unitB)
            const indices = pairKeyIndices.get(key)
            if (indices) indices.push(i)
            else pairKeyIndices.set(key, [i])
        })

        this.count = pairs.length
        this.pairKeyIndices = pairKeyIndices
    }
}

export interface CrossLinkRestraint extends PairRestraint {
    readonly restraintType: 'harmonic' | 'upper bound' | 'lower bound'
    readonly distanceThreshold: number
    readonly psi: number
    readonly sigma1: number
    readonly sigma2: number
}

export interface PredictedContactRestraint extends PairRestraint {
    readonly distance_lower_limit: number
    readonly distance_upper_limit: number
    readonly probability: number
    readonly restraint_type: 'lower bound' | 'upper bound' | 'lower and upper bound'
    readonly model_granularity: 'by-residue' | 'by-feature' | 'by-atom'
}

export interface DistanceRestraint extends PairRestraint {
    readonly upper_limit: number
    readonly upper_limit_esd: number
    readonly lower_limit: number
    readonly lower_limit_esd: number
    readonly probability: number
    readonly restraint_type: 'lower bound' | 'upper bound' | 'lower and upper bound'
    readonly granularity: 'by-residue' | 'by-atom'
}