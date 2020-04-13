/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { BitFlags } from '../../../../mol-util';
import { SecondaryStructureType } from '../../../../mol-model/structure/model/types';
import { IntAdjacencyGraph } from '../../../../mol-math/graph';
import { Unit } from '../../../../mol-model/structure/structure';
import { ProteinInfo } from './protein-info';

export interface DSSPContext {
    params: {
        oldDefinition: boolean
        oldOrdering: boolean
    },
    getResidueFlag: (f: DSSPType) => SecondaryStructureType,
    getFlagName: (f: DSSPType) => string,

    unit: Unit.Atomic
    proteinInfo: ProteinInfo
    /** flags for each residue */
    flags: Uint32Array
    hbonds: DsspHbonds,

    torsionAngles: { phi: Float32Array, psi: Float32Array },
    ladders: Ladder[],
    bridges: Bridge[]
}

export { DSSPType };
type DSSPType = BitFlags<DSSPType.Flag>
namespace DSSPType {
    export const is: (t: DSSPType, f: Flag) => boolean = BitFlags.has;
    export const create: (f: Flag) => DSSPType = BitFlags.create;
    export const enum Flag {
        _ = 0x0,
        H = 0x1,
        B = 0x2,
        E = 0x4,
        G = 0x8,
        I = 0x10,
        S = 0x20,
        T = 0x40,
        T3 = 0x80,
        T4 = 0x100,
        T5 = 0x200,
        T3S = 0x400, // marks 3-turn start
        T4S = 0x800,
        T5S = 0x1000
    }
}

export type DsspHbonds = IntAdjacencyGraph<number, { readonly energies: ArrayLike<number> }>

export interface Ladder {
    previousLadder: number,
    nextLadder: number,
    firstStart: number,
    secondStart: number,
    secondEnd: number,
    firstEnd: number,
    type: BridgeType
}

export const enum BridgeType {
    PARALLEL = 0x0,
    ANTI_PARALLEL = 0x1
}

export class Bridge {
    partner1: number;
    partner2: number;
    type: BridgeType;

    constructor(p1: number, p2: number, type: BridgeType) {
        this.partner1 = Math.min(p1, p2);
        this.partner2 = Math.max(p1, p2);
        this.type = type;
    }
}