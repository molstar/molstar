/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../../model/types'

interface GroupBonds {
    /**
     * Where bonds for atom A start and end.
     * Start offset at idx, end at idx + 1
     */
    offset: ArrayLike<number>,
    neighbor: ArrayLike<number>,

    order: ArrayLike<number>,
    flags: ArrayLike<BondType.Flag>,

    count: number
}

namespace GroupBonds {
    export function createEmpty(): GroupBonds {
        return { offset: [], neighbor: [], order: [], flags: [], count: 0 }
    }
    export function isCovalent(flags: number) {
        return (flags & BondType.Flag.Covalent) !== 0;
    }
    export interface StructConnEntry {
        flags: BondType.Flag,
        order: number,
        distance: number,
        partners: { residueIndex: number, atomIndex: number, symmetry: string }[]
    }
    export interface StructConn {
        getResidueEntries(residueAIndex: number, residueBIndex: number): ReadonlyArray<StructConnEntry>
        getAtomEntries(atomIndex: number): ReadonlyArray<StructConnEntry>
    }
    export interface ComponentBondInfoEntry {
        map: Map<string, Map<string, { order: number, flags: number }>>
    }
    export interface ComponentBondInfo {
        entries: Map<string, ComponentBondInfoEntry>
    }
}

export { GroupBonds }