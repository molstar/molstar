/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BondType } from '../../../model/types'
import { IntGraph } from 'mol-math/graph';

type IntraUnitBonds = IntGraph<{ readonly order: ArrayLike<number>, readonly flags: ArrayLike<BondType.Flag> }>

namespace IntraUnitBonds {
    export function createEmpty(): IntraUnitBonds {
        return IntGraph.create([], [], [], 0, { flags: [], order: [] });
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

export { IntraUnitBonds }