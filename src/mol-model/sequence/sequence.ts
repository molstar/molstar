/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// import { Column } from 'mol-data/db'
// TODO

interface Sequence {

}

namespace Sequence {
    export const enum Kind {
        Protein = 'protein',
        RNA = 'RNA',
        DNA = 'DNA',
        Generic = 'generic'
    }

    export type ProteinAlphabet = 'X'
    export type RnaAlphabet = 'X'
    export type DnaAlphabet = 'C'
    export type GenericAlphabet = 'X'

    export interface Base<K extends Kind, Alphabet extends string> {
        readonly kind: K,
        readonly sequence: ReadonlyArray<Alphabet>
    }
}

export { Sequence }
