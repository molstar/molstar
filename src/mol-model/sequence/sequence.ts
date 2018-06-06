/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AminoAlphabet, NuclecicAlphabet, getProteinOneLetterCode, getRnaOneLetterCode, getDnaOneLetterCode } from './constants';
import { Column } from 'mol-data/db'

// TODO add mapping support to other sequence spaces, e.g. uniprot
// TODO sequence alignment (take NGL code as starting point)

type Sequence = Sequence.Protein | Sequence.DNA | Sequence.RNA | Sequence.Generic

namespace Sequence {
    export const enum Kind {
        Protein = 'protein',
        RNA = 'RNA',
        DNA = 'DNA',
        Generic = 'generic'
    }

    export interface Base<K extends Kind, Alphabet extends string> {
        readonly kind: K,
        readonly offset: number,
        readonly sequence: ArrayLike<Alphabet>
    }

    export interface Protein extends Base<Kind.Protein, AminoAlphabet> { }
    export interface RNA extends Base<Kind.RNA, NuclecicAlphabet> { }
    export interface DNA extends Base<Kind.DNA, NuclecicAlphabet> { }
    export interface Generic extends Base<Kind.Generic, 'X'> { }

    export function create(kind: Kind, sequence: string, offset: number = 0): Sequence {
        return { kind: kind as any, sequence: sequence as any, offset };
    }

    export function getSequenceString(seq: Sequence) {
        return seq.sequence as string;
    }

    function determineKind(names: Column<string>) {
        for (let i = 0, _i = Math.min(names.rowCount, 10); i < _i; i++) {
            const name = names.value(i) || '';
            if (getProteinOneLetterCode(name) !== 'X') return { kind: Kind.Protein, code: getProteinOneLetterCode };
            if (getRnaOneLetterCode(name) !== 'X') return { kind: Kind.RNA, code: getRnaOneLetterCode };
            if (getDnaOneLetterCode(name) !== 'X') return { kind: Kind.DNA, code: getDnaOneLetterCode };
        }
        return { kind: Kind.Generic, code: (v: string) => 'X' };
    }

    export function ofResidueNames(residueName: Column<string>, seqId: Column<number>): Sequence {
        if (seqId.rowCount === 0) throw new Error('cannot be empty');

        const { kind, code } = determineKind(residueName);
        return new Impl(kind, residueName, seqId, code) as Sequence;
    }

    class Impl implements Base<any, any> {
        private _offset = 0;
        private _seq: string | undefined = void 0;

        get offset() {
            if (this._seq !== void 0) return this._offset;
            this.create();
            return this._offset;
        }

        get sequence(): any {
            if (this._seq !== void 0) return this._seq;
            this.create();
            return this._seq;
        }

        private create() {
            let maxSeqId = 0, minSeqId = Number.MAX_SAFE_INTEGER;
            for (let i = 0, _i = this.seqId.rowCount; i < _i; i++) {
                const id = this.seqId.value(i);
                if (maxSeqId < id) maxSeqId = id;
                if (id < minSeqId) minSeqId = id;
            }

            const count = maxSeqId - minSeqId + 1;
            const sequenceArray = new Array(maxSeqId + 1);
            for (let i = 0; i < count; i++) {
                sequenceArray[i] = 'X';
            }

            for (let i = 0, _i = this.seqId.rowCount; i < _i; i++) {
                sequenceArray[this.seqId.value(i) - minSeqId] = this.code(this.residueName.value(i) || '');
            }

            this._seq = sequenceArray.join('');
            this._offset = minSeqId - 1;
        }

        constructor(public kind: Kind, private residueName: Column<string>, private seqId: Column<number>, private code: (name: string) => string) {

        }
    }
}

export { Sequence }
