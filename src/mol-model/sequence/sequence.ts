/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AminoAlphabet, NuclecicAlphabet, getProteinOneLetterCode, getRnaOneLetterCode, getDnaOneLetterCode } from './constants';
import { Column } from '../../mol-data/db'

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
        readonly labels: ArrayLike<string>
        /** maps seqId to list of compIds */
        readonly microHet: ReadonlyMap<number, string[]>
    }

    export interface Protein extends Base<Kind.Protein, AminoAlphabet> { }
    export interface RNA extends Base<Kind.RNA, NuclecicAlphabet> { }
    export interface DNA extends Base<Kind.DNA, NuclecicAlphabet> { }
    export interface Generic extends Base<Kind.Generic, 'X' | '-'> { }

    export function create<K extends Kind, Alphabet extends string>(kind: K, sequence: Alphabet[], labels: string[], microHet: Map<number, string[]>, offset: number = 0): Base<K, Alphabet> {
        return { kind: kind, sequence: sequence, labels, microHet, offset };
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

    function modCode(code: (name: string) => string, map: ReadonlyMap<string, string>): (name: string) => string {
        return n => {
            const ret = code(n);
            if (ret !== 'X' || !map.has(n)) return ret;
            return code(map.get(n)!);
        }
    }

    export function ofResidueNames(residueName: Column<string>, seqId: Column<number>, modifiedMap?: ReadonlyMap<string, string>): Sequence {
        if (seqId.rowCount === 0) throw new Error('cannot be empty');

        const { kind, code } = determineKind(residueName);

        if (!modifiedMap || modifiedMap.size === 0) {
            return new Impl(kind, residueName, seqId, code) as Sequence;
        }
        return new Impl(kind, residueName, seqId, modCode(code, modifiedMap)) as Sequence;
    }

    class Impl<K extends Kind, Alphabet extends string> implements Base<K, Alphabet> {
        private _offset = 0;
        private _seq: ArrayLike<Alphabet> | undefined = void 0;
        private _labels: ArrayLike<string> | undefined = void 0;
        private _microHet: ReadonlyMap<number, string[]> | undefined = void 0;

        get offset() {
            if (this._seq !== void 0) return this._offset;
            this.create();
            return this._offset;
        }

        get sequence(): ArrayLike<Alphabet> {
            if (this._seq !== void 0) return this._seq;
            this.create();
            return this._seq!;
        }

        get labels(): ArrayLike<string> {
            if (this._labels !== void 0) return this._labels;
            this.create();
            return this._labels!;
        }

        get microHet(): ReadonlyMap<number, string[]> {
            if (this._microHet !== void 0) return this._microHet;
            this.create();
            return this._microHet!;
        }

        private create() {
            let maxSeqId = 0, minSeqId = Number.MAX_SAFE_INTEGER;
            for (let i = 0, _i = this.seqId.rowCount; i < _i; i++) {
                const id = this.seqId.value(i);
                if (maxSeqId < id) maxSeqId = id;
                if (id < minSeqId) minSeqId = id;
            }

            const count = maxSeqId - minSeqId + 1;
            const sequenceArray = new Array<string>(maxSeqId + 1);
            const labels = new Array<string[]>(maxSeqId + 1);
            for (let i = 0; i < count; i++) {
                sequenceArray[i] = '-';
                labels[i] = [];
            }

            const compIds = new Array<string[]>(maxSeqId + 1);
            for (let i = minSeqId; i <= maxSeqId; ++i) {
                compIds[i] = [];
            }

            for (let i = 0, _i = this.seqId.rowCount; i < _i; i++) {
                const seqId = this.seqId.value(i)
                const idx = seqId - minSeqId;
                const name = this.residueName.value(i);
                const code = this.code(name);
                // in case of MICROHETEROGENEITY `sequenceArray[idx]` may already be set
                if (!sequenceArray[idx] || sequenceArray[idx] === '-') {
                    sequenceArray[idx] = code;
                }
                labels[idx].push(code === 'X' ? name : code);
                compIds[seqId].push(name);
            }

            const microHet = new Map()
            for (let i = minSeqId; i <= maxSeqId; ++i) {
                if (compIds[i].length > 1) microHet.set(i, compIds[i])
            }

            this._seq = sequenceArray.join('') as unknown as ArrayLike<Alphabet>;
            this._labels = labels.map(l => l.length > 1 ? `(${l.join('|')})` : l.join(''));
            this._microHet = microHet
            this._offset = minSeqId - 1;
        }

        constructor(public kind: K, private residueName: Column<string>, private seqId: Column<number>, private code: (name: string) => string) {

        }
    }
}

export { Sequence }
