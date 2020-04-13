/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AminoAlphabet, NuclecicAlphabet, getProteinOneLetterCode, getRnaOneLetterCode, getDnaOneLetterCode } from './constants';
import { Column } from '../../mol-data/db';

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
        readonly length: number,
        readonly offset: number,

        readonly code: Column<Alphabet>
        readonly label: Column<string>

        readonly seqId: Column<number>
        readonly compId: Column<string>

        /** maps seqId to list of compIds */
        readonly microHet: ReadonlyMap<number, string[]>
    }

    export interface Protein extends Base<Kind.Protein, AminoAlphabet> { }
    export interface RNA extends Base<Kind.RNA, NuclecicAlphabet> { }
    export interface DNA extends Base<Kind.DNA, NuclecicAlphabet> { }
    export interface Generic extends Base<Kind.Generic, 'X' | '-'> { }

    export function create<K extends Kind, Alphabet extends string>(kind: K, code: Column<Alphabet>, label: Column<string>, seqId: Column<number>, compId: Column<string>, microHet: Map<number, string[]>, offset: number = 0): Base<K, Alphabet> {
        const length = code.rowCount;
        return { kind, code, label, seqId, compId, microHet, offset, length };
    }

    export function getSequenceString(seq: Sequence) {
        const array = seq.code.toArray();
        return (array instanceof Array ? array : Array.from(array)).join('');
    }

    function determineKind(names: Column<string>) {
        for (let i = 0, _i = Math.min(names.rowCount, 10); i < _i; i++) {
            const name = names.value(i) || '';
            if (getProteinOneLetterCode(name) !== 'X') return Kind.Protein;
            if (getRnaOneLetterCode(name) !== 'X') return Kind.RNA;
            if (getDnaOneLetterCode(name) !== 'X') return Kind.DNA;
        }
        return Kind.Generic;
    }

    function codeProvider(kind: Kind, map?: ReadonlyMap<string, string>) {
        let code: (name: string) => string;
        switch (kind) {
            case Kind.Protein: code = getProteinOneLetterCode; break;
            case Kind.DNA: code = getDnaOneLetterCode; break;
            case Kind.RNA: code = getRnaOneLetterCode; break;
            case Kind.Generic: code = () => 'X'; break;
            default: throw new Error(`unknown kind '${kind}'`);
        }
        if (map && map.size > 0) {
            return (name: string) => {
                const ret = code(name);
                if (ret !== 'X' || !map.has(name)) return ret;
                return code(map.get(name)!);
            };
        }
        return code;
    }

    export function ofResidueNames(compId: Column<string>, seqId: Column<number>): Sequence {
        if (seqId.rowCount === 0) throw new Error('cannot be empty');

        const kind = determineKind(compId);
        return new ResidueNamesImpl(kind, compId, seqId) as Sequence;
    }

    class ResidueNamesImpl<K extends Kind, Alphabet extends string> implements Base<K, Alphabet> {
        private _offset = 0;
        private _length = 0;
        private _microHet: ReadonlyMap<number, string[]> | undefined = void 0;
        private _code: Column<Alphabet> | undefined = undefined
        private _label: Column<string> | undefined = undefined

        private codeFromName: (name: string) => string

        get code(): Column<Alphabet> {
            if (this._code !== void 0) return this._code;
            this.create();
            return this._code!;
        }

        get label(): Column<string> {
            if (this._label !== void 0) return this._label;
            this.create();
            return this._label!;
        }

        get offset() {
            if (this._code !== void 0) return this._offset;
            this.create();
            return this._offset;
        }

        get length() {
            if (this._code !== void 0) return this._length;
            this.create();
            return this._length;
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
                const seqId = this.seqId.value(i);
                const idx = seqId - minSeqId;
                const name = this.compId.value(i);
                const code = this.codeFromName(name);
                // in case of MICROHETEROGENEITY `sequenceArray[idx]` may already be set
                if (!sequenceArray[idx] || sequenceArray[idx] === '-') {
                    sequenceArray[idx] = code;
                }
                labels[idx].push(code === 'X' ? name : code);
                compIds[seqId].push(name);
            }

            const microHet = new Map();
            for (let i = minSeqId; i <= maxSeqId; ++i) {
                if (compIds[i].length > 1) microHet.set(i, compIds[i]);
            }

            this._code = Column.ofStringArray(sequenceArray) as Column<Alphabet>;
            this._label = Column.ofLambda({
                value: i => {
                    const l = labels[i];
                    return l.length > 1 ? `(${l.join('|')})` : l.join('');
                },
                rowCount: labels.length,
                schema: Column.Schema.str
            });
            this._microHet = microHet;
            this._offset = minSeqId - 1;
            this._length = count;
        }

        constructor(public kind: K, public compId: Column<string>, public seqId: Column<number>) {

            this.codeFromName = codeProvider(kind);
        }
    }

    export function ofSequenceRanges(seqIdBegin: Column<number>, seqIdEnd: Column<number>): Sequence {
        const kind = Kind.Generic;

        return new SequenceRangesImpl(kind, seqIdBegin, seqIdEnd) as Sequence;
    }

    class SequenceRangesImpl<K extends Kind, Alphabet extends string> implements Base<K, Alphabet> {
        public offset: number
        public length: number
        public code: Column<Alphabet>
        public label: Column<string>
        public seqId: Column<number>
        public compId: Column<string>
        public microHet: ReadonlyMap<number, string[]>

        constructor(public kind: K, private seqIdStart: Column<number>, private seqIdEnd: Column<number>) {
            let maxSeqId = 0, minSeqId = Number.MAX_SAFE_INTEGER;
            for (let i = 0, _i = this.seqIdStart.rowCount; i < _i; i++) {
                const idStart = this.seqIdStart.value(i);
                const idEnd = this.seqIdEnd.value(i);
                if (idStart < minSeqId) minSeqId = idStart;
                if (maxSeqId < idEnd) maxSeqId = idEnd;
            }

            const count = maxSeqId - minSeqId + 1;

            this.code = Column.ofConst('X', count, Column.Schema.str) as Column<Alphabet>;
            this.label = Column.ofConst('', count, Column.Schema.str);
            this.seqId = Column.ofLambda({
                value: row => row + minSeqId + 1,
                rowCount: count,
                schema: Column.Schema.int
            });
            this.compId = Column.ofConst('', count, Column.Schema.str);

            this.offset = minSeqId - 1;
            this.length = count;
        }
    }
}

export { Sequence };
