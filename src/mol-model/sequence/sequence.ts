/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { AminoAlphabet, NuclecicAlphabet, getProteinOneLetterCode, getRnaOneLetterCode, getDnaOneLetterCode } from './constants';
import { Column } from '../../mol-data/db';

// TODO add mapping support to other sequence spaces, e.g. uniprot

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

        /** One letter code */
        readonly code: Column<Alphabet>
        readonly label: Column<string>

        readonly seqId: Column<number>
        /** Component id */
        readonly compId: Column<string>

        /** returns index for given seqId */
        readonly index: (seqId: number) => number

        /** maps seqId to list of compIds */
        readonly microHet: ReadonlyMap<number, string[]>
    }

    export interface Protein extends Base<Kind.Protein, AminoAlphabet> { }
    export interface RNA extends Base<Kind.RNA, NuclecicAlphabet> { }
    export interface DNA extends Base<Kind.DNA, NuclecicAlphabet> { }
    export interface Generic extends Base<Kind.Generic, 'X' | '-'> { }

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
        public length: number
        public code: Column<Alphabet>
        public label: Column<string>
        public seqId: Column<number>
        public compId: Column<string>
        public microHet: ReadonlyMap<number, string[]> = new Map()

        private indexMap: Map<number, number>
        index(seqId: number) {
            return this.indexMap.get(seqId)!;
        }

        constructor(public kind: K, compId: Column<string>, seqId: Column<number>) {
            const codeFromName = codeProvider(kind);
            const codes: string[] = [];
            const compIds: string[] = [];
            const seqIds: number[] = [];
            const microHet = new Map<number, string[]>();

            let idx = 0;
            const indexMap = new Map<number, number>();
            for (let i = 0, il = seqId.rowCount; i < il; ++i) {
                const seq_id = seqId.value(i);

                if (!indexMap.has(seq_id)) {
                    indexMap.set(seq_id, idx);
                    const comp_id = compId.value(i);
                    compIds[idx] = comp_id;
                    seqIds[idx] = seq_id;
                    codes[idx] = codeFromName(comp_id);
                    idx += 1;
                } else {
                    // micro-heterogeneity
                    if (!microHet.has(seq_id)) {
                        microHet.set(seq_id, [compIds[indexMap.get(seq_id)!], compId.value(i)]);
                    } else {
                        microHet.get(seq_id)!.push(compId.value(i));
                    }
                }
            }

            const labels: string[] = [];
            for (let i = 0, il = idx; i < il; ++i) {
                const mh = microHet.get(seqIds[i]);
                if (mh) {
                    const l = mh.map(id => {
                        const c = codeFromName(id);
                        return c === 'X' ? id : c;
                    });
                    labels[i] = `(${l.join('|')})`;
                } else {
                    labels[i] = codes[i] === 'X' ? compIds[idx] : codes[i];
                }
            }

            this.length = idx;
            this.code = Column.ofStringArray(codes) as Column<Alphabet>;
            this.compId = Column.ofStringArray(compIds);
            this.seqId = Column.ofIntArray(seqIds);
            this.label = Column.ofStringArray(labels);
            this.microHet = microHet;
            this.indexMap = indexMap;
        }
    }

    export function ofSequenceRanges(seqIdBegin: Column<number>, seqIdEnd: Column<number>): Sequence {
        const kind = Kind.Generic;

        return new SequenceRangesImpl(kind, seqIdBegin, seqIdEnd) as Sequence;
    }

    class SequenceRangesImpl<K extends Kind, Alphabet extends string> implements Base<K, Alphabet> {
        public length: number
        public code: Column<Alphabet>
        public label: Column<string>
        public seqId: Column<number>
        public compId: Column<string>
        public microHet: ReadonlyMap<number, string[]> = new Map()

        private minSeqId: number
        index(seqId: number) {
            return seqId - this.minSeqId;
        }

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

            this.length = count;
            this.minSeqId = minSeqId;
        }
    }
}

export { Sequence };
