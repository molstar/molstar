/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StructureQuery } from '../../mol-model/structure/query';
import { Expression } from '../language/expression';
import { MolScriptBuilder as MS } from '../language/builder';
import { compile } from '../runtime/query/base';

// TODO: make this into a separate "language"?

type ResidueListSelectionEntry =
    | { kind: 'single', asym_id: string; seq_id: number; ins_code?: string }
    | { kind: 'range', asym_id: string; seq_id_beg: number; seq_id_end: number; }

function entriesToQuery(xs: ResidueListSelectionEntry[], kind: 'auth' | 'label') {
    const groups: Expression[] = [];

    const asym_id_key = kind === 'auth' ? 'auth_asym_id' as const : 'label_asym_id' as const;
    const seq_id_key = kind === 'auth' ? 'auth_seq_id' as const : 'label_seq_id' as const;

    for (const x of xs) {
        if (x.kind === 'range') {
            groups.push(MS.struct.generator.atomGroups({
                'chain-test': MS.core.rel.eq([MS.ammp(asym_id_key), x.asym_id]),
                'residue-test': MS.core.rel.inRange([MS.ammp(seq_id_key), x.seq_id_beg, x.seq_id_end])
            }));
        } else {
            const ins_code = (x.ins_code ?? '').trim();

            groups.push(MS.struct.generator.atomGroups({
                'chain-test': MS.core.rel.eq([MS.ammp(asym_id_key), x.asym_id]),
                'residue-test': MS.core.logic.and([
                    MS.core.rel.eq([MS.ammp(seq_id_key), x.seq_id]),
                    MS.core.rel.eq([MS.ammp('pdbx_PDB_ins_code'), ins_code])
                ])
            }));
        }
    }

    const query = MS.struct.combinator.merge(groups);

    return compile(query) as StructureQuery;
}

function parseRange(c: string, s: string[], e: number): ResidueListSelectionEntry | undefined {
    if (!c || s.length === 0 || Number.isNaN(+s[0])) return;
    if (Number.isNaN(e)) {
        return { kind: 'single', asym_id: c, seq_id: +s[0], ins_code: s[1] };
    }
    return { kind: 'range', asym_id: c, seq_id_beg: +s[0], seq_id_end: e };
}

function parseInsCode(e?: string) {
    if (!e) return [];
    return e.split(':');
}

function parseResidueListSelection(input: string): ResidueListSelectionEntry[] {
    return input.split(',') // A 1-3, B 3 => [A 1-3, B 3]
        .map(e => e.trim().split(/\s+|[-]/g).filter(e => !!e)) // [A 1-3, B 3] => [[A, 1, 3], [B, 3]]
        .map(e => parseRange(e[0], parseInsCode(e[1]), +e[2]))
        .filter(e => !!e) as ResidueListSelectionEntry[];
}

// parses a list of residue ranges, e.g. A 10-100, B 30, C 12:i
export function compileResidueListSelection(input: string, idType: 'auth' | 'label') {
    const entries = parseResidueListSelection(input);
    return entriesToQuery(entries, idType);
}