/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AlignmentOp, globalAlign } from '@molstar/core/algorithms/sequence/needleman-wunsch';
import { categoryColumn } from '@molstar/core/repr/tree';
import { encodeSeq, inlineSeq } from '@molstar/seqx-core/data';
import { type Cif, cifRows } from '../files/cif.js';

// The record rules of p53's structure (examples.md, p53; structure-syncing.md, Datasets), from PDB 2OCJ's mmCIF as
// committed (`fixtures/structures/2OCJ.cif.gz`), by `scripts/record.ts p53-structure`: the entity's letters with chain
// A's secondary structure as Seq*'s `seq` file, and the fake structure viewer's own chain table (every chain instance,
// its residues with their label and author numbers, observed or not), read from the same file and held apart from
// Seq*'s nodes. Here, in the package's sources, so a unit test holds the rules; the script imports them.

/** The entity's system in the `seq` file. */
export const StructureSystem = '2OCJ-1';

/** A residue of a chain instance in the fake viewer's table. */
export interface ChainResidue {
  /** `label_seq_id` (1-based), null for a non-polymer's residue. */
  readonly label: number | null;
  readonly auth: string;
  /** `pdbx_PDB_ins_code`, empty for none. */
  readonly ins: string;
  readonly comp: string;
  readonly observed: boolean;
}

/** A chain instance: `label_asym_id`, `auth_asym_id`, its entity and operators, its residues in order. */
export interface ChainInstance {
  readonly label_asym_id: string;
  readonly auth_asym_id: string;
  readonly entity: string;
  readonly operators: readonly string[];
  readonly residues: readonly ChainResidue[];
}

/** The fake viewer's chain table: the entry's id and its chain instances of the deposited model, waters left out. */
export interface ChainTable {
  readonly entry: string;
  readonly chains: readonly ChainInstance[];
}

const value = (s: string | undefined) => (s === undefined || s === '.' || s === '?' ? '' : s);

/**
 * The chain table: polymer chains from `pdbx_poly_seq_scheme` (every residue of the entity, observed where `atom_site`
 * of the first model has it), non-polymers from `pdbx_nonpoly_scheme` (observed: deposited), waters left out; one
 * operator, the deposited model's (`1_555`).
 */
export function chainTable(cif: Cif, entry: string): ChainTable {
  const water = new Set(
    cifRows(cif, '_entity')
      .filter((e) => e.type === 'water')
      .map((e) => e.id),
  );
  const modeled = new Set<string>();
  for (const a of cifRows(cif, '_atom_site'))
    if (a.pdbx_PDB_model_num === '1' && a.label_seq_id !== '.') modeled.add(`${a.label_asym_id} ${a.label_seq_id}`);
  const chains = new Map<string, { auth: string; entity: string; residues: ChainResidue[] }>();
  const chainOf = (asym: string, auth: string, entity: string) => {
    let c = chains.get(asym);
    if (!c) chains.set(asym, (c = { auth, entity, residues: [] }));
    return c;
  };
  for (const r of cifRows(cif, '_pdbx_poly_seq_scheme'))
    chainOf(r.asym_id, r.pdb_strand_id, r.entity_id).residues.push({
      label: Number(r.seq_id),
      auth: value(r.auth_seq_num) || value(r.pdb_seq_num),
      ins: value(r.pdb_ins_code),
      comp: r.mon_id,
      observed: modeled.has(`${r.asym_id} ${r.seq_id}`),
    });
  for (const r of cifRows(cif, '_pdbx_nonpoly_scheme'))
    if (!water.has(r.entity_id))
      chainOf(r.asym_id, r.pdb_strand_id, r.entity_id).residues.push({
        label: null,
        auth: value(r.auth_seq_num) || value(r.pdb_seq_num),
        ins: value(r.pdb_ins_code),
        comp: r.mon_id,
        observed: true,
      });
  const order = cifRows(cif, '_struct_asym').map((a) => a.id);
  return {
    entry,
    chains: order
      .filter((asym) => chains.has(asym))
      .map((asym) => {
        const c = chains.get(asym)!;
        return {
          label_asym_id: asym,
          auth_asym_id: c.auth,
          entity: c.entity,
          operators: ['1_555'],
          residues: c.residues,
        };
      }),
  };
}

/** The order of a residue row's values in the committed table. */
export const ResidueFields = ['label', 'auth', 'ins', 'comp', 'observed'] as const;

/**
 * The chain table as committed: JSON, a chain instance a line, its residues rows of `ResidueFields` (`observed` 1 or
 * 0), so the file stays small and reads back with `readChainTable`.
 */
export function chainTableJson(t: ChainTable): string {
  const chain = (c: ChainInstance) =>
    JSON.stringify({
      ...c,
      residues: c.residues.map((r) => [r.label, r.auth, r.ins, r.comp, r.observed ? 1 : 0]),
    });
  return `{"entry":${JSON.stringify(t.entry)},"fields":${JSON.stringify(ResidueFields)},"chains":[\n${t.chains.map(chain).join(',\n')}\n]}\n`;
}

/** The committed table read back. */
export function readChainTable(text: string): ChainTable {
  const raw = JSON.parse(text) as {
    entry: string;
    chains: (Omit<ChainInstance, 'residues'> & { residues: [number | null, string, string, string, 0 | 1][] })[];
  };
  return {
    entry: raw.entry,
    chains: raw.chains.map((c) => ({
      ...c,
      residues: c.residues.map(([label, auth, ins, comp, observed]) => ({
        label,
        auth,
        ins,
        comp,
        observed: observed === 1,
      })),
    })),
  };
}

/**
 * The entity's `seq` file: its canonical letters (`entity_poly.pdbx_seq_one_letter_code_can`, the one polymer entity)
 * and chain A's helices and strands (`struct_conf`, `struct_sheet_range`, by `label_seq_id`) as ranges with a
 * `category`.
 */
export function structureSeq(cif: Cif): string {
  const [poly] = cifRows(cif, '_entity_poly');
  const letters = poly.pdbx_seq_one_letter_code_can.replace(/\s+/g, '');
  const rows: { start: number; end: number; kind: string }[] = [];
  for (const r of cifRows(cif, '_struct_conf'))
    if (r.beg_label_asym_id === 'A' && r.conf_type_id.startsWith('HELX'))
      rows.push({ start: Number(r.beg_label_seq_id) - 1, end: Number(r.end_label_seq_id), kind: 'helix' });
  for (const r of cifRows(cif, '_struct_sheet_range'))
    if (r.beg_label_asym_id === 'A')
      rows.push({ start: Number(r.beg_label_seq_id) - 1, end: Number(r.end_label_seq_id), kind: 'strand' });
  rows.sort((a, b) => a.start - b.start || a.end - b.end);
  const data = inlineSeq(
    { name: StructureSystem, length: letters.length, alphabet: 'protein' },
    {
      letters: { form: 'letters', start: 0, letters },
      ranges: {
        form: 'list',
        start: Int32Array.from(rows, (e) => e.start),
        end: Int32Array.from(rows, (e) => e.end),
        id: categoryColumn(rows.map((e) => `${e.kind}${e.start + 1}`)),
        fields: [{ name: 'category', column: categoryColumn(rows.map((e) => e.kind)) }],
      },
    },
  );
  return `${encodeSeq(data)}\n`;
}

const Codes = 'ACDEFGHIKLMNPQRSTVWYX';
const code = (s: string) => Uint8Array.from(s.toUpperCase(), (c) => Math.max(0, Codes.indexOf(c)));

/**
 * The runs of aligned positions of `a` and `b` in their global alignment (core's Needleman–Wunsch: +4 a match, −1 a
 * mismatch, gaps 8 and 1, the ends free): what a `segments` map of the pair holds (`from` on `a`, `to` on `b`).
 */
export function alignedSegments(a: string, b: string): { from: number[]; to: number[]; length: number[] } {
  const { ops } = globalAlign(code(a), code(b), {
    match: 4,
    mismatch: -1,
    gapOpen: 8,
    gapExtend: 1,
    leadingGaps: 'free',
    trailingGaps: 'free',
    gapTies: 'open',
  });
  const out = { from: [] as number[], to: [] as number[], length: [] as number[] };
  let i = 0;
  let j = 0;
  let run = false;
  for (const op of ops) {
    if (op === AlignmentOp.Pair) {
      if (run) out.length[out.length.length - 1]++;
      else {
        out.from.push(i);
        out.to.push(j);
        out.length.push(1);
      }
      run = true;
      i++;
      j++;
    } else {
      run = false;
      if (op === AlignmentOp.OnlyA) i++;
      else j++;
    }
  }
  return out;
}
