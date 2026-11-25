/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import { encodeSeq, inlineSeq, type SeqData } from '@molstar/seqx-core/data';
import { type Cif, cifRows } from '../files/cif.js';

// What the record script keeps of PDB 1EHZ's mmCIF (`scripts/record.ts trna`): chain A's sequence, its modified residues
// and the entry's own base-pair table (`ndb_struct_na_base_pair`), nothing renamed. Kept here, in the package's
// sources, so a unit test holds the rules; the script imports them.

/**
 * The six tertiary pairs of the structure (1-based): those outside the cloverleaf's four stems that hold the L-shape
 * (Kim et al. 1974). Four of them cross the stems' pairs, 26–44 crosses none (it encloses the anticodon stem), 54–58
 * crosses only tertiary pairs.
 */
export const TertiaryPairs: readonly (readonly [number, number])[] = [
  [8, 14],
  [15, 48],
  [18, 55],
  [19, 56],
  [26, 44],
  [54, 58],
];

/** A pair's kind by its Saenger class (the entry's `hbond_type_28`): 19 and 20 Watson–Crick, 28 wobble, else non-canonical. */
export function kindOfSaenger(c: string): 'canonical' | 'wobble' | 'non-canonical' {
  return c === '19' || c === '20' ? 'canonical' : c === '28' ? 'wobble' : 'non-canonical';
}

export interface Modified {
  /** 0-based. */
  readonly position: number;
  readonly code: string;
  readonly parent: string;
}

export interface EntryPair {
  /** 0-based, `i` before `j` in the chain. */
  readonly i: number;
  readonly j: number;
  readonly kind: 'canonical' | 'wobble' | 'non-canonical';
  readonly saenger: string;
  readonly tertiary: boolean;
}

export interface Trimmed {
  readonly letters: string;
  readonly modified: readonly Modified[];
  readonly pairs: readonly EntryPair[];
}

/** Chain A's canonical sequence (parent letters), its modified residues and its base pairs, from the entry's text. */
export function trimEntry(cif: Cif, chain = 'A'): Trimmed {
  const poly = cifRows(cif, '_entity_poly')[0];
  const letters = poly.pdbx_seq_one_letter_code_can.replaceAll(/\s/g, '');
  const modified = cifRows(cif, '_pdbx_struct_mod_residue')
    .filter((r) => r.label_asym_id === chain)
    .map((r) => ({ position: Number(r.label_seq_id) - 1, code: r.label_comp_id, parent: r.parent_comp_id }))
    .sort((a, b) => a.position - b.position);
  const pairs = cifRows(cif, '_ndb_struct_na_base_pair')
    .filter((r) => r.i_label_asym_id === chain && r.j_label_asym_id === chain)
    .map((r): EntryPair => {
      const a = Number(r.i_label_seq_id);
      const b = Number(r.j_label_seq_id);
      const [i, j] = a < b ? [a, b] : [b, a];
      return {
        i: i - 1,
        j: j - 1,
        kind: kindOfSaenger(r.hbond_type_28),
        saenger: r.hbond_type_28,
        tertiary: TertiaryPairs.some(([x, y]) => x === i && y === j),
      };
    })
    .sort((a, b) => a.i - b.i);
  return { letters, modified, pairs };
}

export const SystemName = '1EHZ-A';

/** The molecule's `seq` file: its letters and the modified residues (a range each, its code as `label`, its parent letter as `category`). */
export function moleculeSeq(t: Trimmed): string {
  const data: SeqData = inlineSeq(
    { name: SystemName, length: t.letters.length, alphabet: 'rna' },
    {
      letters: { form: 'letters', start: 0, letters: t.letters },
      ranges: {
        form: 'list',
        start: Int32Array.from(t.modified, (m) => m.position),
        end: Int32Array.from(t.modified, (m) => m.position + 1),
        id: categoryColumn(t.modified.map((m) => `${m.code}${m.position + 1}`)),
        fields: [
          { name: 'label', column: categoryColumn(t.modified.map((m) => m.code)) },
          { name: 'category', column: categoryColumn(t.modified.map((m) => m.parent)) },
        ],
      },
    },
  );
  return `${encodeSeq(data)}\n`;
}

/** The base pairs' `seq` file: a range on the first base, its partner's position and the pair's kind, class and flag. */
export function pairsSeq(t: Trimmed): string {
  const data: SeqData = inlineSeq(
    { name: SystemName, length: t.letters.length, alphabet: 'rna' },
    {
      ranges: {
        form: 'list',
        start: Int32Array.from(t.pairs, (p) => p.i),
        end: Int32Array.from(t.pairs, (p) => p.i + 1),
        id: categoryColumn(t.pairs.map((p) => `${p.i + 1}-${p.j + 1}`)),
        fields: [
          { name: 'partner', column: Int32Array.from(t.pairs, (p) => p.j) },
          { name: 'type', column: categoryColumn(t.pairs.map((p) => p.kind)) },
          { name: 'saenger', column: categoryColumn(t.pairs.map((p) => p.saenger)) },
          { name: 'tertiary', column: categoryColumn(t.pairs.map((p) => (p.tertiary ? 'yes' : 'no'))) },
        ],
      },
    },
  );
  return `${encodeSeq(data)}\n`;
}
