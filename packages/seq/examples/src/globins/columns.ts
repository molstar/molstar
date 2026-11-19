/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rowsSystemOf } from '@molstar/seqx-core/resolve';
import { FamilySpecies, globinFamily } from './data.js';

// The rows' columns (Features/Rows): the Globins family with data columns on both sides of its rows, one of each kind:
// its species (the rows' own `OS` field) as text and its class as a chip from a file joined by row key, its identity
// to the ancestor as a number, its residues as a bar and its score from the joined file as a number. A title's click
// orders the rows by its field; a narrow view hides the right ones first.

/** The joined file: a record a row, keyed by the row's label (`id`), its class and a score. */
export const ColumnsFile = 'fixtures/globins/family-columns.tsv';

/** Each species' class (what the chip colors). */
const ClassOf: Readonly<Record<(typeof FamilySpecies)[number], string>> = {
  'Homo sapiens': 'Mammalia',
  'Pan troglodytes': 'Mammalia',
  'Mus musculus': 'Mammalia',
  'Bos taurus': 'Mammalia',
  'Gallus gallus': 'Aves',
  'Xenopus laevis': 'Amphibia',
  'Danio rerio': 'Actinopterygii',
  'Physeter catodon': 'Mammalia',
};

/**
 * The joined file's text (what `fixtures/globins/family-columns.tsv` holds, compared by a test): `id`, `class` and a
 * `score` (the row's residues less its index, a stand-in for a search's bit score); every fifth row left out, so a
 * missing record is an empty cell.
 */
export function columnsTsv(): string {
  const f = globinFamily();
  const lines = ['id\tclass\tscore'];
  f.label.forEach((label, k) => {
    if (k % 5 === 4) return;
    const residues = f.letters[k].replace(/-/g, '').length;
    lines.push(`${label}\t${ClassOf[f.species[k] as (typeof FamilySpecies)[number]]}\t${residues * 2 - k}`);
  });
  return `${lines.join('\n')}\n`;
}

/** What the view needs: the records runtime (the default viewer has it; a bare context does not). */
export function registry(): readonly SeqRegistryEntry[] {
  return [{ kind: 'format-runtime', runtime: RecordsRuntime }];
}

/** The family with a column of each kind on each side, and a file joined by row key. */
export function build(): Node {
  const b = createSeqSpec({ title: 'Globins: row columns' });
  const f = globinFamily();
  const aln = b.data.alignment('globins', {
    alphabet: 'protein',
    // an inline table's other columns are the rows' own fields
    rows: { label: [...f.label], letters: [...f.letters], OS: [...f.species] } as {
      label: string[];
      letters: string[];
    },
    reference: f.label[0],
  });
  // a records file whose system is the rows system: its fields are joined to the rows by their key (`id`)
  b.data.source('columns', {
    format: 'tsv',
    uri: ColumnsFile,
    system: rowsSystemOf('/data/globins'),
    types: { score: 'int32' },
  });
  b.views.view('main', (v) => {
    const s = v.section('family', { of: aln });
    s.track('ruler', { label: 'Columns' }).layer('ruler', { type: 'ruler' });
    s.track('rows', { label: '' }).layer('rows', {
      type: 'rows',
      data: aln,
      columns: [
        { field: 'OS', title: 'Species', width: 104 },
        { field: 'class', as: 'chip', title: 'Class', width: 92 },
        { field: 'identity', side: 'right', as: 'number', title: 'Identity' },
        { field: 'length', side: 'right', as: 'bar', title: 'Residues' },
        { field: 'score', side: 'right', as: 'number', title: 'Score' },
      ],
    });
  });
  return b.build();
}
