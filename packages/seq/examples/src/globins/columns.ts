/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rowsSystemOf } from '@molstar/seqx-core/resolve';
import { ColumnsFile, familyRegistry, ReferenceRow, SeedFile } from './data.js';

// The rows' columns (Features/Rows): the Globins seed with data columns on both sides of its rows, one of each kind:
// its species and its class (UniProt's, from a file joined by row key) as text and as a chip, its identity to the
// reference row as a number, its residues as a bar and its UniProt entry's whole length as a number. A title's click
// orders the rows by its field; a narrow view hides the right ones first.

export { ColumnsFile };

/** What the view needs: the Stockholm and records runtimes (the default viewer has the latter; a bare context neither). */
export const registry = familyRegistry;

/** The family with a column of each kind on each side, and a file joined by row key. */
export function build(): Node {
  const b = createSeqSpec({ title: 'Globins: row columns' });
  const file = b.data.source('seed', { format: 'stockholm', uri: SeedFile });
  const aln = b.data.alignment('globins', { alphabet: 'protein', source: file, reference: ReferenceRow });
  // a records file whose system is the rows system: its fields are joined to the rows by their key (`id`)
  b.data.source('columns', {
    format: 'tsv',
    uri: ColumnsFile,
    system: rowsSystemOf('/data/globins'),
    types: { uniprot_length: 'int32' },
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
        { field: 'uniprot_length', side: 'right', as: 'number', title: 'UniProt' },
      ],
    });
  });
  return b.build();
}
