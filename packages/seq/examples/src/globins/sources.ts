/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { createSeqSpec } from '@molstar/seqx-core/lang';
import { rowsSystemOf } from '@molstar/seqx-core/resolve';
import { ColumnsFile, ReferenceRow, SeedFile } from './data.js';
import { SeedAttribution, UniprotAttribution } from './recording.js';

export { SeedAttribution, UniprotAttribution };

// The family's sources, one helper for every view of the seed: the Stockholm file and the rows' columns file joined
// by row key, each with its attribution and licence (carried into an exported spec). CC BY 4.0 asks that changes be
// indicated: the UniProt file's attribution says what was derived.

/** The seed (`/data/seed`) as an alignment (`/data/globins`, reference row set) and its columns file (`/data/columns`). */
export function familySources(b: ReturnType<typeof createSeqSpec>) {
  const file = b.data.source('seed', {
    format: 'stockholm',
    uri: SeedFile,
    attribution: SeedAttribution,
    licence: 'CC0-1.0',
  });
  const aln = b.data.alignment('globins', { alphabet: 'protein', source: file, reference: ReferenceRow });
  // a records file whose system is the rows system: its fields are joined to the rows by their key (`id`)
  b.data.source('columns', {
    format: 'tsv',
    uri: ColumnsFile,
    system: rowsSystemOf('/data/globins'),
    types: { uniprot_length: 'int32' },
    attribution: UniprotAttribution,
    licence: 'CC-BY-4.0',
  });
  return { file, aln };
}
