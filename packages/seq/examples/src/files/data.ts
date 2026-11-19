/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { encodeMsgPack } from '@molstar/core/data/msgpack';
import { Domains } from '../p53/data.js';
import { PlasmidFiles, PlasmidName } from '../plasmid/derive.js';

// The Files dataset: p53's domains as a table in every file form the records runtime reads (TSV, CSV, JSON and
// MessagePack, an array of records and an object of columns), and the plasmid's features as GFF3 (the Plasmid
// dataset's file). The tables are generated from the code's own data and committed (`fixtures/files/`), compared with
// their generators by a test (`SEQ_RECORDS=update` writes them again); the views read the committed files by location,
// as a published file is.
// The tables count from 1 and name their last position, as a file does (the runtime's defaults), so no spec needs a
// transform to read them.

/** The columns of the domains table, in file order. */
export const DomainColumns = ['id', 'label', 'start', 'end', 'category'] as const;

/** The domains as a table: positions from 1, the last one named (the file's convention). */
export function domainRows() {
  return Domains.id.map((id, i) => ({
    id,
    label: Domains.label[i],
    start: Domains.start[i] + 1,
    end: Domains.end[i],
    category: Domains.category[i],
  }));
}

/** The domains as columns (an object of equal-length arrays). */
export function domainColumns() {
  const rows = domainRows();
  return Object.fromEntries(DomainColumns.map((c) => [c, rows.map((r) => r[c])]));
}

const cell = (v: string | number, delimiter: string) => {
  const s = String(v);
  return s.includes(delimiter) || s.includes('"') || s.includes('\n') ? `"${s.replaceAll('"', '""')}"` : s;
};

function delimited(delimiter: string): string {
  const lines = [DomainColumns.join(delimiter)];
  for (const r of domainRows()) lines.push(DomainColumns.map((c) => cell(r[c], delimiter)).join(delimiter));
  return `${lines.join('\n')}\n`;
}

/**
 * The columns as JSON, a column a line: as the repository's formatter writes a short array (the files are checked in
 * as it formats them, so a format check leaves them as they are).
 */
function soaJson(): string {
  const lines = Object.entries(domainColumns()).map(
    ([name, column]) => `  "${name}": [${column.map((v) => JSON.stringify(v)).join(', ')}]`,
  );
  return `{\n${lines.join(',\n')}\n}\n`;
}

/** A fixture of the domains table by file name under `fixtures/files/`: its form and its content. */
export const DomainFiles = {
  'domains.tsv': { format: 'tsv', content: () => delimited('\t') },
  'domains.csv': { format: 'csv', content: () => delimited(',') },
  'domains-aos.json': { format: 'json', content: () => `${JSON.stringify(domainRows(), null, 2)}\n` },
  'domains-soa.json': { format: 'json', content: soaJson },
  'domains-aos.msgpack': { format: 'msgpack', content: () => encodeMsgPack(domainRows()) },
  'domains-soa.msgpack': { format: 'msgpack', content: () => encodeMsgPack(domainColumns()) },
} as const satisfies Record<string, { format: string; content: () => string | Uint8Array }>;

export type DomainFile = keyof typeof DomainFiles;

/**
 * The plasmid's GFF3 (the Plasmid dataset's rotated variant: its features derived in code, a region row saying it is
 * circular, bla and a primer site ending past the length), the sequence's name in it (its `seqid`), and the table files'
 * folder.
 */
export const PlasmidGffFile: string = PlasmidFiles.rotatedFeatures;
export const PlasmidSeqId: string = PlasmidName;
export const FilesFolder = 'fixtures/files/';
