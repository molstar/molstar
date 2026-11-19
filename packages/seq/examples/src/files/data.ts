/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { encodeMsgPack } from '@molstar/core/data/msgpack';
import { Domains } from '../p53/data.js';
import { PlasmidFeatures, PlasmidLength, plasmidLetters, RotatedOrigin, rotate } from '../plasmid/data.js';

// The Files dataset: p53's domains as a table in every file form the records runtime reads (TSV, CSV, JSON and
// MessagePack, an array of records and an object of columns), and the plasmid's features as GFF3. The files are
// generated from the code's own data and committed (`fixtures/files/`), compared with their generators by a test
// (`SEQ_RECORDS=update` writes them again); the views read the committed files by location, as a published file is.
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

/** The sequence's name in the plasmid's GFF3 (its `seqid`) and the rotated construct's origin. */
export const PlasmidSeqId = 'pUC19c';

/** SO terms for the stand-in's feature kinds. */
const SoTerms: Readonly<Record<string, string>> = {
  CDS: 'CDS',
  misc: 'misc_feature',
  origin: 'rep_origin',
  promoter: 'promoter',
};

/** A feature of the rotated construct: where it starts on it (0-based) and how long it is. */
function onRotated(start: number, end: number) {
  const s = rotate(start, RotatedOrigin);
  return { s, length: end - start };
}

/**
 * The plasmid's features on the rotated construct as NCBI writes a ring's: GFF3 with one `region` row (`Is_circular`),
 * a feature that crosses the origin ending past the length (one row, not two).
 */
export function plasmidGff3(): string {
  const lines = ['##gff-version 3', `##sequence-region ${PlasmidSeqId} 1 ${PlasmidLength}`];
  const row = (type: string, s: number, e: number, strand: string, attributes: string) =>
    [PlasmidSeqId, 'seqx', type, s, e, '.', strand, '.', attributes].join('\t');
  lines.push(row('region', 1, PlasmidLength, '+', `ID=region-${PlasmidSeqId};Name=${PlasmidSeqId};Is_circular=true`));
  for (const f of PlasmidFeatures) {
    const { s, length } = onRotated(f.start, f.end);
    lines.push(
      row(
        SoTerms[f.category] ?? f.category,
        s + 1,
        s + length,
        f.strand,
        `ID=${f.id};Name=${f.label};product=${f.label.replaceAll(',', '%2C')}`,
      ),
    );
  }
  // the primer site across the rotated origin, on the plus strand (generated, EC-11)
  const p = onRotated(RotatedOrigin - 30, RotatedOrigin + 25);
  lines.push(row('primer_bind', p.s + 1, p.s + p.length, '+', 'ID=primer;Name=primer site;note=across the origin'));
  // what follows `##FASTA` is the file's sequence, which the reader ignores (the view's sequence is its own)
  const letters = plasmidLetters();
  lines.push(
    '##FASTA',
    `>${PlasmidSeqId}`,
    (letters.slice(RotatedOrigin) + letters.slice(0, RotatedOrigin)).slice(0, 60),
  );
  return `${lines.join('\n')}\n`;
}

/** The GFF3 file's location, and the table files' folder. */
export const PlasmidGffFile = 'fixtures/files/plasmid.gff3';
export const FilesFolder = 'fixtures/files/';
