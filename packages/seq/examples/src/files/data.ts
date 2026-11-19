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

// One table, two sequences (files.md, "Its systems, by name"): a `target` column names each row's sequence, an empty
// target is the source's own system. The rows are generated here (p53's domains from the Domains data, CDK2's loops as
// illustrative bounds) and committed like the domains table.

/** The many-sequences table's file name (under the tables' folder) and its columns. */
export const ManySequencesFile = 'many-sequences.tsv';
export const ManyColumns = ['target', 'id', 'label', 'start', 'end', 'category'] as const;

/** CDK2 (UniProt P24941, 298 residues; the Complex dataset's chain). */
export const Cdk2 =
  'MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL';

/** The source's own system (a row with no target is on it) and the two named ones. */
export const ManySystems = { own: 'p53', p53: 'p53', cdk2: 'cdk2' } as const;

/**
 * The rows, in file order, positions from 1 with the last one named: p53's TAD, DBD and TET named by their target, its
 * NLS with an empty one (the source's own system, p53), and three CDK2 elements (illustrative bounds).
 */
export function manyRows() {
  const p53 = (id: string, target: string) => {
    const i = Domains.id.indexOf(id);
    return {
      target,
      id,
      label: Domains.label[i],
      start: Domains.start[i] + 1,
      end: Domains.end[i],
      category: Domains.category[i],
    };
  };
  return [
    p53('TAD', ManySystems.p53),
    { target: ManySystems.cdk2, id: 'GRL', label: 'Gly-rich loop', start: 12, end: 19, category: 'loop' },
    p53('DBD', ManySystems.p53),
    { target: ManySystems.cdk2, id: 'PSTAIRE', label: 'PSTAIRE helix', start: 46, end: 56, category: 'helix' },
    p53('NLS', ''),
    { target: ManySystems.cdk2, id: 'TLOOP', label: 'T-loop', start: 145, end: 172, category: 'loop' },
    p53('TET', ManySystems.p53),
  ];
}

/** The many-sequences table as TSV. */
export function manyContent(): string {
  const lines = [ManyColumns.join('\t')];
  for (const r of manyRows()) lines.push(ManyColumns.map((c) => cell(r[c], '\t')).join('\t'));
  return `${lines.join('\n')}\n`;
}
