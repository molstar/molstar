/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import {
  type AssetReader,
  DefaultOpenLimits,
  type FormatRuntime,
  type OpenedSource,
  SeqRuntime,
  FastaRuntime,
} from '@molstar/seqx-core/data';
import { AlignmentFileRuntime } from '@molstar/seqx-core/data/alignment-readers';
import { RecordsRuntime } from '@molstar/seqx-core/data/records';
import { createSeqSpec, type SeqSpecBuilder } from '@molstar/seqx-core/lang';
import { MmcifRuntime } from '@molstar/seqx-core/runtimes/mmcif';

// A user's own file dropped onto a view (examples.md, Showcases; the formats built): the file's name says its format,
// the runtime of the format opens it (the same runtimes a spec's source node is read by) and its coverage says what the
// file holds: sequences (`seq`, FASTA), the units of a structure (mmCIF), an alignment (Stockholm, aligned FASTA) or rows
// of ranges (the records: a table with `start` and `end` columns). The spec the view then shows names the file by a
// location the application's assets port answers from the dropped bytes; nothing leaves the page. GenBank is not read
// (decided: the fifth slice at the earliest) and says so.

/** The runtimes a viewer that takes dropped files lists (the default registry has `seq` and the records). */
export const DropRuntimes: readonly FormatRuntime[] = [
  SeqRuntime,
  RecordsRuntime,
  FastaRuntime,
  AlignmentFileRuntime,
  MmcifRuntime,
];

/** The most sequences or units a dropped file's view shows: one section each. */
export const MaxSections = 8;

/** The alphabet a source names, narrowed to the ones a node takes (protein when none). */
const alphabetOf = (a: string | undefined): 'dna' | 'rna' | 'protein' => (a === 'dna' || a === 'rna' ? a : 'protein');

/** The location the assets port answers a dropped file at. */
export const droppedUri = (name: string) => `dropped/${encodeURIComponent(name)}`;

const Extensions: ReadonlyArray<readonly [RegExp, string]> = [
  [/\.seq$/i, 'seq'],
  [/\.(fa|fasta|faa|fna)$/i, 'fasta'],
  [/\.(sto|stk|stockholm)$/i, 'stockholm'],
  [/\.(afa|aln|aligned\.fasta)$/i, 'aligned-fasta'],
  [/\.tsv$/i, 'tsv'],
  [/\.csv$/i, 'csv'],
  [/\.json$/i, 'json'],
  [/\.(msgpack|mpk)$/i, 'msgpack'],
  [/\.(gff|gff3)$/i, 'gff3'],
  [/\.(cif|mmcif)(\.gz)?$/i, 'mmcif'],
];

/** The format a file name says, or why it is refused. */
export function formatOfName(name: string): { readonly format: string } | { readonly refused: string } {
  if (/\.(gb|gbk|genbank)(\.gz)?$/i.test(name)) return { refused: 'GenBank files are not read yet' };
  for (const [re, format] of Extensions) if (re.test(name)) return { format };
  return {
    refused: `The format of ${name} is not known: Seq* reads .seq, FASTA, Stockholm, aligned FASTA, TSV, CSV, JSON, MessagePack, GFF3 and mmCIF files`,
  };
}

/** What a dropped file made of it: the spec to show, and a line saying what the file held. */
export type Dropped =
  | { readonly spec: Node; readonly format: string; readonly summary: string }
  | { readonly refused: string };

const readerOf = (bytes: Uint8Array): AssetReader => ({
  readText: async () => new TextDecoder().decode(bytes),
  readBytes: async () => bytes,
});

/** The opened source of the dropped bytes, or the refusal (a file that does not parse says why). */
async function openDropped(
  name: string,
  format: string,
  bytes: Uint8Array,
): Promise<OpenedSource | { refused: string }> {
  const runtime = DropRuntimes.find((r) => r.formats.includes(format));
  if (!runtime) return { refused: `No runtime reads ${format}` };
  try {
    return await runtime.open(
      { address: '/data/file', format, uri: droppedUri(name) },
      { assets: readerOf(bytes), signal: new AbortController().signal, limits: DefaultOpenLimits },
    );
  } catch (e) {
    const issues = (e as { issues?: readonly { message: string }[] }).issues;
    return {
      refused: `${name} is not read as ${format}: ${issues?.[0]?.message ?? (e instanceof Error ? e.message : e)}`,
    };
  }
}

const isRefusal = (r: OpenedSource | { refused: string }): r is { refused: string } => 'refused' in r;

/** One section a system: its ruler, its letters and, where the source has them, its values and ranges. */
function sequenceSections(b: SeqSpecBuilder, opened: OpenedSource, format: string, name: string) {
  const source = '/data/file' as never;
  const systems = opened.coverage.systems.slice(0, MaxSections);
  const withRanges = new Set(opened.coverage.shapes.filter((s) => s.shape === 'ranges').map((s) => s.system));
  const withValues = new Set(opened.coverage.shapes.filter((s) => s.shape === 'values').map((s) => s.system));
  const parts = systems.map((s, i) => {
    const id = `s${i + 1}`;
    const sequence = b.data.sequence(id, { alphabet: alphabetOf(s.alphabet), source, name: s.name });
    const ranges = withRanges.has(s.name) ? b.data.ranges(`${id}-ranges`, { target: sequence, source }) : undefined;
    const values = withValues.has(s.name) ? b.data.values(`${id}-values`, { target: sequence, source }) : undefined;
    return { id, label: s.name, sequence, ranges, values };
  });
  b.views.view('main', (v) => {
    for (const p of parts)
      v.section(p.id, { of: p.sequence }, (s) => {
        s.track('ruler', { label: p.label }).layer('ruler', { type: 'ruler' });
        // a structure's units carry values (its B-factors, a prediction's pLDDT), never letters
        if (format !== 'mmcif') s.track('letters', { label: 'Sequence' }).layer('letters', { type: 'letters' });
        if (p.values) s.track('values', { label: 'Values' }).layer('values', { type: 'values', data: p.values });
        if (p.ranges) s.track('ranges', { label: 'Ranges' }).layer('ranges', { type: 'ranges', data: p.ranges });
      });
  });
  const more = opened.coverage.systems.length - systems.length;
  return `${format}: ${systems.length} ${format === 'mmcif' ? 'chains' : 'sequences'}${more > 0 ? ` of ${opened.coverage.systems.length} shown` : ''}`;
}

/**
 * The spec a dropped file is shown by, and what it held: the file is opened with the runtime of its format (so a file
 * that does not parse is refused here, with the runtime's reason) and its coverage decides the sections.
 */
export async function specOfDropped(name: string, bytes: Uint8Array): Promise<Dropped> {
  const named = formatOfName(name);
  if ('refused' in named) return named;
  const { format } = named;
  const opened = await openDropped(name, format, bytes);
  if (isRefusal(opened)) return opened;
  const b = createSeqSpec({ title: name });
  const common = { format, uri: droppedUri(name) } as const;
  if (format === 'stockholm' || format === 'aligned-fasta') {
    const file = b.data.source('file', common);
    const alphabet = alphabetOf(opened.coverage.systems.find((s) => s.alphabet)?.alphabet);
    const aln = b.data.alignment('aln', { alphabet, source: file });
    b.views.view('main', (v) =>
      v.section('alignment', { of: aln }, (s) => {
        s.track('logo', { label: 'Conservation', height: 60 }).layer('logo', { type: 'logo', data: aln });
        s.track('rows', { label: '' }).layer('rows', { type: 'rows', data: aln });
      }),
    );
    dispose(opened);
    return { spec: b.build(), format, summary: `${format}: an alignment` };
  }
  if (format === 'seq' || format === 'fasta' || format === 'mmcif') {
    b.data.source('file', common);
    const summary = sequenceSections(b, opened, format, name);
    if (!opened.coverage.systems.length) {
      dispose(opened);
      return { refused: `${name} holds no sequence to show` };
    }
    dispose(opened);
    return { spec: b.build(), format, summary };
  }
  // the records: rows of ranges (`start`, `end`, and what the file's defaults read) on a sequence the rows span
  const read = await opened.read(
    { system: opened.coverage.systems[0]?.name ?? '', start: 0, end: 2 ** 30, level: 1, shapes: ['ranges'] },
    new AbortController().signal,
  );
  dispose(opened);
  const ends = read.ok && read.tables.ranges?.form === 'list' ? read.tables.ranges.end : undefined;
  if (!ends?.length) return { refused: `${name} holds no rows of ranges (a table with start and end columns)` };
  const length = Math.max(...ends);
  const system = b.data.source('file', { ...common, system: 'target' });
  const target = b.data.sequence('target', { alphabet: 'protein', letters: 'X'.repeat(length) });
  const ranges = b.data.ranges('rows', { target, source: system });
  b.views.view('main', (v) =>
    v.section('target', { of: target }, (s) => {
      s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
      s.track('ranges', { label: name }).layer('ranges', { type: 'ranges', data: ranges });
    }),
  );
  return { spec: b.build(), format, summary: `${format}: ${ends.length} rows on positions 1 to ${length}` };
}

function dispose(opened: OpenedSource) {
  try {
    opened.dispose();
  } catch {
    // a runtime's own cleanup is not the drop's to report
  }
}
