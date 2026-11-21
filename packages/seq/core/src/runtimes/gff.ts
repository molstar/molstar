/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { SourceIssues } from '../data/issues.js';
import {
  abortError,
  type Coverage,
  type FailureKind,
  type FormatRuntime,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type SourceEntry,
  type SystemInfo,
} from '../data/protocol.js';
import { wholeLevel } from '../data/seq-format.js';
import { SourceOpenError } from '../data/seq-runtime.js';
import {
  type Assembly,
  assembleTranscripts,
  DefaultTranscriptTypes,
  type GffItem,
  mergeIntervals,
  transcriptTest,
} from './gff-assemble.js';
import { GffCauses } from './gff-issues.js';
import { type GffFormat, type GffRow, MaxAttributes, parseGff } from './gff-parse.js';
import { fieldInfos, type ItemField, NameIndex, rangesOf, rowsIn, type SystemItems, systemItems } from './gff-read.js';
import { decodeText, textBytes } from './internal/text.js';

// GFF3 and GTF (_spec/README.md "Data", GFF3 and GTF), the opt-in runtime `@molstar/seqx-core/runtimes/gff`: a file read
// whole within `limits` (plain, gzip or BGZF by its magic), its lines assembled into transcripts by their type
// (`assembleTranscripts`), served as `ranges` by the source's `features` (`transcripts`, `genes` or `rows`), each item's
// pieces with its whole extent; `##sequence-region` lengths (EC-188); a name index over GFF3's `Name`, `ID`, `Alias` or
// GTF's `gene_name`, `gene_id`, `transcript_name`, `transcript_id` (EC-228). Positions as the formats write them (base
// 1, ends inclusive), converted at the edge: the source's `base` and `ends` are not read (EC-195). By region (a tabix
// `index`) is to come: an `index` is said and the file read whole.

export {
  type AssembleOptions,
  type Assembly,
  assembleTranscripts,
  DefaultTranscriptTypes,
  type GffItem,
  PieceTypes,
} from './gff-assemble.js';
export { GffIssues, GtfIssues } from './gff-issues.js';
export type { GffFormat, GffRow } from './gff-parse.js';
export { parseGff } from './gff-parse.js';

export const GffFormats = ['gff3', 'gtf'] as const;

/** What a source's items are: its transcripts (the default), its genes (transcripts merged) or the file's lines. */
export type GffFeatures = 'transcripts' | 'genes' | 'rows';

export interface GffOptions {
  /** The types read as transcripts (any type ending in `_transcript` too). */
  readonly transcriptTypes: readonly string[];
}

export const DefaultGffOptions: GffOptions = { transcriptTypes: DefaultTranscriptTypes };

/** A whole-sequence feature's types: a feature of one over a whole sequence describes it (Ensembl writes one each). */
const Describing = new Set(['region', 'chromosome', 'supercontig']);

/** The fields `findNames` searches, by format. */
const NameFields: Readonly<Record<GffFormat, readonly string[]>> = {
  gff3: ['Name', 'ID', 'Alias'],
  gtf: ['gene_name', 'gene_id', 'transcript_name', 'transcript_id'],
};

/** Attributes not served as fields: the identity and the hierarchy (`id` and the item's pieces say them). */
const NotFields = new Set(['ID', 'Parent']);

const IndexHint = 'A smaller file, or a host’s adapter reading it by region';

function refuse(
  address: string,
  code: 'seq.source.unreadable' | 'seq.source.malformed' | 'seq.source.too_large' | 'seq.source.no_data',
  message: string,
  hint?: string,
): never {
  throw new SourceOpenError([issueOf(SourceIssues, code, message, { address, ...(hint && { hint }) })]);
}

function failure(kind: FailureKind, message: string): RegionAnswer & { ok: false } {
  const code =
    kind === 'not-served'
      ? 'seq.source.not_served'
      : kind === 'too-large'
        ? 'seq.source.too_large'
        : 'seq.source.failed';
  return { ok: false, kind, issue: issueOf(SourceIssues, code, message) };
}

/** An item's label: GFF3's `Name`, GTF's `transcript_name` or `gene_name`, else its id. */
function labelOf(it: GffItem): string {
  const a = it.attributes;
  return a.Name ?? (it.role === 'gene' ? a.gene_name : (a.transcript_name ?? a.gene_name)) ?? it.id;
}

/** The items of the file's lines (`features: rows`): rows sharing an `ID` one item, its pieces those rows. */
function rowItems(rows: readonly GffRow[]): (GffItem & { readonly row: GffRow })[] {
  const out: (GffItem & { readonly row: GffRow })[] = [];
  const byId = new Map<string, number>();
  const seen = new Map<string, number>();
  for (const r of rows) {
    if (r.id !== undefined) {
      const k = byId.get(r.id);
      if (k !== undefined) {
        const it = out[k];
        out[k] = {
          ...it,
          start: Math.min(it.start, r.start),
          end: Math.max(it.end, r.end),
          pieces: mergeIntervals([...it.pieces, r.start, r.end]),
        };
        continue;
      }
      byId.set(r.id, out.length);
    }
    const base = r.id ?? `${r.type}:${r.start + 1}-${r.end}`;
    const k = (seen.get(base) ?? 0) + 1;
    seen.set(base, k);
    out.push({
      id: k === 1 ? base : `${base}#${k}`,
      seqid: r.seqid,
      type: r.type,
      role: 'feature',
      start: r.start,
      end: r.end,
      strand: r.strand,
      pieces: [r.start, r.end],
      attributes: r.attributes,
      line: r.line,
      row: r,
    });
  }
  return out;
}

/** The items a `features` mode serves, from the assembly (or the rows). */
function itemsOf(features: GffFeatures, a: Assembly, rows: readonly GffRow[]): GffItem[] {
  if (features === 'rows') return rowItems(rows);
  if (features === 'transcripts') return [...a.transcripts, ...a.features];
  const genes = new Set(a.genes.map((g) => g.id));
  return [...a.genes, ...a.transcripts.filter((t) => t.gene === undefined || !genes.has(t.gene)), ...a.features];
}

/** The fields of the items: the fixed ones, then the attributes in order of first use, at most `MaxAttributes`. */
function itemFields(features: GffFeatures, items: readonly GffItem[], causes: GffCauses): ItemField[] {
  const nan = (v: number | undefined) => (v === undefined ? Number.NaN : v);
  const fixed: ItemField[] = [
    { name: 'label', type: 'category', value: labelOf },
    { name: 'type', type: 'category', value: (it) => it.type },
  ];
  if (features === 'rows') {
    const row = (it: GffItem) => (it as GffItem & { row: GffRow }).row;
    fixed.push(
      { name: 'source', type: 'category', value: (it) => row(it).source },
      { name: 'score', type: 'float32', value: (it) => row(it).score },
      { name: 'phase', type: 'int8', value: (it) => row(it).phase },
    );
  } else
    fixed.push(
      { name: 'gene', type: 'category', value: (it) => it.gene ?? null },
      { name: 'gene_name', type: 'category', value: (it) => it.geneName ?? null },
      { name: 'thick_start', type: 'float64', value: (it) => nan(it.thickStart) },
      { name: 'thick_end', type: 'float64', value: (it) => nan(it.thickEnd) },
      { name: 'cds_phase', type: 'int8', value: (it) => it.cdsPhase ?? -1 },
    );
  const taken = new Set([...fixed.map((f) => f.name), 'strand']);
  const keys: string[] = [];
  const seen = new Set<string>();
  let past = 0;
  for (const it of items)
    for (const k in it.attributes) {
      if (seen.has(k) || NotFields.has(k) || taken.has(k)) continue;
      seen.add(k);
      if (keys.length < MaxAttributes) keys.push(k);
      else past++;
    }
  if (past) causes.add('seq.source.gff.attributes', `More than ${MaxAttributes} attribute names: ${past} not fields`);
  return [
    ...fixed,
    ...keys.map((k) => ({ name: k, type: 'category' as const, value: (it: GffItem) => it.attributes[k] ?? null })),
  ];
}

async function open(source: SourceEntry, ctx: OpenContext, options: GffOptions): Promise<OpenedSource> {
  const at = source.address;
  const format = source.format as GffFormat;
  if (!(GffFormats as readonly string[]).includes(format))
    refuse(at, 'seq.source.no_data', `Not GFF3 or GTF: ${format}`);
  const features: GffFeatures = ((f) => (f === 'genes' || f === 'rows' ? f : 'transcripts'))(
    (source as { readonly features?: string }).features,
  );
  const causes = new GffCauses();
  if (source.base !== undefined || source.ends !== undefined)
    causes.add('seq.source.gff.transform', `${format} fixes base 1 and inclusive ends`);
  if (source.index !== undefined) causes.add('seq.source.gff.index', `${format} is read whole`);

  let bytes: Uint8Array;
  if (typeof source.value === 'string') bytes = new TextEncoder().encode(source.value);
  else if (source.uri !== undefined) {
    try {
      bytes = await ctx.assets.readBytes(source.uri, ctx.signal);
    } catch {
      if (ctx.signal.aborted) throw abortError();
      refuse(at, 'seq.source.unreadable', `The file ${source.uri} could not be read`);
    }
  } else refuse(at, 'seq.source.no_data', `A ${format} source needs a uri or a text value`);
  const text = await textBytes(bytes, ctx.limits.maxBytes, ctx.signal);
  if (!text.ok)
    refuse(at, text.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed', text.message, IndexHint);
  if (ctx.signal.aborted) throw abortError();
  const isTranscript = transcriptTest(options.transcriptTypes);
  const file = parseGff(decodeText(text.bytes).text, format, isTranscript, causes, ctx.limits.maxRows);
  if (typeof file === 'string') refuse(at, 'seq.source.too_large', file, IndexHint);
  if (file.fastaBytes) causes.add('seq.source.gff.fasta', `${file.fastaBytes} characters after ##FASTA ignored`);

  // a whole-sequence feature describes its sequence (its length when no directive gives one), never an item
  const lengths = file.lengths;
  const rows = file.rows.filter((r) => {
    if (!Describing.has(r.type) || r.start !== 0) return true;
    const known = lengths.get(r.seqid);
    if (known !== undefined && known !== r.end) return true;
    if (known === undefined) lengths.set(r.seqid, r.end);
    causes.add(
      'seq.source.gff.region',
      `A ${r.type} feature over ${r.seqid} read as its description`,
      `line ${r.line}`,
    );
    return false;
  });

  const assembly = assembleTranscripts(rows, { transcriptTypes: options.transcriptTypes });
  for (const [type, count] of assembly.promoted)
    for (let k = 0; k < count; k++) causes.add('seq.source.gff.type', `\`${type}\` read as a transcript`);
  for (const l of assembly.cycles) causes.add('seq.source.gff.parent', 'A cycle of `Parent`', `line ${l}`);
  for (const l of assembly.deep) causes.add('seq.source.gff.parent', `A \`Parent\` chain deeper than 16`, `line ${l}`);
  if (format === 'gff3') {
    for (const l of assembly.orphans) causes.add('seq.source.gff.parent', 'A `Parent` no line has', `line ${l}`);
  } else if (assembly.synthetic && source.index !== undefined)
    causes.add('seq.source.gtf.extent', 'No transcript lines: transcripts made from their rows');
  for (const l of assembly.outside) causes.add('seq.source.gff.extent', 'A child outside its parent', `line ${l}`);

  const items = itemsOf(features, assembly, rows);
  const fields = itemFields(features, items, causes);
  const bySystem = new Map<string, GffItem[]>();
  for (const name of lengths.keys()) bySystem.set(name, []);
  for (const it of items) {
    let l = bySystem.get(it.seqid);
    if (!l) bySystem.set(it.seqid, (l = []));
    l.push(it);
  }
  const systems = new Map<string, SystemItems>();
  let extent = 1;
  for (const [name, list] of bySystem) {
    const s = systemItems(list);
    systems.set(name, s);
    extent = Math.max(extent, lengths.get(name) ?? 0, s.maxEnd);
  }
  const whole = wholeLevel(extent);
  const infos: SystemInfo[] = [...systems.keys()].map((name) =>
    lengths.has(name) ? { name, length: lengths.get(name)! } : { name },
  );
  const served = new Set(items);
  const names = new NameIndex(
    features === 'rows' ? items : [...assembly.genes, ...assembly.transcripts, ...assembly.features],
    NameFields[format],
    (it) => served.has(it),
  );
  const coverage: Coverage = {
    systems: infos,
    shapes: infos.map((s) => ({ system: s.name, shape: 'ranges' as const, fields: fieldInfos(fields) })),
    levels: whole > 1 ? [1, whole] : [1],
    names: NameFields[format],
  };
  const maxRows = ctx.limits.maxRows;
  return {
    coverage,
    issues: causes.issues(at),
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      const sys = systems.get(request.system);
      if (!sys) return Promise.resolve(failure('not-served', `No sequence ${request.system}`));
      if (!coverage.levels.includes(request.level))
        return Promise.resolve(failure('not-served', `No level ${request.level}`));
      if (!request.shapes.includes('ranges'))
        return Promise.resolve(failure('not-served', `The ${format} source serves ranges only`));
      const got = rowsIn(sys, request.start, request.end);
      if (got.length / 3 > maxRows)
        return Promise.resolve(failure('too-large', `${got.length / 3} rows in one region`));
      return Promise.resolve({ ok: true as const, tables: { ranges: rangesOf(sys, got, fields) } });
    },
    findNames(text, opts, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      return Promise.resolve(names.find(text, Math.max(0, Math.floor(opts.max ?? 20)), opts.prefix ?? false));
    },
    dispose() {},
  };
}

/** A GFF3 and GTF runtime with its own options (`DefaultGffOptions` for those not given). */
export function gffRuntime(over: Partial<GffOptions> = {}): FormatRuntime {
  const options = { ...DefaultGffOptions, ...over };
  return { name: 'gff', formats: GffFormats, open: (source, ctx) => open(source, ctx, options) };
}

/**
 * The runtime of `gff3` (in place of the records preset when a host composes it with `withAdapter`) and `gtf`, with the
 * default transcript types: opt in by listing it.
 */
export const GffRuntime: FormatRuntime = gffRuntime();
