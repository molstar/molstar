/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  abortError,
  type Coverage,
  type FormatRuntime,
  type OpenContext,
  type OpenedSource,
  type SourceEntry,
  type SystemInfo,
} from '../data/protocol.js';
import type { TabixLimits } from '@molstar/io/common/tabix';
import { wholeLevel } from '../data/seq-format.js';
import { assembleTranscripts, DefaultTranscriptTypes, type GffItem, transcriptTest } from './gff-assemble.js';
import { GffCauses } from './gff-issues.js';
import { type GffFormat, parseGff } from './gff-parse.js';
import {
  describes,
  fieldInfos,
  type GffFeatures,
  itemFields,
  itemsOf,
  NameIndex,
  rangesOf,
  rowsIn,
  type SystemItems,
  systemItems,
} from './gff-read.js';
import { failure, refuse } from './internal/answers.js';
import { decodeText, textBytes } from './internal/text.js';

// GFF3 and GTF (_spec/README.md "Data", GFF3 and GTF), the opt-in runtime `@molstar/seqx-core/runtimes/gff`: a file read
// whole within `limits` (plain, gzip or BGZF by its magic), its lines assembled into transcripts by their type
// (`assembleTranscripts`), served as `ranges` by the source's `features` (`transcripts`, `genes` or `rows`), each item's
// pieces with its whole extent; `##sequence-region` lengths (EC-188); a name index over GFF3's `Name`, `ID`, `Alias` or
// GTF's `gene_name`, `gene_id`, `transcript_name`, `transcript_id` (EC-228). Positions as the formats write them (base
// 1, ends inclusive), converted at the edge: the source's `base` and `ends` are not read (EC-195). With a tabix `index`
// read by region (`gff-indexed.ts`, loaded then): each tile's lines assembled alone, by type.

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

export type { GffFeatures } from './gff-read.js';

export interface GffOptions {
  /** The types read as transcripts (any type ending in `_transcript` too). */
  readonly transcriptTypes: readonly string[];
  /** The coarsest level listed by region (4,096: 4 Mb tiles); a level is listed when each of its tiles fits. */
  readonly maxLevel: number;
  /** Bins a tile (the store's `tileBins`, 1,024): a level's tiles are `level × tileBins` positions. */
  readonly tileBins: number;
  /** io's tabix bounds for a source read by region (`DefaultTabixLimits` for those not given). */
  readonly tabix?: Partial<TabixLimits>;
}

export const DefaultGffOptions: GffOptions = {
  transcriptTypes: DefaultTranscriptTypes,
  maxLevel: 4096,
  tileBins: 1024,
};

/** The fields `findNames` searches, by format. */
const NameFields: Readonly<Record<GffFormat, readonly string[]>> = {
  gff3: ['Name', 'ID', 'Alias'],
  gtf: ['gene_name', 'gene_id', 'transcript_name', 'transcript_id'],
};

const IndexHint = 'A smaller file, or a host’s adapter reading it by region';

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
  if (source.index !== undefined && source.uri !== undefined)
    return (await import('./gff-indexed.js')).openIndexed(
      { source, ctx, options, causes, format, features },
      source.uri,
      source.index,
    );
  if (source.index !== undefined) causes.add('seq.source.gff.index', `${format} with an index needs a uri: read whole`);

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
    if (!describes(r, lengths)) return true;
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
