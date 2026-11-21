/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { SourceIssues } from '../data/issues.js';
import {
  abortError,
  type Coverage,
  type FailureKind,
  type FormatRuntime,
  type OpenContext,
  type OpenedSource,
  type PairsList,
  type RegionAnswer,
  type RegionRequest,
  type SourceEntry,
  type SystemInfo,
} from '../data/protocol.js';
import { wholeLevel } from '../data/seq-format.js';
import { BedCauses } from './bed-issues.js';
import { indexPairs, PairsBuilder, pairsIn, pairsList, pairsSchema, rowBytes } from './bed-pairs.js';
import { type BedLine, type DeclaredType, fieldInfos, MaxPosition } from './bed-parse.js';
import {
  answerOf,
  type BedFormat,
  BedFormats,
  type BedOptions,
  bedgraphIntervals,
  bedItems,
  dataCells,
  failure,
  IndexHint,
  type Opened,
  refuse,
  refuseField,
  schemaOf,
  spanLines,
} from './bed-read.js';
import { wholeLines } from './internal/whole-lines.js';
import { textBytes } from './internal/text.js';

// The BED family (_spec/README.md "Data", BED family), the opt-in runtime `@molstar/seqx-core/runtimes/bed`: BED (BED3
// to BED12, `bed N+M` by `fields`), bedGraph, BEDPE and chrom.sizes. BED and bedGraph read whole, or by region at level
// 1 through a tabix `index` (`bed-indexed.ts`, loaded then), by one line parser either way; BEDPE read whole and indexed
// by both ends; chrom.sizes a coverage of lengths. Positions as the formats write them (base 0, ends exclusive): the
// source's `base` and `ends` are not read (EC-195).

export { BedIssues } from './bed-issues.js';
export { type BedFormat, BedFormats, type BedOptions } from './bed-read.js';

export const DefaultBedOptions: BedOptions = { maxPairsBytes: 16 * 1024 * 1024 };

/** A pairs read's answer: the pairs list, or why there is none (the protocol's failure kinds). */
export type PairsAnswer =
  | { readonly ok: true; readonly pairs: PairsList }
  | { readonly ok: false; readonly kind: FailureKind; readonly issue: Issue };

/**
 * An opened BEDPE: its systems (no shape yet: the protocol's `pairs` shape is to come), every pair with an end in a
 * region by `pairs`, and what it holds (`summaryBytes`, counted by the store as a pyramid's).
 */
export interface BedpeSource extends OpenedSource {
  readonly summaryBytes: number;
  /** Rows with an end meeting `[start, end)` of `system`, each once, in file order (EC-208). */
  pairs(system: string, start: number, end: number, signal: AbortSignal): Promise<PairsAnswer>;
}

/** The text a whole read parses: inline `value` or the file at `uri`, inflated by its magic within `maxBytes`. */
async function wholeText(o: Opened): Promise<BedLine[]> {
  const { source, ctx, causes } = o;
  const at = source.address;
  let bytes: Uint8Array;
  if (typeof source.value === 'string') bytes = new TextEncoder().encode(source.value);
  else if (source.uri !== undefined) {
    try {
      bytes = await ctx.assets.readBytes(source.uri, ctx.signal);
    } catch {
      if (ctx.signal.aborted) throw abortError();
      refuse(at, 'seq.source.unreadable', `The file ${source.uri} could not be read`);
    }
  } else refuse(at, 'seq.source.no_data', `A ${o.format} source needs a uri or a text value`);
  const text = await textBytes(bytes, ctx.limits.maxBytes, ctx.signal);
  if (!text.ok)
    refuse(at, text.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed', text.message, IndexHint);
  const lines = wholeLines(text.bytes, {
    meta: '',
    maxLines: ctx.limits.maxRows,
    ...(o.options.tabix?.maxLineBytes !== undefined && { maxLineBytes: o.options.tabix.maxLineBytes }),
  });
  if (!lines.ok)
    refuse(at, lines.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed', lines.message, IndexHint);
  const l = lines.value.lines;
  return spanLines(l.bytes, l.lineStart, l.lineEnd, (i) => `line ${l.lineNumber[i]}`, causes);
}

/** A source over what a whole read holds: a request answered from it, by system and level. */
function wholeSource(coverage: Coverage, issues: Issue[], answer: (r: RegionRequest) => RegionAnswer): OpenedSource {
  return {
    coverage,
    issues,
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      return Promise.resolve(
        !coverage.systems.some((s) => s.name === request.system)
          ? failure('not-served', `No sequence ${request.system}`)
          : !coverage.levels.includes(request.level)
            ? failure('not-served', `No level ${request.level}`)
            : answer(request),
      );
    },
    dispose() {},
  };
}

/** The data lines' count and the first one's columns. */
function firstLine(text: BedLine[]) {
  let columns = 0;
  let count = 0;
  for (const r of dataCells(text, new BedCauses())) {
    if (!count) columns = r.cells.length;
    count++;
  }
  return { columns, count };
}

async function openWhole(o: Opened): Promise<OpenedSource> {
  const { source, ctx, causes, format } = o;
  const text = await wholeText(o);
  const maxRows = ctx.limits.maxRows;
  const issues = () => causes.issues(source.address);
  const rows = () => dataCells(text, causes);

  if (format === 'bed') {
    const schema = schemaOf(source, Math.max(3, firstLine(text).columns));
    const systems = bedItems(rows(), schema, causes, maxRows);
    if (typeof systems === 'string') refuse(source.address, 'seq.source.too_large', systems, IndexHint);
    let extent = 1;
    for (const items of systems.values()) extent = Math.max(extent, items.maxEnd);
    // read whole, items are served at the level whose one tile holds the longest system too (the inline rule)
    const whole = wholeLevel(extent);
    const names = [...systems.keys()];
    const fields = fieldInfos(schema);
    const coverage: Coverage = {
      systems: names.map((name) => ({ name })),
      shapes: names.map((system) => ({ system, shape: 'ranges' as const, fields })),
      levels: whole > 1 ? [1, whole] : [1],
    };
    return wholeSource(coverage, issues(), (r) => answerOf('bed', systems.get(r.system), r, maxRows));
  }
  if (format === 'bedgraph') {
    const systems = bedgraphIntervals(rows(), causes);
    const names = [...systems.keys()];
    const coverage: Coverage = {
      systems: names.map((name) => ({ name })),
      shapes: names.map((system) => ({ system, shape: 'values' as const })),
      levels: [1],
    };
    return wholeSource(coverage, issues(), (r) => answerOf('bedgraph', systems.get(r.system), r, maxRows));
  }
  if (format === 'chrom-sizes') {
    const systems: SystemInfo[] = [];
    const refused: string[] = [];
    for (const { cells, place } of rows()) {
      const [name, length] = cells;
      if (length !== undefined && /^\d+$/.test(length) && Number(length) > MaxPosition) refused.push(name);
      else if (cells.length < 2 || !/^\d+$/.test(length))
        causes.add('seq.source.bed.line', 'Not a name and a whole length', place);
      else if (systems.some((s) => s.name === name))
        causes.add('seq.source.bed.line', 'A name again (the first kept)', place);
      else systems.push({ name, length: Number(length) });
    }
    const all = issues();
    if (refused.length)
      all.push(
        issueOf(
          SourceIssues,
          'seq.source.too_large',
          `Longer than 2^31 − 1, not read: ${refused.slice(0, 10).join(', ')}`,
          {
            address: source.address,
          },
        ),
      );
    return wholeSource({ systems, shapes: [], levels: [1] }, all, () => failure('not-served', 'Lengths only'));
  }
  return openPairs(o, text);
}

function openPairs(o: Opened, text: BedLine[]): BedpeSource {
  const { source, ctx, causes, options } = o;
  const { columns, count } = firstLine(text);
  const schema = pairsSchema(
    columns || 6,
    source.fields,
    source.types as Readonly<Record<string, DeclaredType>> | undefined,
  );
  if ('error' in schema) refuseField(source.address, schema.error);
  const check = (bytes: number) => {
    if (bytes > options.maxPairsBytes)
      refuse(
        source.address,
        'seq.source.too_large',
        `The BEDPE source holds ${bytes} bytes with its two orders, over ${options.maxPairsBytes}`,
        'A format indexed by both ends (pairix), or a levelled pairs format',
      );
  };
  // the floor of what the lines will hold, refused before any is parsed
  check(count * rowBytes(schema, 1));
  const builder = new PairsBuilder(schema);
  for (const { cells, place } of dataCells(text, causes)) builder.add(cells, place, causes);
  const bytes = builder.rows * rowBytes(schema, builder.names.length);
  check(bytes);
  const rows = builder.finish();
  const index = indexPairs(rows);
  const coverage: Coverage = { systems: rows.names.map((name) => ({ name })), shapes: [], levels: [1] };
  return {
    ...wholeSource(coverage, causes.issues(source.address), () => failure('not-served', 'Pairs only (`pairs`)')),
    summaryBytes: bytes,
    pairs(system, start, end, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      const s = rows.names.indexOf(system);
      if (s < 0) return Promise.resolve(failure('not-served', `No sequence ${system}`));
      const found = pairsIn(rows, index, s, start, end);
      return Promise.resolve(
        found.length > ctx.limits.maxRows
          ? failure('too-large', `${found.length} pairs in one region`)
          : { ok: true, pairs: pairsList(rows, found) },
      );
    },
  };
}

async function open(source: SourceEntry, ctx: OpenContext, options: BedOptions): Promise<OpenedSource> {
  const format = source.format as BedFormat;
  if (!BedFormats.includes(format)) refuse(source.address, 'seq.source.no_data', `Not a BED format: ${format}`);
  const causes = new BedCauses();
  if (source.base !== undefined || source.ends !== undefined)
    causes.add('seq.source.bed.transform', `${format} fixes base 0 and exclusive ends`);
  const o: Opened = { source, ctx, options, causes, format };
  if (source.index !== undefined) {
    if (format === 'bedpe' || format === 'chrom-sizes') causes.add('seq.source.bed.index', `${format} is read whole`);
    else if (source.uri !== undefined)
      return (await import('./bed-indexed.js')).openIndexed(o, source.uri, source.index);
  }
  return openWhole(o);
}

/** A runtime of the BED family with its own bounds (`DefaultBedOptions` for those not given). */
export function bedRuntime(over: Partial<BedOptions> = {}): FormatRuntime {
  const options = { ...DefaultBedOptions, ...over };
  return { name: 'bed', formats: BedFormats, open: (source, ctx) => open(source, ctx, options) };
}

/** The runtime of `bed`, `bedgraph`, `bedpe` and `chrom-sizes` with the default bounds: opt in by listing it. */
export const BedRuntime: FormatRuntime = bedRuntime();
