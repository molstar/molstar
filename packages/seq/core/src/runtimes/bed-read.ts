/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import type { TabixLimits } from '@molstar/io/common/tabix';
import type { Issue } from '@molstar/core/repr/tree';
import {
  type FailureKind,
  type OpenContext,
  type OpenedSource,
  type PairsList,
  type RangesList,
  type RegionAnswer,
  type RegionRequest,
  type SourceEntry,
} from '../data/protocol.js';
import { SourceOpenError } from '../data/seq-runtime.js';
import { BedCauses, BedIssues } from './bed-issues.js';
import {
  type BedLine,
  type BedSchema,
  bedSchema,
  cellsOf,
  type DeclaredType,
  type Intervals,
  IntervalsBuilder,
  type Items,
  ItemsBuilder,
  isMetaLine,
  binsIn,
  piecesIn,
  rangesList,
  valuesIn,
} from './bed-parse.js';
import { failure } from './internal/answers.js';
import { densityOf } from './internal/levels.js';

const EmptyIntervals: Intervals = {
  start: new Int32Array(0),
  end: new Int32Array(0),
  value: new Float32Array(0),
  maxSpan: 1,
};
import { decodeText } from './internal/text.js';

// What the BED family's whole and region reads share (_spec/README.md "Data", BED family): data lines into items or
// intervals, and the answer of a request from them, so a region read equals the whole read clipped.

export const BedFormats = ['bed', 'bedgraph', 'bedpe', 'chrom-sizes'] as const;
export type BedFormat = (typeof BedFormats)[number];

/** What the runtime may hold and read, whatever a file says (security.md: a file is untrusted input). */
/**
 * Region-reads' `summaries.maxBytes`: what a source holds besides the store's entries (a pairs pyramid, a BEDPE's rows
 * and orders), 16 MiB. The store's setting when it has one (RG); until then the BED runtime's `maxPairsBytes` default.
 */
export const SummariesMaxBytes = 16 * 1024 * 1024;

export interface BedOptions {
  /** A BEDPE's rows with both orders, in bytes: region-reads' `summaries.maxBytes` (16 MiB, 600,000 six-column rows). */
  readonly maxPairsBytes: number;
  /** The coarsest level listed by region (4,096: 4 Mb tiles); a level is listed when each of its tiles fits. */
  readonly maxLevel: number;
  /** Bins a tile (the store's `tileBins`, 1,024): a level's tiles are `level × tileBins` positions. */
  readonly tileBins: number;
  /** io's tabix bounds for a source read by region (`DefaultTabixLimits` for those not given). */
  readonly tabix?: Partial<TabixLimits>;
}

/** An open in progress: the source, its context, the runtime's options and the causes said. */
export interface Opened {
  readonly source: SourceEntry;
  readonly ctx: OpenContext;
  readonly options: BedOptions;
  readonly causes: BedCauses;
  readonly format: BedFormat;
}

/** Refuses a source whose `fields` name no column (`seq.source.bed.field`). */
export function refuseField(address: string, message: string): never {
  throw new SourceOpenError([
    issueOf(BedIssues, 'seq.source.bed.field', message, {
      address,
      hint: BedIssues.codes['seq.source.bed.field'].hint,
    }),
  ]);
}

/** The schema of a BED source's first line, or the source refused. */
export function schemaOf(source: SourceEntry, columns: number): BedSchema {
  const s = bedSchema(columns, source.fields, source.types as Readonly<Record<string, DeclaredType>> | undefined);
  return 'error' in s ? refuseField(source.address, s.error) : s;
}

/** Data lines (meta and blank lines read past), each split into cells; whitespace-split lines said once. */
export function* dataCells(lines: Iterable<BedLine>, causes: BedCauses): Generator<{ cells: string[]; place: string }> {
  let spaces = false;
  for (const l of lines) {
    if (!l.text || isMetaLine(l.text)) continue;
    const c = cellsOf(l.text);
    if (c.spaces && !spaces) {
      spaces = true;
      causes.add('seq.source.bed.whitespace', 'Lines without a tab, split on spaces', l.place);
    }
    if (c.cells.length) yield { cells: c.cells, place: l.place };
  }
}

/** A text's lines from byte spans, U+FFFD said once. */
export function spanLines(
  bytes: Uint8Array,
  starts: ArrayLike<number>,
  ends: ArrayLike<number>,
  place: (i: number) => string,
  causes: BedCauses,
): BedLine[] {
  const out: BedLine[] = [];
  let replaced = false;
  for (let i = 0; i < starts.length; i++) {
    const d = decodeText(bytes, starts[i], ends[i]);
    if (d.replaced && !replaced) {
      replaced = true;
      causes.add('seq.source.bed.encoding', 'Replacement characters (U+FFFD)', place(i));
    }
    out.push({ text: d.text, place: place(i) });
  }
  return out;
}

type Rows = Iterable<{ cells: string[]; place: string }>;

/** A BED's items by system; past `maxRows` pieces (BED12's blocks counted each) the message why. */
export function bedItems(
  rows: Rows,
  schema: BedSchema,
  causes: BedCauses,
  maxRows: number,
): Map<string, Items> | string {
  const systems = new Map<string, ItemsBuilder>();
  let pieces = 0;
  for (const { cells, place } of rows) {
    if (cells.length < 3) {
      causes.add('seq.source.bed.line', 'Fewer than 3 columns', place);
      continue;
    }
    let b = systems.get(cells[0]);
    if (!b) systems.set(cells[0], (b = new ItemsBuilder(schema)));
    const before = b.pieces;
    // the rows left bound a line's blocks before they are read (a line of a million `0,` blocks is refused unsplit)
    const got = b.add(cells, place, causes, maxRows - pieces);
    pieces += b.pieces - before;
    if (got === 'too-large' || pieces > maxRows) return `More than ${maxRows} rows (BED12’s blocks counted each)`;
  }
  const out = new Map<string, Items>();
  for (const [name, b] of systems) out.set(name, b.finish());
  return out;
}

export function bedgraphIntervals(rows: Rows, causes: BedCauses): Map<string, Intervals> {
  const systems = new Map<string, IntervalsBuilder>();
  for (const { cells, place } of rows) {
    let b = systems.get(cells[0]);
    if (!b) systems.set(cells[0], (b = new IntervalsBuilder()));
    b.add(cells, place, causes);
  }
  const out = new Map<string, Intervals>();
  for (const [name, b] of systems) out.set(name, b.finish());
  return out;
}

/**
 * The answer of a BED system's items or a bedGraph's intervals over a request (whole and region reads alike). With
 * `coarse` (by region), a level above 1 not asked as items is a BED's density or a bedGraph's bins, counted from the
 * region's own records.
 */
export function answerOf(
  format: 'bed' | 'bedgraph',
  data: Items | Intervals | undefined,
  request: RegionRequest,
  maxRows: number,
  coarse = false,
): RegionAnswer {
  const { start, end } = request;
  const shape = format === 'bed' ? 'ranges' : 'values';
  if (!request.shapes.includes(shape)) return failure('not-served', `The ${format} source serves ${shape} only`);
  const bins = coarse && request.level > 1 && !request.items;
  if (format === 'bed') {
    const rows = data ? piecesIn(data as Items, start, end) : [];
    if (rows.length > maxRows) return failure('too-large', `${rows.length} rows in one region`);
    const list: RangesList = data
      ? rangesList(data as Items, rows)
      : { form: 'list', start: new Int32Array(0), end: new Int32Array(0) };
    return { ok: true, tables: { ranges: bins ? densityOf(list, request) : list } };
  }
  if (bins)
    return { ok: true, tables: { values: binsIn((data as Intervals) ?? EmptyIntervals, start, end, request.level) } };
  if (end - start > maxRows) return failure('too-large', `${end - start} positions in one region`);
  return {
    ok: true,
    tables: {
      values: data ? valuesIn(data as Intervals, start, end) : { form: 'positions', start, value: new Float32Array(0) },
    },
  };
}

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

/** The data lines' count and the first one's columns. */
export function dataLines(text: BedLine[]) {
  let columns = 0;
  let count = 0;
  for (const r of dataCells(text, new BedCauses())) {
    if (!count) columns = r.cells.length;
    count++;
  }
  return { columns, count };
}
