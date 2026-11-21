/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import type { AssetRangeReader } from '@molstar/core/util/assets';
import { type BgzfReader, bgzfReader } from '@molstar/io/common/bgzf';
import {
  DefaultTabixLimits,
  readHeaderLines,
  readRegion,
  readTabixIndex,
  regionBytes,
  type TabixIndex,
  type TabixLimits,
  type TabixLines,
} from '@molstar/io/common/tabix';
import { SourceIssues } from '../../data/issues.js';
import { type AssetReader, abortError, type FailureKind, isAbortError, type OpenLimits } from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';

// A tab-line source read whole or by region through tabix (_spec/README.md "Data", Genome runtimes' internals): lines
// as bytes either way, so one line parser serves both and "a region read equals the whole read clipped" is a property
// of the lines alone. Shared by the BED, GFF and VCF runtimes, not published (`./runtimes/internal/*`).

export type { TabixLines } from '@molstar/io/common/tabix';

/** The lines of a file read whole: its header lines (`meta`, the first `skip`) and the rest, with 1-based numbers. */
export interface WholeLines {
  readonly header: LineSpans;
  readonly lines: LineSpans;
}

export interface LineSpans {
  readonly bytes: Uint8Array;
  readonly lineStart: Uint32Array;
  readonly lineEnd: Uint32Array;
  /** Each line's number in the file, 1-based (an issue's place). */
  readonly lineNumber: Uint32Array;
}

export interface LineOptions {
  /** A header line's first character (`#`); empty: none. */
  readonly meta: string;
  /** Lines skipped at the top as header. */
  readonly skip?: number;
  /** One line's bound (`DefaultTabixLimits.maxLineBytes`, 16 MiB). */
  readonly maxLineBytes?: number;
  /** Lines a file may hold (`limits.maxRows`). */
  readonly maxLines?: number;
}

export type LinesResult<T> =
  | { readonly ok: true; readonly value: T }
  | { readonly ok: false; readonly kind: 'too-large' | 'malformed'; readonly message: string };

/**
 * A text read whole split into lines (`\n`, a `\r` before it dropped; blank lines left out): the header lines at the top
 * (those starting with `meta`, and the first `skip`), then the rest, a `meta` line among them kept as a comment's
 * place is the format's. A line over `maxLineBytes` is `malformed`, more than `maxLines` lines `too-large`; the spans
 * are counted before they are allocated.
 */
export function wholeLines(bytes: Uint8Array, options: LineOptions): LinesResult<WholeLines> {
  const maxLine = options.maxLineBytes ?? DefaultTabixLimits.maxLineBytes;
  const maxLines = options.maxLines ?? Number.POSITIVE_INFINITY;
  const skip = options.skip ?? 0;
  const meta = options.meta ? options.meta.charCodeAt(0) : -1;
  /** Calls `line` for each line kept (the blank ones left out but a skipped one), the header ones first. */
  const scan = (line: (s: number, e: number, number: number, header: boolean) => boolean | undefined) => {
    let inHeader = true;
    for (let s = 0, number = 1; s < bytes.length; number++) {
      let e = bytes.indexOf(10, s);
      if (e < 0) e = bytes.length;
      const next = e + 1;
      if (e > s && bytes[e - 1] === 13) e--;
      if (inHeader && !(number <= skip || (e > s && bytes[s] === meta))) inHeader = false;
      if ((e > s || (inHeader && number <= skip)) && line(s, e, number, inHeader) === false) return;
      s = next;
    }
  };
  // counted first, so the spans are allocated once and a bound is refused before they are
  let headerCount = 0;
  let count = 0;
  let failure: LinesResult<WholeLines> | undefined;
  scan((s, e, number, header) => {
    if (e - s > maxLine) {
      failure = { ok: false, kind: 'malformed', message: `Line ${number} is over ${maxLine} bytes` };
      return false;
    }
    if (header) headerCount++;
    else if (count - headerCount >= maxLines) {
      failure = { ok: false, kind: 'too-large', message: `The file has more than ${maxLines} lines` };
      return false;
    }
    count++;
  });
  if (failure) return failure;
  const starts = new Uint32Array(count);
  const ends = new Uint32Array(count);
  const numbers = new Uint32Array(count);
  let n = 0;
  scan((s, e, number) => {
    starts[n] = s;
    ends[n] = e;
    numbers[n] = number;
    n++;
  });
  const spans = (from: number, to: number): LineSpans => ({
    bytes,
    lineStart: starts.subarray(from, to),
    lineEnd: ends.subarray(from, to),
    lineNumber: numbers.subarray(from, to),
  });
  return { ok: true, value: { header: spans(0, headerCount), lines: spans(headerCount, n) } };
}

/** A region's lines, or why there are none: a `FailureKind` the runtime answers with, and a message. */
export type RegionLines =
  | { readonly ok: true; readonly lines: TabixLines }
  | { readonly ok: false; readonly kind: FailureKind; readonly message: string };

/** A tab-line file by region through its tabix index (TBI or CSI). */
export interface IndexedLines {
  readonly index: TabixIndex;
  /** Whether the file ends with BGZF's EOF block (false: it may be cut short; the runtime says it, the file is read). */
  readonly eof: boolean;
  /** The header lines at the top of the file. */
  header(signal: AbortSignal): Promise<RegionLines>;
  /** Lines meeting `[start, end)` (0-based, half-open) of the sequence `name`, each with its extent by the index. */
  region(name: string, start: number, end: number, signal: AbortSignal): Promise<RegionLines>;
  /** A region's compressed bytes from the index, nothing read (`regionBytes`): a tile's read cost. */
  cost(name: string, start: number, end: number): number;
  dispose(): void;
}

/** A range reader over bytes in memory (a host without range reads: the file read whole, bounded). */
function bytesRanges(bytes: Uint8Array): AssetRangeReader {
  return {
    size: bytes.length,
    read: (ranges, signal) => {
      if (signal?.aborted) return Promise.reject(abortError());
      return Promise.resolve(
        ranges.map((r) => bytes.subarray(r.offset, Math.min(bytes.length, r.offset + Math.max(0, r.length)))),
      );
    },
    dispose() {},
  };
}

/** A failure of io's readers or a range read as an answer's kind: io's `kind`, a range read `failed`. */
function failureOf(e: unknown): { kind: FailureKind; message: string } {
  const name = (e as Error)?.name;
  const message = (e as Error)?.message ?? `${e}`;
  if (name === 'TabixError' || name === 'BgzfError') {
    return { kind: (e as { kind: 'malformed' | 'too-large' }).kind, message };
  }
  return { kind: 'failed', message };
}

async function guarded<T>(
  signal: AbortSignal,
  f: () => Promise<T>,
): Promise<T | { failure: { kind: FailureKind; message: string } }> {
  try {
    return await f();
  } catch (e) {
    if (isAbortError(e) || signal.aborted) throw abortError();
    return { failure: failureOf(e) };
  }
}

function openError(
  address: string,
  code: 'seq.source.unreadable' | 'seq.source.malformed' | 'seq.source.too_large',
  message: string,
) {
  return new SourceOpenError([issueOf(SourceIssues, code, message, { address })]);
}

export interface IndexedLinesContext {
  readonly assets: AssetReader;
  readonly signal: AbortSignal;
  readonly limits: OpenLimits;
  /** The source's address (an issue's place). */
  readonly address: string;
  /** io's tabix bounds (defaults: `DefaultTabixLimits`). */
  readonly tabix?: Partial<TabixLimits>;
}

/**
 * Opens a BGZF file by range (`assets.openRange`; a host without it reads the file whole within `limits.maxBytes`) and
 * its index (`readBytes`, within `maxIndexBytes`), and checks the EOF block once. A file or an index that cannot be read
 * is a `SourceOpenError`: `seq.source.unreadable`, a malformed index `seq.source.malformed`, one past its bound
 * `seq.source.too_large`.
 */
export async function openIndexedLines(uri: string, indexUri: string, ctx: IndexedLinesContext): Promise<IndexedLines> {
  const { assets, signal, address } = ctx;
  const tabix = { ...DefaultTabixLimits, ...ctx.tabix };
  let ranges: AssetRangeReader;
  try {
    if (assets.openRange) ranges = await assets.openRange(uri, signal);
    else {
      const bytes = await assets.readBytes(uri, signal);
      if (bytes.length > ctx.limits.maxBytes) {
        throw openError(
          address,
          'seq.source.too_large',
          `The file ${uri} is over ${ctx.limits.maxBytes} bytes; a host with range reads reads it by region`,
        );
      }
      ranges = bytesRanges(bytes);
    }
  } catch (e) {
    if (isAbortError(e) || signal.aborted) throw abortError();
    if (e instanceof SourceOpenError) throw e;
    throw openError(address, 'seq.source.unreadable', `The file ${uri} could not be read`);
  }
  let index: TabixIndex;
  let reader: BgzfReader | undefined;
  try {
    let indexBytes: Uint8Array;
    try {
      indexBytes = await assets.readBytes(indexUri, signal);
    } catch (e) {
      if (isAbortError(e) || signal.aborted) throw abortError();
      throw openError(address, 'seq.source.unreadable', `The index ${indexUri} could not be read`);
    }
    try {
      index = await readTabixIndex(indexBytes, tabix, signal);
    } catch (e) {
      if (isAbortError(e) || signal.aborted) throw abortError();
      const f = failureOf(e);
      throw openError(
        address,
        f.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed',
        `The index ${indexUri} cannot be read: ${f.message}`,
      );
    }
    reader = bgzfReader(ranges, tabix);
    let eof: boolean;
    try {
      eof = await reader.hasEof(signal);
    } catch (e) {
      if (isAbortError(e) || signal.aborted) throw abortError();
      throw openError(address, 'seq.source.unreadable', `The file ${uri} could not be read`);
    }
    const r = reader;
    const answer = (got: TabixLines | { failure: { kind: FailureKind; message: string } }): RegionLines =>
      'failure' in got ? { ok: false, ...got.failure } : { ok: true, lines: got };
    return {
      index,
      eof,
      header: async (s) => answer(await guarded(s, () => readHeaderLines(r, index, tabix, s))),
      region: async (name, start, end, s) =>
        answer(await guarded(s, () => readRegion(r, index, name, start, end, tabix, s))),
      cost: (name, start, end) => regionBytes(index, name, start, end),
      dispose() {
        r.dispose();
        ranges.dispose();
      },
    };
  } catch (e) {
    reader?.dispose();
    ranges.dispose();
    throw e;
  }
}
