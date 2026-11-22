/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import type { AssetRangeReader } from '@molstar/core/util/assets';
import { type BgzfReader, bgzfReader, VirtualOffset } from '@molstar/io/common/bgzf';
import {
  DefaultTabixLimits,
  regionChunks,
  readHeaderLines,
  readRegion,
  readTabixIndex,
  regionBytes,
  type TabixIndex,
  type TabixLimits,
  type TabixLines,
} from '@molstar/io/common/tabix';
import { SourceIssues } from '../../data/issues.js';
import { indexedExtent, indexSpan, LevelCandidates, tabixLevels } from './levels.js';
import { type AssetReader, abortError, type FailureKind, isAbortError, type OpenLimits } from '../../data/protocol.js';
import { SourceOpenError } from '../../data/seq-runtime.js';

// A tab-line source read whole or by region through tabix (_spec/README.md "Data", Genome runtimes' internals): lines
// as bytes either way, so one line parser serves both and "a region read equals the whole read clipped" is a property
// of the lines alone. Shared by the BED, GFF and VCF runtimes, not published (`./runtimes/internal/*`).

export type { TabixLines } from '@molstar/io/common/tabix';

export {
  type LineOptions,
  type LineSpans,
  type LinesResult,
  type WholeLines,
  wholeLines,
} from './whole-lines.js';

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
  /**
   * The file's first data line (its bytes, without the newline; `#`, `track` and `browser` lines read past), or undefined
   * when the index lists none: what a format
   * whose columns are told by its first line (BED's `N`) reads at open. Within `maxLineBytes`.
   */
  head(
    signal: AbortSignal,
  ): Promise<{ readonly ok: true; readonly line?: Uint8Array } | Exclude<RegionLines, { ok: true }>>;
  /** A region's compressed bytes from the index, nothing read (`regionBytes`): a tile's read cost. */
  cost(name: string, start: number, end: number): number;
  /** io's bound on a region's compressed bytes (`maxRegionBytes`): a level's every tile within it. */
  readonly maxRegionBytes: number;
  dispose(): void;
}

/**
 * The levels a source by region lists (`tabixLevels`): each indexed sequence to its length when the file says one, else
 * to the end of its last 4 Mb tile with a cost in the index.
 */
export function indexedLevels(
  lines: Pick<IndexedLines, 'index' | 'cost' | 'maxRegionBytes'>,
  lengths: ReadonlyMap<string, number | undefined>,
  options: { readonly maxLevel: number; readonly tileBins: number },
): number[] {
  const { index } = lines;
  const span = indexSpan(index.minShift, index.depth);
  const cost = (name: string, start: number, end: number) => lines.cost(name, start, end);
  const coarsest = LevelCandidates[LevelCandidates.length - 1] * options.tileBins;
  const systems = index.names.map((name) => ({
    name,
    extent: Math.min(span, lengths.get(name) ?? indexedExtent(cost, name, coarsest, span)),
  }));
  return tabixLevels(cost, systems, { ...options, maxRegionBytes: lines.maxRegionBytes });
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

/** Whether a line is read past before the data: the index's meta lines (`#`), and `track` and `browser` lines. */
function isHeadLine(bytes: Uint8Array, meta: number): boolean {
  if (bytes[0] === meta) return true;
  const word = (w: string) =>
    bytes.length >= w.length &&
    [...w].every((c, i) => bytes[i] === c.charCodeAt(0)) &&
    (bytes.length === w.length || bytes[w.length] === 9 || bytes[w.length] === 32);
  return word('track') || word('browser');
}

/**
 * The first data line: from the least chunk start the index lists (htslib's first line), block by block (a block's end
 * told by its `BSIZE`, read from the file's bytes), the meta, `track` and `browser` lines read past, each line within
 * `maxLine` bytes and the lines read past within `maxHead` together.
 */
async function firstLine(
  reader: BgzfReader,
  ranges: AssetRangeReader,
  index: TabixIndex,
  maxLine: number,
  maxHead: number,
  signal: AbortSignal,
): Promise<Uint8Array | undefined> {
  let first = Number.POSITIVE_INFINITY;
  // the index's whole span (a .tbi's 2^29), never past int32
  const span = indexSpan(index.minShift, index.depth);
  for (const name of index.names) {
    const c = regionChunks(index, name, 0, span)[0];
    if (c && c.begin < first) first = c.begin;
  }
  if (!Number.isFinite(first)) return undefined;
  const meta = index.meta ? index.meta.charCodeAt(0) : -1;
  const malformed = (why: string) => Object.assign(new Error(why), { name: 'TabixError', kind: 'malformed' });
  /** The line being read, its parts across blocks. */
  let parts: Uint8Array[] = [];
  let length = 0;
  let passed = 0;
  const line = () => {
    const out = new Uint8Array(length);
    for (let i = 0, o = 0; i < parts.length; o += parts[i].length, i++) out.set(parts[i], o);
    return out.length && out[out.length - 1] === 13 ? out.subarray(0, out.length - 1) : out;
  };
  for (let at = first; ; ) {
    const block = VirtualOffset.block(at);
    if (block + 18 > ranges.size) break;
    const [h] = await ranges.read([{ offset: block, length: 18 }], signal);
    if (h.length < 18) break;
    const next = block + (h[16] | (h[17] << 8)) + 1;
    const bytes = await reader.read(at, VirtualOffset.of(next, 0), signal);
    if (bytes.length === 0) break;
    for (let s = 0; s < bytes.length; ) {
      const nl = bytes.indexOf(10, s);
      const e = nl < 0 ? bytes.length : nl;
      parts.push(bytes.subarray(s, e));
      length += e - s;
      if (length > maxLine) throw malformed(`A line is over ${maxLine} bytes`);
      if (nl < 0) break;
      const l = line();
      if (l.length && !isHeadLine(l, meta)) return l;
      passed += length + 1;
      if (passed > maxHead) throw malformed(`No data line in the first ${maxHead} bytes`);
      parts = [];
      length = 0;
      s = nl + 1;
    }
    if (next >= ranges.size) break;
    at = VirtualOffset.of(next, 0);
  }
  const l = line();
  return l.length && !isHeadLine(l, meta) ? l : undefined;
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
      head: async (s) => {
        const got = await guarded(s, () => firstLine(r, ranges, index, tabix.maxLineBytes, tabix.maxHeaderBytes, s));
        return got instanceof Uint8Array || got === undefined
          ? { ok: true, ...(got && { line: got }) }
          : { ok: false, ...got.failure };
      },
      cost: (name, start, end) => regionBytes(index, name, start, end),
      maxRegionBytes: tabix.maxRegionBytes,
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
