/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { BgzfReader } from './bgzf.js';
import {
  BgzfError,
  type BgzfLimits,
  BgzfRangeReader,
  DefaultBgzfLimits,
  inflateBgzf,
  isBgzf,
  VirtualOffset,
} from './internal/bgzf-blocks.js';
import { binCount, binFirst, binParent, metaBin, reg2bins } from './internal/tabix-bins.js';
import { type Extent, lineExtent, MaxPosition, sameBytes } from './internal/tabix-lines.js';

// Tabix (`_spec/bgzf.md`): TBI and CSI indexes of a BGZF text file sorted by position, and a region's lines as bytes
// with each line's extent by the index's rule (htslib's `tbx_parse1`). Every count and offset is checked against the
// limits and the index's own bytes before it allocates.

export interface TabixIndex {
  readonly kind: 'tbi' | 'csi';
  readonly format: 'generic' | 'sam' | 'vcf';
  /** The 0x10000 flag (BED): positions 0-based, ends exclusive. */
  readonly zeroBased: boolean;
  /** 1-based column numbers; `end` 0: none. */
  readonly columns: { readonly seq: number; readonly begin: number; readonly end: number };
  /** A header line's first character (`#`); empty when the index names none. */
  readonly meta: string;
  /** Lines skipped at the top. */
  readonly skip: number;
  /** The sequences, in index order. */
  readonly names: readonly string[];
  /** 14 for a .tbi. */
  readonly minShift: number;
  /** 5 for a .tbi. */
  readonly depth: number;
}

export interface TabixLimits extends BgzfLimits {
  /** An index inflated (64 MiB). */
  readonly maxIndexBytes: number;
  /** A region's chunks, compressed, together (8 MiB). */
  readonly maxRegionBytes: number;
  /** One line (16 MiB: a VCF line of 100,000 samples at 160 bytes each). */
  readonly maxLineBytes: number;
  /** The header lines (16 MiB). */
  readonly maxHeaderBytes: number;
}

export const DefaultTabixLimits: TabixLimits = {
  ...DefaultBgzfLimits,
  maxIndexBytes: 64 * 1024 * 1024,
  maxRegionBytes: 8 * 1024 * 1024,
  maxLineBytes: 16 * 1024 * 1024,
  maxHeaderBytes: 16 * 1024 * 1024,
};

/** An index or a region that cannot be read (`malformed`), or a bound passed (`too-large`); by `name` across copies. */
export class TabixError extends Error {
  constructor(
    readonly kind: 'malformed' | 'too-large',
    message: string,
  ) {
    super(message);
    this.name = 'TabixError';
  }
}

export interface Chunk {
  readonly begin: VirtualOffset;
  readonly end: VirtualOffset;
}

/** Lines meeting a region, as bytes: each line's span in `bytes` and its extent by the index's columns (tbx_parse1). */
export interface TabixLines {
  readonly bytes: Uint8Array;
  readonly lineStart: Uint32Array;
  readonly lineEnd: Uint32Array;
  /** 0-based, half-open; −1 for a header line. */
  readonly begin: Int32Array;
  readonly end: Int32Array;
}

/** A reference's bins (chunks as begin, end pairs; CSI's `loffset`) and TBI's linear index. */
interface Reference {
  readonly bins: Map<number, { readonly loffset: number; readonly chunks: Float64Array }>;
  readonly linear: Float64Array;
}

interface IndexData {
  readonly refs: readonly Reference[];
  readonly byName: ReadonlyMap<string, number>;
  /** The first line's virtual offset: where the header ends; undefined when no line is indexed. */
  readonly firstLine: number | undefined;
}

/** The parsed bins of each index, kept off the public object. */
const Data = new WeakMap<TabixIndex, IndexData>();

function dataOf(index: TabixIndex): IndexData {
  const d = Data.get(index);
  if (!d) throw new TabixError('malformed', 'not an index read by readTabixIndex');
  return d;
}

function limitsOf(limits?: Partial<TabixLimits>): TabixLimits {
  return { ...DefaultTabixLimits, ...limits };
}

/** A bounds-checked little-endian cursor over the inflated index. */
class Cursor {
  private readonly view: DataView;
  at = 0;
  constructor(readonly bytes: Uint8Array) {
    this.view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  }
  get left() {
    return this.bytes.length - this.at;
  }
  need(n: number, what: string) {
    if (!(n >= 0) || n > this.left) throw new TabixError('malformed', `the index ends inside ${what}`);
  }
  i32(what: string) {
    this.need(4, what);
    const v = this.view.getInt32(this.at, true);
    this.at += 4;
    return v;
  }
  u32(what: string) {
    this.need(4, what);
    const v = this.view.getUint32(this.at, true);
    this.at += 4;
    return v;
  }
  voffset(what: string) {
    this.need(8, what);
    let v: number;
    try {
      v = VirtualOffset.read(this.view, this.at);
    } catch (e) {
      throw new TabixError('malformed', `${what}: ${(e as Error).message}`);
    }
    this.at += 8;
    return v;
  }
  /** A count of items of at least `size` bytes each, checked against the bytes left. */
  count(size: number, what: string, max = Infinity) {
    const n = this.i32(what);
    if (n < 0 || n > max || n * size > this.left) throw new TabixError('malformed', `${what} ${n} cannot be`);
    return n;
  }
}

const Formats = ['generic', 'sam', 'vcf'] as const;

/** The tabix header (TBI's, or a CSI's aux): format, columns, meta, skip, names. */
function readConf(c: Cursor, nRef: number | undefined) {
  const format = c.i32('the format');
  const type = format & 0xffff;
  if (type > 2 || format & ~0x1ffff)
    throw new TabixError('malformed', `an index format ${format} (SAM, VCF or generic)`);
  const seq = c.i32('the columns');
  const begin = c.i32('the columns');
  const end = c.i32('the columns');
  if (seq < 1 || begin < 1 || end < 0 || seq > 0xffff || begin > 0xffff || end > 0xffff) {
    throw new TabixError('malformed', `columns ${seq}, ${begin}, ${end} cannot be`);
  }
  const metaCode = c.i32('the meta character');
  if (metaCode < 0 || metaCode > 255) throw new TabixError('malformed', `a meta character ${metaCode}`);
  const skip = c.i32('the lines skipped');
  if (skip < 0) throw new TabixError('malformed', `${skip} lines skipped`);
  const lnm = c.i32('the names');
  c.need(lnm, 'the names');
  const raw = c.bytes.subarray(c.at, c.at + lnm);
  c.at += lnm;
  const names: string[] = [];
  const decoder = new TextDecoder();
  for (let s = 0, i = 0; i < raw.length; i++) {
    if (raw[i] === 0) {
      names.push(decoder.decode(raw.subarray(s, i)));
      s = i + 1;
    } else if (i === raw.length - 1) {
      throw new TabixError('malformed', 'the names do not end with a NUL');
    }
  }
  if (nRef !== undefined && names.length !== nRef) {
    throw new TabixError('malformed', `${names.length} names for ${nRef} references`);
  }
  return {
    format: Formats[type],
    zeroBased: (format & 0x10000) !== 0,
    columns: { seq, begin, end },
    meta: metaCode ? String.fromCharCode(metaCode) : '',
    skip,
    names,
  };
}

function readReferences(c: Cursor, nRef: number, kind: 'tbi' | 'csi', minShift: number, depth: number) {
  const refs: Reference[] = [];
  const bins = binCount(depth);
  const meta = metaBin(depth);
  const maxWindows = 2 ** (3 * depth);
  let firstLine: number | undefined;
  for (let r = 0; r < nRef; r++) {
    const nBin = c.count(kind === 'tbi' ? 8 : 16, `reference ${r}'s bins`, bins + 1);
    const map = new Map<number, { loffset: number; chunks: Float64Array }>();
    for (let b = 0; b < nBin; b++) {
      const bin = c.u32('a bin');
      if (bin >= bins && bin !== meta) throw new TabixError('malformed', `a bin ${bin} past the scheme's ${bins}`);
      const loffset = kind === 'csi' ? c.voffset('a bin offset') : 0;
      const nChunk = c.count(16, `bin ${bin}'s chunks`);
      const chunks = new Float64Array(2 * nChunk);
      for (let k = 0; k < nChunk; k++) {
        const begin = c.voffset('a chunk');
        const end = c.voffset('a chunk');
        // the pseudo-bin's second pair is counts (mapped, unmapped), not offsets
        if (begin > end && bin !== meta)
          throw new TabixError('malformed', `a chunk of bin ${bin} ends before it begins`);
        chunks[2 * k] = begin;
        chunks[2 * k + 1] = end;
      }
      if (bin === meta) continue; // the pseudo-bin: the reference's offsets and counts, not lines
      if (map.has(bin)) throw new TabixError('malformed', `bin ${bin} of reference ${r} listed twice`);
      map.set(bin, { loffset, chunks });
      for (let k = 0; k < nChunk; k++) {
        if (firstLine === undefined || chunks[2 * k] < firstLine) firstLine = chunks[2 * k];
      }
    }
    let linear = new Float64Array(0);
    if (kind === 'tbi') {
      const nIntv = c.count(8, `reference ${r}'s linear index`, maxWindows);
      linear = new Float64Array(nIntv);
      for (let k = 0; k < nIntv; k++) linear[k] = c.voffset('a linear index offset');
    }
    refs.push({ bins: map, linear });
  }
  return { refs, firstLine };
}

/** A `.tbi` or `.csi` (BGZF, told by its magic), every count and offset checked against `limits` and its own bytes. */
export async function readTabixIndex(
  bytes: Uint8Array,
  limits?: Partial<TabixLimits>,
  signal?: AbortSignal,
): Promise<TabixIndex> {
  const l = limitsOf(limits);
  let raw: Uint8Array;
  if (isBgzf(bytes)) {
    try {
      raw = await inflateBgzf(bytes, { maxInflatedBytes: l.maxIndexBytes }, signal);
    } catch (e) {
      if (e instanceof BgzfError) throw new TabixError(e.kind === 'too-large' ? 'too-large' : 'malformed', e.message);
      throw e;
    }
  } else {
    if (bytes.length > l.maxIndexBytes) {
      throw new TabixError('too-large', `an index of ${bytes.length} bytes, over the bound of ${l.maxIndexBytes}`);
    }
    raw = bytes;
  }
  signal?.throwIfAborted();
  const c = new Cursor(raw);
  c.need(4, 'the magic');
  const magic = String.fromCharCode(raw[0], raw[1], raw[2], raw[3]);
  let index: TabixIndex;
  let data: { refs: Reference[]; firstLine: number | undefined };
  c.at = 4;
  if (magic === 'TBI\x01') {
    const nRef = c.count(8, 'the references');
    const conf = readConf(c, nRef);
    data = readReferences(c, nRef, 'tbi', 14, 5);
    index = { kind: 'tbi', ...conf, minShift: 14, depth: 5 };
  } else if (magic === 'CSI\x01') {
    const minShift = c.i32('min_shift');
    const depth = c.i32('depth');
    if (minShift < 1 || minShift > 32 || depth < 1 || depth > 10 || minShift + 3 * depth > 52) {
      throw new TabixError('malformed', `a CSI scheme of min_shift ${minShift}, depth ${depth}`);
    }
    const lAux = c.i32('the aux data');
    c.need(lAux, 'the aux data');
    if (lAux < 28) throw new TabixError('malformed', 'a CSI without a tabix header (not a text file index)');
    const aux = new Cursor(raw.subarray(c.at, c.at + lAux));
    c.at += lAux;
    const conf = readConf(aux, undefined);
    const nRef = c.count(4, 'the references');
    if (conf.names.length !== nRef) {
      throw new TabixError('malformed', `${conf.names.length} names for ${nRef} references`);
    }
    data = readReferences(c, nRef, 'csi', minShift, depth);
    index = { kind: 'csi', ...conf, minShift, depth };
  } else {
    throw new TabixError('malformed', 'not a tabix index (no TBI or CSI magic)');
  }
  const byName = new Map<string, number>();
  index.names.forEach((n, i) => byName.set(n, i));
  Data.set(index, { refs: data.refs, byName, firstLine: data.firstLine });
  return index;
}

/** The smallest virtual offset a region's lines can start at (htslib's `min_off`). */
function minOffset(index: TabixIndex, ref: Reference, start: number): number {
  const window = Math.floor(start / 2 ** index.minShift);
  if (ref.linear.length > 0) {
    const n = ref.linear.length;
    let min = window >= n ? ref.linear[n - 1] : ref.linear[window];
    if (min === 0) {
      for (let i = Math.min(window, n) - 1; i >= 0; i--) {
        if (ref.linear[i] !== 0) {
          min = ref.linear[i];
          break;
        }
      }
    }
    return min;
  }
  if (index.kind !== 'csi') return 0;
  // CSI: the smallest bin holding `start` that the index has, up the tree: its `loffset`
  let bin = binFirst(index.depth) + window;
  while (bin > 0 && !ref.bins.has(bin)) bin = binParent(bin);
  return ref.bins.get(bin)?.loffset ?? 0;
}

/** The chunks that may hold lines meeting `[start, end)` (0-based, half-open): reg2bins, the linear index, merged. */
export function regionChunks(index: TabixIndex, name: string, start: number, end: number): readonly Chunk[] {
  const data = dataOf(index);
  const tid = data.byName.get(name);
  if (tid === undefined) return [];
  start = Math.max(0, start);
  if (!(end > start)) return [];
  const ref = data.refs[tid];
  const min = minOffset(index, ref, start);
  const found: Chunk[] = [];
  for (const bin of reg2bins(start, end, index.minShift, index.depth)) {
    const b = ref.bins.get(bin);
    if (!b) continue;
    for (let k = 0; k < b.chunks.length; k += 2) {
      if (b.chunks[k + 1] > min) found.push({ begin: Math.max(b.chunks[k], 0), end: b.chunks[k + 1] });
    }
  }
  found.sort((a, b) => a.begin - b.begin || a.end - b.end);
  // overlapping chunks, and chunks meeting in one block, read once (each line once)
  const merged: { begin: number; end: number }[] = [];
  for (const c of found) {
    const last = merged[merged.length - 1];
    if (last && (c.begin <= last.end || VirtualOffset.block(c.begin) === VirtualOffset.block(last.end))) {
      if (c.end > last.end) last.end = c.end;
    } else {
      merged.push({ begin: c.begin, end: c.end });
    }
  }
  return merged;
}

/** The compressed bytes of chunks: each from its first block's start to its last block's start (the index's offsets). */
function chunksBytes(chunks: readonly Chunk[]) {
  let total = 0;
  for (const c of chunks) total += VirtualOffset.block(c.end) - VirtualOffset.block(c.begin);
  return total;
}

/**
 * The compressed bytes of the chunks `regionChunks` gives, without reading them (a level's cost): each chunk from its
 * first block's start to its last block's start, as the index's offsets say (the last block's own bytes, at most
 * 64 KiB a chunk, are not in the index).
 */
export function regionBytes(index: TabixIndex, name: string, start: number, end: number): number {
  return chunksBytes(regionChunks(index, name, start, end));
}

function readBounded(reader: BgzfReader, begin: number, end: number, max: number, signal?: AbortSignal) {
  if (reader instanceof BgzfRangeReader) return reader.readWithin(begin, end, max, signal);
  return reader.read(begin, end, signal).then((b) => {
    if (b.length > max) throw new TabixError('too-large', `a read of ${b.length} bytes, over the bound of ${max}`);
    return b;
  });
}

/** Lines of `bytes` from `from` to `to` (`\n`, a `\r` before it dropped): starts and ends, each line within `max`. */
function splitLines(bytes: Uint8Array, from: number, to: number, max: number, starts: number[], ends: number[]) {
  let s = from;
  while (s < to) {
    let e = bytes.indexOf(10, s);
    if (e < 0 || e >= to) e = to;
    const next = e + 1;
    if (e > s && bytes[e - 1] === 13) e--;
    if (e - s > max) throw new TabixError('malformed', `a line of ${e - s} bytes, over the bound of ${max}`);
    starts.push(s);
    ends.push(e);
    s = next;
  }
}

/**
 * Lines meeting a region, as bytes, each with its extent by the index's rule; a line from a chunk's edge outside the
 * region is dropped. A line naming another sequence (an index for another file), or without a position, is
 * `malformed`; past `maxRegionBytes` (compressed, from the index) `too-large` before anything is read.
 */
export async function readRegion(
  reader: BgzfReader,
  index: TabixIndex,
  name: string,
  start: number,
  end: number,
  limits?: Partial<TabixLimits>,
  signal?: AbortSignal,
): Promise<TabixLines> {
  const l = limitsOf(limits);
  const chunks = regionChunks(index, name, start, end);
  const cost = chunksBytes(chunks);
  if (cost > l.maxRegionBytes) {
    throw new TabixError('too-large', `the region's chunks are ${cost} bytes, over the bound of ${l.maxRegionBytes}`);
  }
  const parts: Uint8Array[] = [];
  let total = 0;
  for (const c of chunks) {
    signal?.throwIfAborted();
    const got = await readBounded(reader, c.begin, c.end, l.maxInflatedBytes - total, signal);
    parts.push(got);
    total += got.length;
  }
  const bytes = parts.length === 1 ? parts[0] : join(parts, total);
  const starts: number[] = [];
  const ends: number[] = [];
  for (let i = 0, o = 0; i < parts.length; o += parts[i].length, i++) {
    splitLines(bytes, o, o + parts[i].length, l.maxLineBytes, starts, ends);
  }
  const nameBytes = new TextEncoder().encode(name);
  const meta = index.meta ? index.meta.charCodeAt(0) : -1;
  const lineStart = new Uint32Array(starts.length);
  const lineEnd = new Uint32Array(starts.length);
  const begin = new Int32Array(starts.length);
  const endOut = new Int32Array(starts.length);
  const x: Extent = { seqStart: 0, seqEnd: 0, begin: 0, end: 0 };
  let n = 0;
  for (let k = 0; k < starts.length; k++) {
    const s = starts[k];
    const e = ends[k];
    if (s === e || bytes[s] === meta) continue;
    if (!lineExtent(bytes, s, e, index, x)) {
      throw new TabixError('malformed', `a line of ${name}'s region has no sequence or position (byte ${s})`);
    }
    if (!sameBytes(bytes, x.seqStart, x.seqEnd, nameBytes)) {
      throw new TabixError('malformed', `a line of another sequence in ${name}'s chunks: the index is not this file's`);
    }
    if (x.begin > MaxPosition || x.end > MaxPosition) {
      throw new TabixError('malformed', `a line of ${name} past position 2^31 − 1`);
    }
    if (!(x.begin < end && x.end > start)) continue;
    lineStart[n] = s;
    lineEnd[n] = e;
    begin[n] = x.begin;
    endOut[n] = x.end;
    n++;
  }
  return {
    bytes,
    lineStart: lineStart.slice(0, n),
    lineEnd: lineEnd.slice(0, n),
    begin: begin.slice(0, n),
    end: endOut.slice(0, n),
  };
}

function join(parts: readonly Uint8Array[], total: number) {
  const out = new Uint8Array(total);
  for (let i = 0, o = 0; i < parts.length; o += parts[i].length, i++) out.set(parts[i], o);
  return out;
}

/**
 * The header lines (those starting with `meta`, and `skip` lines) at the top of the file: the bytes before the first
 * indexed line (or, with no line indexed, blocks from the start until a line that is not a header), within
 * `maxHeaderBytes`.
 */
export async function readHeaderLines(
  reader: BgzfReader,
  index: TabixIndex,
  limits?: Partial<TabixLimits>,
  signal?: AbortSignal,
): Promise<TabixLines> {
  const l = limitsOf(limits);
  const { firstLine } = dataOf(index);
  let bytes: Uint8Array;
  if (firstLine !== undefined) {
    bytes = await readBounded(reader, 0, firstLine, l.maxHeaderBytes, signal);
  } else if (reader instanceof BgzfRangeReader) {
    const parts: Uint8Array[] = [];
    let total = 0;
    for (let at: number | undefined = 0; at !== undefined; ) {
      const block = await reader.blockAt(at, signal);
      if (!block) break;
      if (total + block.bytes.length > l.maxHeaderBytes) {
        throw new TabixError('too-large', `the header passes the bound of ${l.maxHeaderBytes} bytes`);
      }
      parts.push(block.bytes);
      total += block.bytes.length;
      at = block.next;
      // stop at a block holding a line that is not a header line
      if (block.bytes.length && !allHeader(join(parts, total), index)) break;
    }
    bytes = join(parts, total);
  } else {
    bytes = new Uint8Array(0);
  }
  const starts: number[] = [];
  const ends: number[] = [];
  splitLines(bytes, 0, bytes.length, l.maxLineBytes, starts, ends);
  const meta = index.meta ? index.meta.charCodeAt(0) : -1;
  let n = 0;
  for (let k = 0; k < starts.length; k++) {
    const header = k < index.skip || (starts[k] < ends[k] && bytes[starts[k]] === meta);
    if (!header) break;
    starts[n] = starts[k];
    ends[n] = ends[k];
    n++;
  }
  const none = new Int32Array(n).fill(-1);
  return {
    bytes,
    lineStart: Uint32Array.from(starts.slice(0, n)),
    lineEnd: Uint32Array.from(ends.slice(0, n)),
    begin: none,
    end: none.slice(),
  };
}

/** Whether every complete line of `bytes` is a header line. */
function allHeader(bytes: Uint8Array, index: TabixIndex) {
  const meta = index.meta ? index.meta.charCodeAt(0) : -1;
  let line = 0;
  for (let s = 0; s < bytes.length; line++) {
    const e = bytes.indexOf(10, s);
    if (e < 0) return true;
    if (line >= index.skip && bytes[s] !== meta) return false;
    s = e + 1;
  }
  return true;
}
