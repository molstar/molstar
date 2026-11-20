/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import type { AssetRangeReader } from '@molstar/core/util/assets';
import { SourceIssues } from '../data/issues.js';
import {
  abortError,
  type AssetReader,
  type Coverage,
  type FormatRuntime,
  isAbortError,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type RegionRequest,
  type SourceEntry,
  type ValuesBins,
  type ValuesPositions,
} from '../data/protocol.js';
import { SourceOpenError } from '../data/seq-runtime.js';

// bigWig (_spec/README.md "Data", bigWig), the opt-in runtime `@molstar/seqx-core/runtimes/bigwig`: a source of values
// read by region and zoom level through range reads of its asset (core's `openRanges`, `util/_spec/assets.md` §10).
// Opening reads the header, the zoom headers and the chromosome tree; a region reads the R-tree index of its level
// (the full data at level 1, a zoom level's summaries at its reduction), then the data blocks it points to, inflated
// and decoded. Levels are the file's own: 1 and its zoom levels' reductions. Little- and big-endian files alike.

/** The assets port with range reads: `openRange` when the host has it (core's `openRanges` under its policy). */
export interface RangeAssetReader extends AssetReader {
  openRange?(uri: string, signal: AbortSignal): Promise<AssetRangeReader>;
}

const BigWigMagic = 0x888ffc26;
const ChromTreeMagic = 0x78ca8c91;
const RTreeMagic = 0x2468ace0;
/** The deepest index tree read (a file of a few billion blocks is about 6). */
const MaxDepth = 16;
/** Decoded blocks a source keeps for a pan back (least recently used go first). */
const BlockCache = 16;

class Malformed extends Error {}

interface Zoom {
  readonly level: number;
  readonly index: number;
}

interface Chrom {
  readonly id: number;
  readonly name: string;
  readonly length: number;
}

interface Block {
  readonly offset: number;
  readonly size: number;
}

/** An unsigned 64-bit integer as a number (exact to 2^53, a file of 8 PB). */
function u64(dv: DataView, o: number, le: boolean): number {
  const lo = dv.getUint32(o + (le ? 0 : 4), le);
  const hi = dv.getUint32(o + (le ? 4 : 0), le);
  return hi * 0x100000000 + lo;
}

const viewOf = (b: Uint8Array) => new DataView(b.buffer, b.byteOffset, b.byteLength);

/** A reader of whole bytes: a host without range reads (a test's, a CLI's `{ reader }`). */
function wholeReader(bytes: Uint8Array): AssetRangeReader {
  return {
    size: bytes.length,
    read: async (ranges) =>
      ranges.map((r) => bytes.subarray(Math.min(r.offset, bytes.length), Math.min(bytes.length, r.offset + r.length))),
    dispose() {},
  };
}

/** Inflates a zlib block, at most `limit` bytes (the header's `uncompressBufSize`). */
async function inflate(bytes: Uint8Array, limit: number): Promise<Uint8Array> {
  const stream = new Blob([bytes as Uint8Array<ArrayBuffer>]).stream().pipeThrough(new DecompressionStream('deflate'));
  const reader = stream.getReader();
  const out = new Uint8Array(limit);
  let n = 0;
  for (;;) {
    let next: ReadableStreamReadResult<Uint8Array>;
    try {
      next = await reader.read();
    } catch {
      throw new Malformed('a data block does not inflate');
    }
    if (next.done) break;
    if (n + next.value.length > limit) {
      await reader.cancel().catch(() => {});
      throw new Malformed(`a data block inflates past ${limit} bytes`);
    }
    out.set(next.value, n);
    n += next.value.length;
  }
  return out.subarray(0, n);
}

/** An opened bigWig: its header, chromosomes and the reads of a region. */
class BigWig {
  /** Decoded blocks by offset, least recently used first; only what arrived (a read in flight is the reader's to share). */
  private readonly blocks = new Map<number, Uint8Array>();

  constructor(
    readonly reader: AssetRangeReader,
    readonly le: boolean,
    readonly zooms: readonly Zoom[],
    readonly fullIndex: number,
    readonly uncompress: number,
    readonly chroms: ReadonlyMap<string, Chrom>,
  ) {}

  /** The bytes of `ranges`, each the length asked (else the file is cut short). */
  async read(ranges: readonly Block[], signal: AbortSignal, what: string): Promise<Uint8Array[]> {
    const out = await this.reader.read(
      ranges.map((r) => ({ offset: r.offset, length: r.size })),
      signal,
    );
    for (let i = 0; i < out.length; i++)
      if (out[i].length < ranges[i].size) throw new Malformed(`${what} runs past the end of the file`);
    return out;
  }

  /** The blocks of an R-tree that overlap `[start, end)` of `chrom`, in file order. */
  async query(index: number, chrom: number, start: number, end: number, signal: AbortSignal): Promise<Block[]> {
    const le = this.le;
    const [head] = await this.read([{ offset: index, size: 48 }], signal, 'an index');
    const dv = viewOf(head);
    if (dv.getUint32(0, le) !== RTreeMagic) throw new Malformed('an index is not an R-tree');
    const blockSize = dv.getUint32(4, le);
    if (blockSize < 1 || blockSize > 1 << 16) throw new Malformed(`an index's block size ${blockSize}`);
    const found: Block[] = [];
    let level = [index + 48];
    for (let depth = 0; level.length; depth++) {
      if (depth > MaxDepth) throw new Malformed('an index deeper than any file needs');
      // a node's items: 32 bytes a leaf's, 24 an inner one's; the longer read, clipped at the end of the file
      const nodes = await this.reader.read(
        level.map((offset) => ({ offset, length: 4 + 32 * blockSize })),
        signal,
      );
      const next: number[] = [];
      for (const node of nodes) {
        if (node.length < 4) throw new Malformed('an index node runs past the end of the file');
        const v = viewOf(node);
        const leaf = v.getUint8(0) === 1;
        const count = v.getUint16(2, le);
        const size = leaf ? 32 : 24;
        if (count > blockSize || 4 + count * size > node.length) throw new Malformed('an index node past its size');
        for (let i = 0, o = 4; i < count; i++, o += size) {
          const sc = v.getUint32(o, le);
          const sb = v.getUint32(o + 4, le);
          const ec = v.getUint32(o + 8, le);
          const eb = v.getUint32(o + 12, le);
          const after = sc > chrom || (sc === chrom && sb >= end);
          const before = ec < chrom || (ec === chrom && eb <= start);
          if (after || before) continue;
          if (leaf) found.push({ offset: u64(v, o + 16, le), size: u64(v, o + 24, le) });
          else next.push(u64(v, o + 16, le));
        }
      }
      level = next;
    }
    return found.sort((a, b) => a.offset - b.offset);
  }

  /** The data blocks inflated, from the source's cache or read (the missing ones asked at once). */
  async blocksOf(found: readonly Block[], signal: AbortSignal): Promise<Uint8Array[]> {
    const missing = found.filter((b) => !this.blocks.has(b.offset));
    const raw = missing.length ? await this.read(missing, signal, 'a data block') : [];
    const inflated = await Promise.all(raw.map((r) => (this.uncompress > 0 ? inflate(r, this.uncompress) : r)));
    const got = new Map(missing.map((b, i) => [b.offset, inflated[i]]));
    const out = found.map((b) => got.get(b.offset) ?? this.blocks.get(b.offset)!);
    found.forEach((b, i) => {
      // most recently used last
      this.blocks.delete(b.offset);
      this.blocks.set(b.offset, out[i]);
      if (this.blocks.size > BlockCache) this.blocks.delete(this.blocks.keys().next().value!);
    });
    return out;
  }

  /** Per-position values of `[start, end)` from the full data's sections (bedGraph, varStep, fixedStep). */
  positions(blocks: readonly Uint8Array[], chrom: number, start: number, end: number): ValuesPositions {
    const le = this.le;
    const value = new Float32Array(end - start).fill(Number.NaN);
    const put = (a: number, b: number, x: number) => {
      for (let p = Math.max(a, start), e = Math.min(b, end); p < e; p++) value[p - start] = x;
    };
    for (const block of blocks) {
      const v = viewOf(block);
      if (block.length < 24) throw new Malformed('a data section shorter than its header');
      if (v.getUint32(0, le) !== chrom) continue;
      const first = v.getUint32(4, le);
      const step = v.getUint32(12, le);
      const span = v.getUint32(16, le);
      const type = v.getUint8(20);
      const count = v.getUint16(22, le);
      const size = type === 1 ? 12 : type === 2 ? 8 : type === 3 ? 4 : 0;
      if (!size) throw new Malformed(`a data section of type ${type}`);
      if (24 + count * size > block.length) throw new Malformed('a data section past its block');
      for (let i = 0, o = 24; i < count; i++, o += size) {
        if (type === 1) put(v.getUint32(o, le), v.getUint32(o + 4, le), v.getFloat32(o + 8, le));
        else if (type === 2) {
          const s = v.getUint32(o, le);
          put(s, s + span, v.getFloat32(o + 4, le));
        } else {
          const s = first + i * step;
          put(s, s + span, v.getFloat32(o, le));
        }
      }
    }
    return { form: 'positions', start, value };
  }

  /**
   * Summaries per bin of `level` from a zoom level's records: a record over several bins (a file whose summaries do
   * not start on the level's multiples) is shared by its overlap with each, its min and max given to each.
   */
  bins(blocks: readonly Uint8Array[], chrom: number, level: number, from: number, n: number): ValuesBins {
    const le = this.le;
    const min = new Float32Array(n).fill(Number.NaN);
    const max = new Float32Array(n).fill(Number.NaN);
    const mean = new Float32Array(n).fill(Number.NaN);
    const count = new Float32Array(n);
    const sum = new Float64Array(n);
    const to = from + n * level;
    for (const block of blocks) {
      if (block.length % 32) throw new Malformed('a zoom block is not whole records');
      const v = viewOf(block);
      for (let o = 0; o < block.length; o += 32) {
        if (v.getUint32(o, le) !== chrom) continue;
        const s = v.getUint32(o + 4, le);
        const e = v.getUint32(o + 8, le);
        const c = v.getUint32(o + 12, le);
        if (e <= s || c === 0 || e <= from || s >= to) continue;
        const lo = v.getFloat32(o + 16, le);
        const hi = v.getFloat32(o + 20, le);
        const total = v.getFloat32(o + 24, le);
        for (let b = Math.max(0, Math.floor((s - from) / level)); b < n && from + b * level < e; b++) {
          const share = (Math.min(e, from + (b + 1) * level) - Math.max(s, from + b * level)) / (e - s);
          if (share <= 0) continue;
          if (count[b] === 0 || lo < min[b]) min[b] = lo;
          if (count[b] === 0 || hi > max[b]) max[b] = hi;
          count[b] += c * share;
          sum[b] += total * share;
        }
      }
    }
    for (let b = 0; b < n; b++) if (count[b] > 0) mean[b] = sum[b] / count[b];
    return { form: 'bins', start: from, min, max, mean, count };
  }
}

/** The header, zoom headers and chromosome tree of a bigWig. */
async function readBigWig(reader: AssetRangeReader, signal: AbortSignal, maxRows: number): Promise<BigWig> {
  const [head] = await reader.read([{ offset: 0, length: 64 }], signal);
  if (head.length < 64) throw new Malformed('shorter than a bigWig header');
  const h = viewOf(head);
  const le = h.getUint32(0, true) === BigWigMagic;
  if (!le && h.getUint32(0, false) !== BigWigMagic) throw new Malformed('not a bigWig (its magic number)');
  const zoomCount = h.getUint16(6, le);
  const chromTree = u64(h, 8, le);
  const fullIndex = u64(h, 24, le);
  const uncompress = h.getUint32(52, le);
  if (uncompress > 1 << 26) throw new Malformed(`a block of ${uncompress} bytes inflated`);

  const [zh] = await reader.read([{ offset: 64, length: 24 * zoomCount }], signal);
  if (zh.length < 24 * zoomCount) throw new Malformed('the zoom headers run past the end of the file');
  const z = viewOf(zh);
  const zooms: Zoom[] = [];
  for (let i = 0; i < zoomCount; i++) {
    const level = z.getUint32(24 * i, le);
    if (level > 1 && !zooms.some((x) => x.level === level)) zooms.push({ level, index: u64(z, 24 * i + 16, le) });
  }
  zooms.sort((a, b) => a.level - b.level);

  const [th] = await reader.read([{ offset: chromTree, length: 32 }], signal);
  if (th.length < 32) throw new Malformed('the chromosome tree runs past the end of the file');
  const t = viewOf(th);
  if (t.getUint32(0, le) !== ChromTreeMagic) throw new Malformed('the chromosome tree (its magic number)');
  const blockSize = t.getUint32(4, le);
  const keySize = t.getUint32(8, le);
  const items = u64(t, 16, le);
  if (blockSize < 1 || blockSize > 1 << 16 || keySize < 1 || keySize > 1024 || items > maxRows)
    throw new Malformed('the chromosome tree’s sizes');
  const chroms = new Map<string, Chrom>();
  const latin1 = new TextDecoder('latin1');
  let level = [chromTree + 32];
  for (let depth = 0; level.length; depth++) {
    if (depth > MaxDepth) throw new Malformed('a chromosome tree deeper than any file needs');
    const nodes = await reader.read(
      level.map((offset) => ({ offset, length: 4 + blockSize * (keySize + 8) })),
      signal,
    );
    const next: number[] = [];
    for (const node of nodes) {
      const v = viewOf(node);
      if (node.length < 4) throw new Malformed('a chromosome tree node past the end of the file');
      const leaf = v.getUint8(0) === 1;
      const count = v.getUint16(2, le);
      if (count > blockSize || 4 + count * (keySize + 8) > node.length)
        throw new Malformed('a chromosome tree node past its size');
      for (let i = 0, o = 4; i < count; i++, o += keySize + 8) {
        if (!leaf) {
          next.push(u64(v, o + keySize, le));
          continue;
        }
        const key = node.subarray(o, o + keySize);
        const zero = key.indexOf(0);
        const name = latin1.decode(zero < 0 ? key : key.subarray(0, zero));
        if (chroms.has(name)) throw new Malformed(`the chromosome ${name} is listed twice`);
        chroms.set(name, { id: v.getUint32(o + keySize, le), name, length: v.getUint32(o + keySize + 4, le) });
        if (chroms.size > items) throw new Malformed('more chromosomes than the tree says');
      }
    }
    level = next;
  }
  return new BigWig(reader, le, zooms, fullIndex, uncompress, chroms);
}

function failure(kind: 'not-served' | 'malformed' | 'too-large', why: string): RegionAnswer {
  const code =
    kind === 'not-served'
      ? 'seq.source.not_served'
      : kind === 'malformed'
        ? 'seq.source.malformed'
        : 'seq.source.too_large';
  return { ok: false, kind, issue: issueOf(SourceIssues, code, `The bigWig source: ${why}`) };
}

/** An opened bigWig as a source of `values` by region: level 1 per position, a zoom level per bin of its width. */
export async function openBigWig(
  reader: AssetRangeReader,
  signal: AbortSignal,
  limits: OpenContext['limits'],
  source: SourceEntry = { address: '', format: 'bigwig' },
): Promise<OpenedSource> {
  let file: BigWig;
  try {
    file = await readBigWig(reader, signal, limits.maxRows);
  } catch (e) {
    reader.dispose();
    if (signal.aborted || isAbortError(e)) throw abortError();
    const malformed = e instanceof Malformed;
    throw new SourceOpenError([
      issueOf(
        SourceIssues,
        malformed ? 'seq.source.malformed' : 'seq.source.unreadable',
        `The bigWig source cannot be read: ${String((e as Error)?.message ?? e)}`,
        { address: source.address },
      ),
    ]);
  }
  const chroms = [...file.chroms.values()];
  const coverage: Coverage = {
    systems: chroms.map((c) => ({ name: c.name, length: c.length })),
    shapes: chroms.map((c) => ({ system: c.name, shape: 'values' as const })),
    levels: [1, ...file.zooms.map((z) => z.level)],
  };

  async function read(request: RegionRequest, signal: AbortSignal): Promise<RegionAnswer> {
    if (signal.aborted) throw abortError();
    const chrom = file.chroms.get(request.system);
    if (!chrom) return failure('not-served', `no chromosome ${request.system}`);
    if (!request.shapes.includes('values')) return failure('not-served', 'it serves values only');
    const L = request.level;
    const zoom = L === 1 ? undefined : file.zooms.find((z) => z.level === L);
    if (L !== 1 && !zoom) return failure('not-served', `no zoom level of ${L}`);
    const start = Math.max(0, Math.min(request.start, chrom.length));
    const end = Math.max(start, Math.min(request.end, chrom.length));
    // at a zoom level: bins of its width from a multiple of it, the last one cut at the chromosome's end
    const from = Math.floor(start / L) * L;
    const n = end > start ? Math.ceil(end / L) - from / L : 0;
    if (n > limits.maxRows) return failure('too-large', `${n} values asked at once`);
    try {
      let blocks: Uint8Array[] = [];
      if (n > 0) {
        const index = zoom ? zoom.index : file.fullIndex;
        const found = await file.query(index, chrom.id, from, from + n * L, signal);
        let bytes = 0;
        for (const b of found) bytes += b.size;
        if (bytes > limits.maxBytes) return failure('too-large', `${bytes} bytes of blocks for one region`);
        blocks = await file.blocksOf(found, signal);
      }
      if (signal.aborted) throw abortError();
      const values = zoom ? file.bins(blocks, chrom.id, L, from, n) : file.positions(blocks, chrom.id, start, end);
      return { ok: true, tables: { values } };
    } catch (e) {
      if (signal.aborted || isAbortError(e)) throw abortError();
      if (e instanceof Malformed) return failure('malformed', e.message);
      throw e;
    }
  }

  return { coverage, read, dispose: () => reader.dispose() };
}

async function open(source: SourceEntry, ctx: OpenContext): Promise<OpenedSource> {
  const uri = source.uri;
  if (uri === undefined)
    throw new SourceOpenError([
      issueOf(SourceIssues, 'seq.source.no_data', 'A bigWig source needs a uri', { address: source.address }),
    ]);
  const assets = ctx.assets as RangeAssetReader;
  let reader: AssetRangeReader;
  try {
    if (assets.openRange) reader = await assets.openRange(uri, ctx.signal);
    else {
      // a host without range reads: the file whole, bounded as any file parsed whole
      const bytes = await assets.readBytes(uri, ctx.signal);
      if (bytes.length > ctx.limits.maxBytes) throw new Error(`larger than ${ctx.limits.maxBytes} bytes`);
      reader = wholeReader(bytes);
    }
  } catch (e) {
    if (isAbortError(e) || ctx.signal.aborted) throw abortError();
    throw new SourceOpenError([
      issueOf(SourceIssues, 'seq.source.unreadable', `The file ${uri} could not be read`, { address: source.address }),
    ]);
  }
  return openBigWig(reader, ctx.signal, ctx.limits, source);
}

/** The runtime of bigWig (`format: 'bigwig'`): opt in by listing it in a viewer's registry. */
export const BigWigRuntime: FormatRuntime = { name: 'bigwig', formats: ['bigwig'], open };
