/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { deflateRawSync } from 'node:zlib';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { crc } from '@molstar/core/util/zip/checksum';
import type { AssetRangeReader } from '@molstar/core/util/assets-range';
import { EofBlock, VirtualOffset } from '../internal/bgzf-blocks.js';
import { reg2bin } from '../internal/tabix-bins.js';

// A test-only BGZF and tabix writer: blocks as bgzip writes them, TBI and CSI as htslib builds them (bins with their
// chunks, the linear index, the metadata pseudo-bin). Fixtures are generated, never downloaded.

export interface Bgzipped {
  readonly bytes: Uint8Array;
  /** Each block's file offset and its first byte's offset in the text. */
  readonly blockOffsets: readonly number[];
  readonly blockStarts: readonly number[];
  /** A text offset's virtual offset (a block's end is the next block's start, as htslib tells it). */
  voffset(textOffset: number): VirtualOffset;
}

function block(chunk: Uint8Array, level = 6): Uint8Array {
  const deflated = deflateRawSync(chunk, { level });
  const out = new Uint8Array(18 + deflated.length + 8);
  const view = new DataView(out.buffer);
  out.set([0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, 66, 67, 2, 0]);
  view.setUint16(16, out.length - 1, true);
  out.set(deflated, 18);
  view.setUint32(18 + deflated.length, crc(chunk, 0, chunk.length) >>> 0, true);
  view.setUint32(22 + deflated.length, chunk.length, true);
  return out;
}

/** `text` in BGZF blocks of `blockSize` bytes (bgzip: 65,280), the EOF block last unless `eof` is false. */
export function bgzip(text: Uint8Array, blockSize = 65280, eof = true): Bgzipped {
  const parts: Uint8Array[] = [];
  const blockOffsets: number[] = [];
  const blockStarts: number[] = [];
  let at = 0;
  for (let s = 0; s < text.length; s += blockSize) {
    const b = block(text.subarray(s, s + blockSize));
    blockOffsets.push(at);
    blockStarts.push(s);
    parts.push(b);
    at += b.length;
  }
  const end = at;
  if (eof) parts.push(EofBlock);
  const bytes = join(parts);
  return {
    bytes,
    blockOffsets,
    blockStarts,
    voffset(p) {
      if (p >= text.length) return VirtualOffset.of(end, 0);
      let lo = 0;
      let hi = blockStarts.length - 1;
      while (lo < hi) {
        const mid = (lo + hi + 1) >> 1;
        if (blockStarts[mid] <= p) lo = mid;
        else hi = mid - 1;
      }
      return VirtualOffset.of(blockOffsets[lo], p - blockStarts[lo]);
    },
  };
}

export function join(parts: readonly Uint8Array[]): Uint8Array {
  const out = new Uint8Array(parts.reduce((n, p) => n + p.length, 0));
  let at = 0;
  for (const p of parts) {
    out.set(p, at);
    at += p.length;
  }
  return out;
}

/** A record's place in the text and its extent (0-based, half-open), as the index's rule gives it. */
export interface IndexedLine {
  readonly seq: number;
  readonly begin: number;
  readonly end: number;
  /** The line's first byte and the byte after its newline, in the text. */
  readonly start: number;
  readonly stop: number;
}

export interface IndexOptions {
  readonly kind: 'tbi' | 'csi';
  readonly format: 0 | 1 | 2;
  readonly zeroBased: boolean;
  readonly columns: readonly [number, number, number];
  readonly meta: string;
  readonly skip: number;
  readonly names: readonly string[];
  readonly minShift?: number;
  readonly depth?: number;
}

class LE {
  private parts: number[] = [];
  i32(v: number) {
    this.parts.push(v & 0xff, (v >>> 8) & 0xff, (v >>> 16) & 0xff, (v >>> 24) & 0xff);
  }
  u64(v: number) {
    this.i32(v % 0x100000000);
    this.i32(Math.floor(v / 0x100000000));
  }
  bytes(b: Uint8Array) {
    for (const x of b) this.parts.push(x);
  }
  done() {
    return new Uint8Array(this.parts);
  }
}

/** The index of `lines` (sorted by sequence, then begin) over the BGZF file `file`, inflated (raw TBI or CSI bytes). */
export function tabixIndexRaw(lines: readonly IndexedLine[], file: Bgzipped, o: IndexOptions): Uint8Array {
  const minShift = o.minShift ?? 14;
  const depth = o.depth ?? 5;
  const meta = (8 ** (depth + 1) - 1) / 7 + 1;
  const refs = o.names.map(() => ({
    bins: new Map<number, number[]>(),
    linear: [] as number[],
    first: -1,
    last: -1,
    count: 0,
  }));
  for (const l of lines) {
    const r = refs[l.seq];
    const begin = file.voffset(l.start);
    const end = file.voffset(l.stop);
    const bin = reg2bin(l.begin, Math.max(l.end, l.begin + 1), minShift, depth);
    let chunks = r.bins.get(bin);
    if (!chunks) r.bins.set(bin, (chunks = []));
    if (chunks.length && chunks[chunks.length - 1] === begin) chunks[chunks.length - 1] = end;
    else chunks.push(begin, end);
    const w0 = Math.floor(l.begin / 2 ** minShift);
    const w1 = Math.floor((Math.max(l.end, l.begin + 1) - 1) / 2 ** minShift);
    for (let w = w0; w <= w1; w++) if (r.linear[w] === undefined) r.linear[w] = begin;
    if (r.first < 0) r.first = begin;
    r.last = end;
    r.count++;
  }
  for (const r of refs) {
    // htslib's update_loff: leading windows take the reference's first offset, the others the previous window's
    for (let w = 0; w < r.linear.length; w++) {
      if (r.linear[w] === undefined) r.linear[w] = w === 0 ? r.first : r.linear[w - 1];
    }
  }
  const conf = new LE();
  conf.i32(o.format | (o.zeroBased ? 0x10000 : 0));
  for (const c of o.columns) conf.i32(c);
  conf.i32(o.meta ? o.meta.charCodeAt(0) : 0);
  conf.i32(o.skip);
  const names = join(o.names.map((n) => join([new TextEncoder().encode(n), new Uint8Array(1)])));
  conf.i32(names.length);
  conf.bytes(names);
  const confBytes = conf.done();

  const w = new LE();
  if (o.kind === 'tbi') {
    w.bytes(new TextEncoder().encode('TBI\x01'));
    w.i32(o.names.length);
    w.bytes(confBytes);
  } else {
    w.bytes(new TextEncoder().encode('CSI\x01'));
    w.i32(minShift);
    w.i32(depth);
    w.i32(confBytes.length);
    w.bytes(confBytes);
    w.i32(o.names.length);
  }
  for (const r of refs) {
    const bins = [...r.bins.keys()].sort((a, b) => a - b);
    w.i32(bins.length + (r.count ? 1 : 0));
    for (const bin of bins) {
      const chunks = r.bins.get(bin)!;
      w.i32(bin);
      if (o.kind === 'csi') {
        // htslib: a bin's loffset is the linear index at its first window
        let level = 0;
        while ((8 ** (level + 1) - 1) / 7 <= bin) level++;
        const bot = (bin - (8 ** level - 1) / 7) * 8 ** (depth - level);
        w.u64(bot < r.linear.length ? r.linear[bot] : 0);
      }
      w.i32(chunks.length / 2);
      for (const c of chunks) w.u64(c);
    }
    if (r.count) {
      w.i32(meta);
      if (o.kind === 'csi') w.u64(0);
      w.i32(2);
      w.u64(r.first);
      w.u64(r.last);
      w.u64(r.count);
      w.u64(0);
    }
    if (o.kind === 'tbi') {
      w.i32(r.linear.length);
      for (const v of r.linear) w.u64(v);
    }
  }
  w.u64(0); // n_no_coor
  return w.done();
}

/** The index of `lines`, bgzipped as tabix writes it. */
export function tabixIndex(lines: readonly IndexedLine[], file: Bgzipped, o: IndexOptions): Uint8Array {
  return bgzip(tabixIndexRaw(lines, file, o)).bytes;
}

/** A range reader over bytes in memory, counting its reads. */
export function memoryRanges(bytes: Uint8Array): AssetRangeReader & { reads: number } {
  const r = {
    reads: 0,
    size: bytes.length,
    read(ranges: readonly { offset: number; length: number }[]) {
      r.reads++;
      return Promise.resolve(ranges.map((g) => bytes.subarray(g.offset, Math.min(g.offset + g.length, bytes.length))));
    },
    dispose() {},
  };
  return r;
}

export interface GeneratedFile {
  readonly text: Uint8Array;
  readonly header: readonly string[];
  readonly lines: readonly (IndexedLine & { readonly text: string })[];
  readonly names: readonly string[];
}

/**
 * A BED-like file (0-based, `end` exclusive; columns 1, 2, 3; `#` header) of `n` records over `names`, sorted, seeded;
 * lengths up to `maxLength`, positions from `offset` (a CSI's reference past 512 Mb).
 */
export function generateBed(
  n: number,
  names: readonly string[],
  seed = 1,
  maxLength = 3000,
  offset = 0,
): GeneratedFile {
  const rng = new PCG(seed);
  const header = ['#track name=generated', '#chrom\tstart\tend\tname'];
  const rows: { seq: number; begin: number; end: number }[] = [];
  for (let i = 0; i < n; i++) {
    const seq = Math.floor(rng.float() * names.length);
    const begin = offset + Math.floor(rng.float() * 2_000_000);
    rows.push({ seq, begin, end: begin + 1 + Math.floor(rng.float() ** 4 * maxLength) });
  }
  rows.sort((a, b) => a.seq - b.seq || a.begin - b.begin);
  return textOf(
    header,
    rows.map((r, i) => ({ ...r, text: `${names[r.seq]}\t${r.begin}\t${r.end}\tfeature_${i}` })),
    names,
  );
}

/**
 * A VCF-like file (columns 1, 2, 0; `#` header) of `n` records: SNVs, deletions (their REF's length), and some with
 * `END=` in INFO (an SV's, and one on an SNV line: the index takes it whatever the alleles).
 */
export function generateVcf(n: number, names: readonly string[], seed = 2): GeneratedFile {
  const rng = new PCG(seed);
  const header = [
    '##fileformat=VCFv4.2',
    '##INFO=<ID=END,Number=1,Type=Integer>',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
  ];
  const rows: { seq: number; begin: number; end: number; text?: string }[] = [];
  for (let i = 0; i < n; i++) {
    const seq = Math.floor(rng.float() * names.length);
    const pos = 1 + Math.floor(rng.float() * 1_000_000);
    const kind = rng.float();
    let ref = 'A';
    let info = 'DP=10';
    let end = pos - 1 + 1;
    if (kind < 0.2) {
      ref = 'ACGTACGTAC'.slice(0, 2 + Math.floor(rng.float() * 8));
      end = pos - 1 + ref.length;
    } else if (kind < 0.3) {
      const svEnd = pos + Math.floor(rng.float() * 50_000);
      info = `SVTYPE=DEL;END=${svEnd}`;
      end = Math.max(svEnd, pos);
      if (svEnd <= pos - 1) end = pos;
    } else if (kind < 0.35) {
      // END= at INFO's start, on an SNV line
      const e = pos + 100;
      info = `END=${e};DP=3`;
      end = e;
    } else if (kind < 0.4) {
      // SVEND= is not END=
      info = `SVEND=${pos + 1000}`;
    }
    rows.push({ seq, begin: pos - 1, end, text: `${pos}\tid${i}\t${ref}\tT\t50\tPASS\t${info}` });
  }
  rows.sort((a, b) => a.seq - b.seq || a.begin - b.begin);
  return textOf(
    header,
    rows.map((r) => ({ seq: r.seq, begin: r.begin, end: r.end, text: `${names[r.seq]}\t${r.text}` })),
    names,
  );
}

function textOf(
  header: readonly string[],
  rows: readonly { seq: number; begin: number; end: number; text: string }[],
  names: readonly string[],
): GeneratedFile {
  const enc = new TextEncoder();
  const parts: Uint8Array[] = [];
  let at = 0;
  for (const h of header) {
    const b = enc.encode(`${h}\n`);
    parts.push(b);
    at += b.length;
  }
  const lines = rows.map((r) => {
    const b = enc.encode(`${r.text}\n`);
    parts.push(b);
    const line = { ...r, start: at, stop: at + b.length };
    at += b.length;
    return line;
  });
  return { text: join(parts), header, lines, names };
}
