/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { crc } from '@molstar/core/util/zip/checksum';
import { deflate } from '@molstar/core/util/zip/zip';

// BGZF and tabix as bgzip and `tabix` write them, for the Genomes dataset's generated files (a VCF, a bedGraph, a GFF3
// and a GTF, each with its `.tbi`): blocks of 65,280 bytes deflated by core's deflate with `native: false` (its pure
// JavaScript, so the bytes are the same on every platform and Node version: what a byte-identity test needs), the EOF
// block last; the TBI index of htslib (bins with their chunks, the linear index of 16 kb windows, the metadata
// pseudo-bin), itself bgzipped. Why not the test harness's writer: it deflates with Node's zlib, whose output may
// change with Node's version, and the examples' sources do not depend on the harness.

/** The 28-byte EOF block bgzip writes last. */
const EofBlock = new Uint8Array([
  0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0,
]);

/** bgzip's block size: 65,280 bytes of text a block. */
const BlockSize = 65280;

/** A text in BGZF blocks, and the virtual offset of a text offset (a block's file offset × 2^16 + an offset in it). */
export interface Bgzipped {
  readonly bytes: Uint8Array;
  voffset(textOffset: number): number;
}

/** One BGZF block of `chunk`: a gzip member with the `BC` extra field, deflated raw (core's zlib stream unwrapped). */
async function block(chunk: Uint8Array): Promise<Uint8Array> {
  // core's own deflate, never the platform's (`native: false`): the platform's bytes change with its zlib
  const zlib = await deflate(SyncRuntimeContext, chunk as Uint8Array<ArrayBuffer>, { level: 6, native: false });
  // a zlib stream: a 2-byte header, the raw deflate, a 4-byte Adler-32
  const raw = zlib.subarray(2, zlib.length - 4);
  const out = new Uint8Array(18 + raw.length + 8);
  const view = new DataView(out.buffer);
  out.set([0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 6, 0, 66, 67, 2, 0]);
  view.setUint16(16, out.length - 1, true);
  out.set(raw, 18);
  view.setUint32(18 + raw.length, crc(chunk, 0, chunk.length) >>> 0, true);
  view.setUint32(22 + raw.length, chunk.length, true);
  return out;
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

/** `text` in BGZF blocks, the EOF block last. */
export async function bgzip(text: Uint8Array): Promise<Bgzipped> {
  const parts: Uint8Array[] = [];
  const offsets: number[] = [];
  const starts: number[] = [];
  let at = 0;
  for (let s = 0; s < text.length; s += BlockSize) {
    const b = await block(text.subarray(s, s + BlockSize));
    offsets.push(at);
    starts.push(s);
    parts.push(b);
    at += b.length;
  }
  const end = at;
  parts.push(EofBlock);
  return {
    bytes: join(parts),
    voffset(p) {
      if (p >= text.length) return end * 0x10000;
      let lo = 0;
      let hi = starts.length - 1;
      while (lo < hi) {
        const mid = (lo + hi + 1) >> 1;
        if (starts[mid] <= p) lo = mid;
        else hi = mid - 1;
      }
      return offsets[lo] * 0x10000 + (p - starts[lo]);
    },
  };
}

/** A line of an indexed file: its sequence's index, its extent (0-based, half-open), its bytes in the text. */
export interface IndexedLine {
  readonly seq: number;
  readonly begin: number;
  readonly end: number;
  /** The line's first byte and the byte after its newline. */
  readonly start: number;
  readonly stop: number;
}

/** A TBI's configuration: tabix's preset (`vcf` 2, generic 0), its columns, comment character and names. */
export interface TabixOptions {
  readonly format: 0 | 2;
  readonly zeroBased: boolean;
  /** The sequence, begin and end columns (1-based; 0 for none). */
  readonly columns: readonly [number, number, number];
  readonly meta: string;
  readonly names: readonly string[];
}

/** The smallest bin holding `[beg, end)` (htslib's `hts_reg2bin`, 14 bits and 5 levels). */
function reg2bin(beg: number, end: number): number {
  end -= 1;
  if (beg >> 14 === end >> 14) return 4681 + (beg >> 14);
  if (beg >> 17 === end >> 17) return 585 + (beg >> 17);
  if (beg >> 20 === end >> 20) return 73 + (beg >> 20);
  if (beg >> 23 === end >> 23) return 9 + (beg >> 23);
  if (beg >> 26 === end >> 26) return 1 + (beg >> 26);
  return 0;
}

class Le {
  private readonly parts: number[] = [];
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

/** The TBI index of `lines` (sorted by sequence, then begin) over `file`, bgzipped as tabix writes it. */
export async function tabix(lines: readonly IndexedLine[], file: Bgzipped, o: TabixOptions): Promise<Uint8Array> {
  const refs = o.names.map(() => ({
    bins: new Map<number, number[]>(),
    linear: [] as number[],
    first: -1,
    last: -1,
    n: 0,
  }));
  for (const l of lines) {
    const r = refs[l.seq];
    const begin = file.voffset(l.start);
    const end = file.voffset(l.stop);
    const e = Math.max(l.end, l.begin + 1);
    const bin = reg2bin(l.begin, e);
    let chunks = r.bins.get(bin);
    if (!chunks) r.bins.set(bin, (chunks = []));
    if (chunks.length && chunks[chunks.length - 1] === begin) chunks[chunks.length - 1] = end;
    else chunks.push(begin, end);
    for (let w = l.begin >> 14; w <= (e - 1) >> 14; w++) if (r.linear[w] === undefined) r.linear[w] = begin;
    if (r.first < 0) r.first = begin;
    r.last = end;
    r.n++;
  }
  // htslib's update_loff: leading windows take the reference's first offset, the others the previous window's
  for (const r of refs)
    for (let w = 0; w < r.linear.length; w++)
      if (r.linear[w] === undefined) r.linear[w] = w ? r.linear[w - 1] : r.first;
  const enc = new TextEncoder();
  const w = new Le();
  w.bytes(enc.encode('TBI\x01'));
  w.i32(o.names.length);
  w.i32(o.format | (o.zeroBased ? 0x10000 : 0));
  for (const c of o.columns) w.i32(c);
  w.i32(o.meta.charCodeAt(0));
  w.i32(0);
  const names = join(o.names.map((n) => join([enc.encode(n), new Uint8Array(1)])));
  w.i32(names.length);
  w.bytes(names);
  for (const r of refs) {
    const bins = [...r.bins.keys()].sort((a, b) => a - b);
    w.i32(bins.length + (r.n ? 1 : 0));
    for (const bin of bins) {
      const chunks = r.bins.get(bin)!;
      w.i32(bin);
      w.i32(chunks.length / 2);
      for (const c of chunks) w.u64(c);
    }
    if (r.n) {
      // the metadata pseudo-bin: the reference's first and last offsets, its mapped and unmapped counts
      w.i32(37450);
      w.i32(2);
      w.u64(r.first);
      w.u64(r.last);
      w.u64(r.n);
      w.u64(0);
    }
    w.i32(r.linear.length);
    for (const v of r.linear) w.u64(v);
  }
  w.u64(0);
  return (await bgzip(w.done())).bytes;
}

/** A sorted text's lines as the index reads them: `extent` gives each line's sequence and extent, or none (a header). */
export async function indexed(
  text: string,
  o: TabixOptions,
  extent: (cells: readonly string[]) => { seq: string; begin: number; end: number } | undefined,
): Promise<{ file: Uint8Array; index: Uint8Array }> {
  const bytes = new TextEncoder().encode(text);
  const file = await bgzip(bytes);
  const lines: IndexedLine[] = [];
  const seqOf = new Map(o.names.map((n, i) => [n, i]));
  // the text is ASCII: a character's offset is its byte's
  for (let start = 0; start < text.length; ) {
    const nl = text.indexOf('\n', start);
    const stop = nl < 0 ? text.length : nl + 1;
    const line = text.slice(start, nl < 0 ? text.length : nl);
    const e = line.startsWith(o.meta) ? undefined : extent(line.split('\t'));
    if (e) lines.push({ seq: seqOf.get(e.seq)!, begin: e.begin, end: e.end, start, stop });
    start = stop;
  }
  return { file: file.bytes, index: await tabix(lines, file, o) };
}
