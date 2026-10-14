/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST4b follow-up: a VCF's estimates use the compressed bytes a record takes in the file itself, measured from the
// first data blocks read, not ClinVar's fixed figure. The files here are built in the test (BGZF blocks and a .tbi):
// one chromosome of records every 100 positions, each carrying `width` hex digits of INFO, so a record's bytes differ
// between the two files by far more than the fixed figure's margin.

import { crc32, deflateRawSync } from "node:zlib";
import type { StreamRequest } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { type RangeReader, vcfProvider } from "../../../src/stream/gmod/index.js";
import { req } from "./cases.js";

const signal = new AbortController().signal;
const N = 3000;
const GAP = 100;
const HEADER = `##fileformat=VCFv4.2\n##contig=<ID=chr1,length=${N * GAP}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n`;

/** One BGZF block: gzip with the `BC` extra field (the block's size - 1). */
function block(data: Uint8Array): Uint8Array {
  const body = deflateRawSync(data);
  const out = new Uint8Array(18 + body.length + 8);
  const dv = new DataView(out.buffer);
  out.set([31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0]);
  dv.setUint16(16, out.length - 1, true);
  out.set(body, 18);
  dv.setUint32(18 + body.length, crc32(data), true);
  dv.setUint32(22 + body.length, data.length, true);
  return out;
}

/** A pseudo-random hex string (incompressible past its alphabet: about half a byte a digit). */
function hex(seed: number, n: number): string {
  let x = seed * 2654435761;
  let s = "";
  for (let i = 0; i < n; i++) {
    x = (Math.imul(x, 1103515245) + 12345) | 0;
    s += ((x >>> 16) & 15).toString(16);
  }
  return s;
}

/** The file (blocks: the header with `first` records, then 100 records each) and its index, one chromosome. */
function build(width: number, first: number) {
  const lines = Array.from({ length: N }, (_, i) => `chr1\t${i * GAP + 1}\t.\tA\tG\t.\t.\tX=${hex(i + 1, width)}\n`);
  const groups = [lines.slice(0, first)];
  for (let i = first; i < N; i += 100) groups.push(lines.slice(i, i + 100));
  const enc = new TextEncoder();
  const blocks = groups.map((g, k) => block(enc.encode((k === 0 ? HEADER : "") + g.join(""))));
  const at = [0];
  for (const b of blocks) at.push(at[at.length - 1] + b.length);
  const file = new Uint8Array(at[at.length - 1]);
  blocks.forEach((b, k) => {
    file.set(b, at[k]);
  });
  // The block of record i, for the linear index's first record of each 16 kb window.
  const blockOf = (i: number) => (i < first ? 0 : 1 + Math.floor((i - first) / 100));
  const windows = Math.ceil((N * GAP) / 16384);
  const linear = Array.from({ length: windows }, (_, k) => Math.ceil((k * 16384) / GAP));
  const bytes = new Uint8Array(36 + 5 + 4 + 8 + 16 + 4 + 8 * windows);
  const dv = new DataView(bytes.buffer);
  dv.setUint32(0, 0x01494254, true);
  dv.setInt32(4, 1, true);
  dv.setInt32(8, 2, true);
  [dv.setInt32(12, 1, true), dv.setInt32(16, 2, true), dv.setInt32(24, 35, true), dv.setInt32(32, 5, true)];
  bytes.set(enc.encode("chr1\0"), 36);
  let o = 41;
  const put = (v: number) => {
    dv.setUint32(o, v % 2 ** 32, true);
    dv.setUint32(o + 4, Math.floor(v / 2 ** 32), true);
    o += 8;
  };
  dv.setInt32(o, 1, true);
  dv.setUint32(o + 4, 0, true);
  dv.setInt32(o + 8, 1, true);
  o += 12;
  put(0);
  put(at[at.length - 1] * 65536);
  dv.setInt32(o, windows, true);
  o += 4;
  for (const i of linear) put(at[blockOf(i)] * 65536);
  return { file, index: block(bytes), size: (i: number) => blocks[i].length };
}

function reader(file: Uint8Array): RangeReader & { reads: [number, number][] } {
  const read = (async (offset: number, length: number) => {
    read.reads.push([offset, length]);
    return file.slice(offset, offset + length);
  }) as RangeReader & { reads: [number, number][] };
  read.reads = [];
  return read;
}

const at = (start: number, stop: number): StreamRequest => ({ ...req("chr1", start, stop), max_bytes: 1 });
/** The estimated count the provider gives for [start, stop): the sum of its bins (a data request gated by bytes). */
async function estimated(p: ReturnType<typeof vcfProvider>, start: number, stop: number): Promise<number> {
  const got = await p.fetch(at(start, stop), signal);
  if (got.kind !== "bins") throw new Error("not bins");
  return Array.from(got.bins.count as ArrayLike<number>).reduce((a, b) => a + b, 0);
}

describe("the VCF's bytes a record, measured from the file", () => {
  it("estimates a file of wide records and a file of narrow ones to within a fifth of their counts", async () => {
    for (const width of [8, 70]) {
      const f = build(width, 0);
      const p = vcfProvider(reader(f.file), f.index);
      await p.fetch(req(null), signal);
      // Records every 100 positions: [60,000, 160,000) holds 1,000.
      const sum = await estimated(p, 60_000, 160_000);
      expect(sum).toBeGreaterThan(800);
      expect(sum).toBeLessThan(1200);
    }
  });

  it("uses the fixed figure until a block past the header is read, then the file's own", async () => {
    // The header's block (with 300 records) nearly fills the first read; the blocks after it lie past it.
    const f = build(70, 300);
    const read = reader(f.file);
    const p = vcfProvider(read, f.index);
    await p.fetch(req(null), signal);
    const fixed = await estimated(p, 60_000, 160_000);
    // 28 bytes a record against about 60: the fixed figure counts about twice the records.
    expect(fixed).toBeGreaterThan(1400);
    // A region read (forced) measures the blocks it holds.
    const items = await p.fetch({ ...at(150_000, 200_000), force: true }, signal);
    expect(items.kind).toBe("items");
    const measured = await estimated(p, 60_000, 160_000);
    expect(measured).toBeGreaterThan(800);
    expect(measured).toBeLessThan(1200);
  });
});
