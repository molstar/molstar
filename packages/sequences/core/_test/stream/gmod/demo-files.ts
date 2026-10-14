/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A small synthetic genome written in the four genome file formats the readers take, whole, so any region of it reads
// (the recorded slices of hg38's files answer only the reads they recorded): two sequences, `chrS1` (24,000 bases) and
// `chrS2` (6,000), as a 2bit (an N run, soft-masked repeats), a BigWig of their GC percent in 50-base windows, a BigBed
// of BED12 genes and a bgzipped VCF with its tabix index. Synthetic, deterministic, made in memory (no fixture file):
// the File formats story and demo-files.test.ts read them. Pure TypeScript (no `node:` imports, the page imports it):
// BGZF blocks hold stored (uncompressed) deflate blocks, BigWig and BigBed blocks are not compressed
// (`uncompressBufSize` 0), all little-endian.

import type { RangeReader } from "../../../src/stream/gmod/index.js";

/** The sequences and their lengths. */
export const DEMO_SYSTEMS: readonly [string, number][] = [
  ["chrS1", 24_000],
  ["chrS2", 6_000],
];
/** chrS1's leading N run and its soft-masked repeats, [start, stop). */
const N_RUN: [number, number] = [0, 400];
const REPEATS: [number, number][] = [
  [3_000, 3_300],
  [9_000, 9_600],
  [17_000, 17_250],
];
/** The GC windows (BigWig), and items a BigWig block holds. */
const GC_WINDOW = 50;
const GC_BLOCK = 120;

/** BED12 genes: chrom, name, strand, exons [start, stop) (the thick part from the first exon's 100th base). */
const GENES: [string, string, "+" | "-", [number, number][]][] = [
  [
    "chrS1",
    "GENA",
    "+",
    [
      [2_000, 2_300],
      [3_800, 4_100],
      [5_200, 5_600],
      [7_000, 7_500],
    ],
  ],
  [
    "chrS1",
    "GENB",
    "-",
    [
      [9_800, 10_400],
      [11_900, 12_200],
      [13_700, 14_200],
    ],
  ],
  ["chrS1", "GENC", "+", [[15_500, 16_800]]],
  [
    "chrS1",
    "GENC",
    "+",
    [
      [15_500, 15_900],
      [18_200, 19_000],
    ],
  ],
  [
    "chrS1",
    "GEND",
    "-",
    [
      [20_000, 20_500],
      [21_100, 21_300],
      [22_000, 22_400],
      [23_000, 23_500],
    ],
  ],
  [
    "chrS2",
    "GENE",
    "+",
    [
      [1_000, 1_500],
      [3_200, 4_000],
    ],
  ],
];
const CLNSIG = ["Benign", "Likely_benign", "Uncertain_significance", "Likely_pathogenic", "Pathogenic"];
/** VCF records a BGZF data block holds. */
const VCF_BLOCK = 25;

/** A deterministic generator (a linear congruential one), 0..1. */
function random(seed: number) {
  let x = seed >>> 0;
  return () => {
    x = (Math.imul(x, 1664525) + 1013904223) >>> 0;
    return x / 2 ** 32;
  };
}

/** A sequence's letters: GC content rising and falling along it, N over the N run, repeats in lower case. */
function letters(name: string, length: number): string {
  const r = random(name === "chrS1" ? 17 : 29);
  let out = "";
  for (let i = 0; i < length; i++) {
    const gc = 0.42 + 0.18 * Math.sin(i / 1_900) + 0.08 * Math.sin(i / 230);
    const u = r();
    const base = u < gc ? (r() < 0.5 ? "G" : "C") : r() < 0.5 ? "A" : "T";
    out += base;
  }
  if (name !== "chrS1") return out;
  out = "N".repeat(N_RUN[1] - N_RUN[0]) + out.slice(N_RUN[1]);
  for (const [a, b] of REPEATS) out = out.slice(0, a) + out.slice(a, b).toLowerCase() + out.slice(b);
  return out;
}

/** A little-endian writer. */
class Writer {
  bytes = new Uint8Array(1024);
  at = 0;
  private room(n: number) {
    if (this.at + n <= this.bytes.length) return;
    const next = new Uint8Array(Math.max(this.bytes.length * 2, this.at + n));
    next.set(this.bytes);
    this.bytes = next;
  }
  private dv() {
    return new DataView(this.bytes.buffer);
  }
  u8(v: number) {
    this.room(1);
    this.bytes[this.at++] = v & 255;
  }
  u16(v: number) {
    this.room(2);
    this.dv().setUint16(this.at, v, true);
    this.at += 2;
  }
  u32(v: number) {
    this.room(4);
    this.dv().setUint32(this.at, v >>> 0, true);
    this.at += 4;
  }
  u64(v: number) {
    this.u32(v % 2 ** 32);
    this.u32(Math.floor(v / 2 ** 32));
  }
  f32(v: number) {
    this.room(4);
    this.dv().setFloat32(this.at, v, true);
    this.at += 4;
  }
  f64(v: number) {
    this.room(8);
    this.dv().setFloat64(this.at, v, true);
    this.at += 8;
  }
  raw(b: Uint8Array) {
    this.room(b.length);
    this.bytes.set(b, this.at);
    this.at += b.length;
  }
  /** Writes `v` as a u32 or u64 at `at`, later (an offset known once what follows is written). */
  put32(at: number, v: number) {
    this.dv().setUint32(at, v >>> 0, true);
  }
  put64(at: number, v: number) {
    this.put32(at, v % 2 ** 32);
    this.put32(at + 4, Math.floor(v / 2 ** 32));
  }
  done() {
    return this.bytes.slice(0, this.at);
  }
}
const ascii = (s: string) => Uint8Array.from(s, (c) => c.charCodeAt(0));

// ---- 2bit ----

const CODE: Record<string, number> = { T: 0, C: 1, A: 2, G: 3 };
/** Runs of `test` over `s`, [start, size]. */
function runs(s: string, test: (c: string) => boolean): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i < s.length; ) {
    if (!test(s[i])) {
      i++;
      continue;
    }
    let e = i;
    while (e < s.length && test(s[e])) e++;
    out.push([i, e - i]);
    i = e;
  }
  return out;
}
function twoBit(seqs: [string, string][]): Uint8Array {
  const w = new Writer();
  w.u32(0x1a412743);
  w.u32(0);
  w.u32(seqs.length);
  w.u32(0);
  const offsets: number[] = [];
  for (const [name] of seqs) {
    w.u8(name.length);
    w.raw(ascii(name));
    offsets.push(w.at);
    w.u32(0);
  }
  for (const [k, [, s]] of seqs.entries()) {
    w.put32(offsets[k], w.at);
    w.u32(s.length);
    const ns = runs(s, (c) => c === "N" || c === "n");
    const mask = runs(s, (c) => c !== c.toUpperCase());
    for (const list of [ns, mask]) {
      w.u32(list.length);
      for (const [a] of list) w.u32(a);
      for (const [, n] of list) w.u32(n);
    }
    w.u32(0);
    for (let i = 0; i < s.length; i += 4) {
      let byte = 0;
      for (let j = 0; j < 4; j++) byte = (byte << 2) | (CODE[(s[i + j] ?? "T").toUpperCase()] ?? 0);
      w.u8(byte);
    }
  }
  return w.done();
}

// ---- BigWig and BigBed (one container) ----

interface Block {
  chrom: number;
  start: number;
  stop: number;
  bytes: Uint8Array;
}
/**
 * A BigWig or BigBed: the header, the total summary, the chromosome B+ tree (one leaf), the data's item count and its
 * blocks, and an R-tree of one leaf node over the blocks. No zoom levels, no autoSql, blocks not compressed.
 */
function bbi(kind: "bigwig" | "bigbed", chroms: readonly [string, number][], blocks: Block[], o: Bbi): Uint8Array {
  const w = new Writer();
  w.u32(kind === "bigwig" ? 0x888ffc26 : 0x8789f2eb);
  w.u16(4);
  w.u16(0); // zoom levels
  const [treeAt, dataAt, indexAt, summaryAt] = [8, 16, 24, 44];
  for (let i = 8; i < 64; i++) w.u8(0);
  w.bytes[32] = o.fieldCount;
  w.bytes[34] = o.definedFieldCount;
  // The total summary: bases covered, min, max, sum, sum of squares.
  w.put64(summaryAt, w.at);
  w.u64(o.covered);
  w.f64(o.min);
  w.f64(o.max);
  w.f64(o.sum);
  w.f64(o.sumSquares);
  // The chromosome tree: its header, then one leaf of (name NUL-padded to the longest, id, size).
  const keySize = Math.max(...chroms.map(([n]) => n.length));
  w.put64(treeAt, w.at);
  w.u32(0x78ca8c91);
  w.u32(chroms.length);
  w.u32(keySize);
  w.u32(8);
  w.u64(chroms.length);
  w.u64(0);
  w.u8(1);
  w.u8(0);
  w.u16(chroms.length);
  for (const [id, [name, size]] of chroms.entries()) {
    w.raw(ascii(name.padEnd(keySize, "\0")));
    w.u32(id);
    w.u32(size);
  }
  // The data: its count (BigBed: items; BigWig: sections), then the blocks.
  w.put64(dataAt, w.at);
  w.u64(o.count);
  const placed = blocks.map((b) => {
    const at = w.at;
    w.raw(b.bytes);
    return { ...b, at };
  });
  // The R-tree: its 48-byte header, then one leaf node over every block.
  w.put64(indexAt, w.at);
  w.u32(0x2468ace0);
  w.u32(256); // block size: items a node may hold
  w.u64(placed.length);
  const [first, last] = [placed[0], placed[placed.length - 1]];
  for (const v of [first.chrom, first.start, last.chrom, last.stop]) w.u32(v);
  w.u64(w.at + 8 + 4 + 4 + 4 + placed.length * 32); // the end of the file: the rest of the header, the leaf
  w.u32(1);
  w.u32(0);
  w.u8(1);
  w.u8(0);
  w.u16(placed.length);
  for (const b of placed) {
    for (const v of [b.chrom, b.start, b.chrom, b.stop]) w.u32(v);
    w.u64(b.at);
    w.u64(b.bytes.length);
  }
  return w.done();
}
interface Bbi {
  fieldCount: number;
  definedFieldCount: number;
  count: number;
  covered: number;
  min: number;
  max: number;
  sum: number;
  sumSquares: number;
}

/** GC percent per window of each sequence (N and nothing else left out), as fixedStep sections of `GC_BLOCK` items. */
function gcWig(seqs: [string, string][]) {
  const blocks: Block[] = [];
  const lines: string[] = [];
  let [covered, min, max, sum, sumSquares, count] = [0, Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, 0, 0, 0];
  for (const [chrom, [name, s]] of seqs.entries()) {
    const values: number[] = [];
    for (let i = 0; i < s.length; i += GC_WINDOW) {
      const window = s.slice(i, i + GC_WINDOW).toUpperCase();
      const gc = [...window].filter((c) => c === "G" || c === "C").length;
      values.push(Math.round((100 * gc) / window.length));
    }
    lines.push(`fixedStep chrom=${name} start=1 step=${GC_WINDOW} span=${GC_WINDOW}`);
    for (const v of values) lines.push(String(v));
    for (let k = 0; k < values.length; k += GC_BLOCK) {
      const part = values.slice(k, k + GC_BLOCK);
      const start = k * GC_WINDOW;
      const stop = Math.min(s.length, (k + part.length) * GC_WINDOW);
      const w = new Writer();
      for (const v of [chrom, start, stop, GC_WINDOW, GC_WINDOW]) w.u32(v);
      w.u8(3); // fixedStep
      w.u8(0);
      w.u16(part.length);
      for (const v of part) w.f32(v);
      blocks.push({ chrom, start, stop, bytes: w.done() });
      count++;
    }
    for (const v of values) {
      covered += GC_WINDOW;
      min = Math.min(min, v);
      max = Math.max(max, v);
      sum += v * GC_WINDOW;
      sumSquares += v * v * GC_WINDOW;
    }
  }
  return { blocks, text: lines.join("\n"), stats: { count, covered, min, max, sum, sumSquares } };
}

/** The genes as BED12 lines and BigBed blocks (`perBlock` items each, in order). */
function genesBed(names: readonly string[], perBlock = 2) {
  const rows = GENES.map(([chrom, name, strand, exons]) => {
    const [start, stop] = [exons[0][0], exons[exons.length - 1][1]];
    const rest = [
      name,
      "0",
      strand,
      String(exons[0][0] + 100),
      String(stop - 100),
      "0",
      String(exons.length),
      `${exons.map(([a, b]) => b - a).join(",")},`,
      `${exons.map(([a]) => a - start).join(",")},`,
    ];
    return { chrom: names.indexOf(chrom), name: chrom, start, stop, rest: rest.join("\t") };
  });
  const blocks: Block[] = [];
  for (let k = 0; k < rows.length; k += perBlock) {
    const part = rows.slice(k, k + perBlock).filter((r) => r.chrom === rows[k].chrom);
    const w = new Writer();
    for (const r of part) {
      w.u32(r.chrom);
      w.u32(r.start);
      w.u32(r.stop);
      w.raw(ascii(`${r.rest}\0`));
    }
    blocks.push({
      chrom: part[0].chrom,
      start: part[0].start,
      stop: Math.max(...part.map((r) => r.stop)),
      bytes: w.done(),
    });
    // A block holds one chromosome's items: the next starts after the last taken.
    k -= perBlock - part.length;
  }
  const text = rows.map((r) => `${r.name}\t${r.start}\t${r.stop}\t${r.rest}`).join("\n");
  return { blocks, text, count: rows.length };
}

// ---- VCF, BGZF and tabix ----

const CRC = (() => {
  const t = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c >>> 0;
  }
  return t;
})();
function crc32(b: Uint8Array): number {
  let c = 0xffffffff;
  for (let i = 0; i < b.length; i++) c = CRC[(c ^ b[i]) & 255] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}
/** One BGZF block: a gzip member with the `BC` field, its data in one stored deflate block (at most 65,280 bytes). */
function bgzfBlock(data: Uint8Array): Uint8Array {
  if (data.length > 65_280) throw new Error("a BGZF block holds at most 65,280 bytes here");
  const w = new Writer();
  for (const v of [31, 139, 8, 4, 0, 0, 0, 0, 0, 255]) w.u8(v);
  w.u16(6);
  w.u8(66);
  w.u8(67);
  w.u16(2);
  const sizeAt = w.at;
  w.u16(0);
  w.u8(1); // the final block, stored
  w.u16(data.length);
  w.u16(~data.length & 0xffff);
  w.raw(data);
  w.u32(crc32(data));
  w.u32(data.length);
  const out = w.done();
  new DataView(out.buffer).setUint16(sizeAt, out.length - 1, true);
  return out;
}
/** The tabix bin of [beg, end) (the SAM specification's `reg2bin`). */
function reg2bin(beg: number, end: number): number {
  const e = end - 1;
  if (beg >> 14 === e >> 14) return ((1 << 15) - 1) / 7 + (beg >> 14);
  if (beg >> 17 === e >> 17) return ((1 << 12) - 1) / 7 + (beg >> 17);
  if (beg >> 20 === e >> 20) return ((1 << 9) - 1) / 7 + (beg >> 20);
  if (beg >> 23 === e >> 23) return ((1 << 6) - 1) / 7 + (beg >> 23);
  if (beg >> 26 === e >> 26) return ((1 << 3) - 1) / 7 + (beg >> 26);
  return 0;
}

interface Variant {
  chrom: string;
  pos: number; // 1-based
  id: string;
  ref: string;
  alt: string;
  clnsig: string;
}
/** About one variant every 400 bases of each sequence past chrS1's N run: SNVs, and every seventh an insertion or deletion. */
function variants(seqs: [string, string][]): Variant[] {
  const out: Variant[] = [];
  const r = random(5);
  let n = 1;
  for (const [name, s] of seqs)
    for (let p = 450 + Math.floor(r() * 200); p < s.length - 10; p += 250 + Math.floor(r() * 300)) {
      const ref1 = s[p].toUpperCase();
      const kind = n % 7 === 0 ? "ins" : n % 7 === 3 ? "del" : "snv";
      const ref = kind === "del" ? s.slice(p, p + 3).toUpperCase() : ref1;
      const alt =
        kind === "snv" ? ("ACGT".replace(ref1, "")[Math.floor(r() * 3)] ?? "A") : kind === "ins" ? `${ref1}TG` : ref1;
      out.push({ chrom: name, pos: p + 1, id: `rsS${n++}`, ref, alt, clnsig: CLNSIG[Math.floor(r() * CLNSIG.length)] });
    }
  return out;
}
/** The VCF text, its BGZF file (the header a block, then `VCF_BLOCK` records a block, then the empty last block) and its index. */
function vcfFiles(seqs: [string, string][]) {
  const header = [
    "##fileformat=VCFv4.2",
    "##source=Seq* demo-files.ts (synthetic)",
    ...seqs.map(([name, s]) => `##contig=<ID=${name},length=${s.length}>`),
    '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance (synthetic)">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
  ].join("\n");
  const list = variants(seqs);
  const line = (v: Variant) => `${v.chrom}\t${v.pos}\t${v.id}\t${v.ref}\t${v.alt}\t.\tPASS\tCLNSIG=${v.clnsig}\n`;
  const blocks: Uint8Array[] = [bgzfBlock(ascii(`${header}\n`))];
  /** Each record's virtual offsets [begin, end): its block's file offset << 16 | its offset in the block's data. */
  const voffs: [number, number][] = [];
  let fileAt = blocks[0].length;
  for (let k = 0; k < list.length; k += VCF_BLOCK) {
    let u = 0;
    const lines = list.slice(k, k + VCF_BLOCK).map(line);
    for (const l of lines) {
      voffs.push([fileAt * 65536 + u, fileAt * 65536 + u + l.length]);
      u += l.length;
    }
    const b = bgzfBlock(ascii(lines.join("")));
    blocks.push(b);
    fileAt += b.length;
  }
  blocks.push(bgzfBlock(new Uint8Array(0)));
  const file = new Uint8Array(blocks.reduce((n, b) => n + b.length, 0));
  let at = 0;
  for (const b of blocks) {
    file.set(b, at);
    at += b.length;
  }
  // The index: per sequence its bins' chunks (a run of records in one bin is one chunk) and its linear index.
  const w = new Writer();
  w.raw(ascii("TBI\x01"));
  w.u32(seqs.length);
  w.u32(2); // VCF
  w.u32(1); // the sequence column
  w.u32(2); // the begin column
  w.u32(0); // the end column
  w.u32(35); // '#'
  w.u32(0); // lines skipped
  const names = seqs.map(([n]) => `${n}\0`).join("");
  w.u32(names.length);
  w.raw(ascii(names));
  for (const [name, s] of seqs) {
    const bins = new Map<number, [number, number][]>();
    const linear = new Array<number>(Math.ceil(s.length / 16384)).fill(-1);
    list.forEach((v, i) => {
      if (v.chrom !== name) return;
      const [beg, end] = [v.pos - 1, v.pos - 1 + v.ref.length];
      const chunks = bins.get(reg2bin(beg, end)) ?? [];
      const last = chunks[chunks.length - 1];
      if (last && last[1] === voffs[i][0]) last[1] = voffs[i][1];
      else chunks.push([...voffs[i]]);
      bins.set(reg2bin(beg, end), chunks);
      for (let k = beg >> 14; k <= (end - 1) >> 14; k++) if (linear[k] < 0) linear[k] = voffs[i][0];
    });
    w.u32(bins.size);
    for (const [bin, chunks] of bins) {
      w.u32(bin);
      w.u32(chunks.length);
      for (const [a, b] of chunks) {
        w.u64(a);
        w.u64(b);
      }
    }
    // A window with no record starts where the next one with a record does.
    for (let k = linear.length - 1, next = voffs[voffs.length - 1][1]; k >= 0; k--) {
      if (linear[k] < 0) linear[k] = next;
      next = linear[k];
    }
    w.u32(linear.length);
    for (const v of linear) w.u64(v);
  }
  const tbi = new Uint8Array([...bgzfBlock(w.done()), ...bgzfBlock(new Uint8Array(0))]);
  const text = `${header}\n${list.map(line).join("")}`;
  return { file, tbi, text, count: list.length };
}

/** The demo genome's files and the text they were written from. */
export interface DemoGenome {
  readonly twobit: Uint8Array;
  readonly bigwig: Uint8Array;
  readonly bigbed: Uint8Array;
  readonly vcf: Uint8Array;
  readonly tbi: Uint8Array;
  /** What each holds, as text: FASTA, fixedStep wiggle, BED12 and the VCF itself. */
  readonly text: { readonly fasta: string; readonly wig: string; readonly bed: string; readonly vcf: string };
}

let made: DemoGenome | null = null;
/** The files, made once. */
export function demoGenome(): DemoGenome {
  if (made) return made;
  const seqs = DEMO_SYSTEMS.map(([name, n]): [string, string] => [name, letters(name, n)]);
  const names = seqs.map(([n]) => n);
  const gc = gcWig(seqs);
  const genes = genesBed(names);
  const vcf = vcfFiles(seqs);
  made = {
    twobit: twoBit(seqs),
    bigwig: bbi("bigwig", DEMO_SYSTEMS, gc.blocks, { fieldCount: 0, definedFieldCount: 0, ...gc.stats }),
    bigbed: bbi("bigbed", DEMO_SYSTEMS, genes.blocks, {
      fieldCount: 12,
      definedFieldCount: 12,
      count: genes.count,
      covered: 0,
      min: 0,
      max: 0,
      sum: 0,
      sumSquares: 0,
    }),
    vcf: vcf.file,
    tbi: vcf.tbi,
    text: {
      fasta: seqs.map(([n, s]) => `>${n}\n${s.replace(/(.{60})/g, "$1\n")}`).join("\n"),
      wig: gc.text,
      bed: genes.text,
      vcf: vcf.text,
    },
  };
  return made;
}

/** A reader over bytes in memory (`length` bytes from `offset`, fewer at the end), each read logged in `reads`. */
export function memoryReader(bytes: Uint8Array): RangeReader & { reads: [number, number][] } {
  const read = (async (offset: number, length: number) => {
    read.reads.push([offset, length]);
    return bytes.slice(offset, offset + length);
  }) as RangeReader & { reads: [number, number][] };
  read.reads = [];
  return read;
}
