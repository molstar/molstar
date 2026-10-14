/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Genome files made by hand to break the readers (the review gate's untrusted input, GI-05 to GI-07): each is an error
// at once, never a hang or a memory blow-up.
import { crc32, deflateRawSync, deflateSync } from "node:zlib";
import { describe, expect, it } from "vitest";
import { bigBedProvider, sqlFields } from "../../../src/stream/gmod/bigbed.js";
import { vcfProvider } from "../../../src/stream/gmod/tabix.js";
import { twoBitProvider } from "../../../src/stream/gmod/twobit.js";
import { itWithin } from "../../model/timing.js";
import { req } from "./cases.js";

/** A little-endian byte writer. */
class W {
  readonly b: number[] = [];
  u8(v: number) {
    this.b.push(v & 255);
  }
  u16(v: number) {
    this.u8(v);
    this.u8(v >>> 8);
  }
  u32(v: number) {
    this.u16(v & 0xffff);
    this.u16(v >>> 16);
  }
  u64(v: number) {
    this.u32(v % 2 ** 32);
    this.u32(Math.floor(v / 2 ** 32));
  }
  bytes(a: readonly number[]) {
    for (const x of a) this.b.push(x);
  }
  get at() {
    return this.b.length;
  }
}

/** A BigBed of one chromosome "chr1" (or, with `leaves`, a root over one leaf per group of 4-letter names, chromosome
 * ids in order) and one block (zlib-compressed with `compressed`, its `uncompressBufSize`); a chromosome tree
 * `treeDepth` levels of nodes that point twice at the next; or an R-tree whose root points twice at itself. With
 * `autoSql`, the file names its fields (`fields` of them, 3 + the rest a record holds). */
function bigBed(o: {
  block: number[];
  treeDepth?: number;
  rtreeCycle?: boolean;
  leaves?: string[][];
  compressed?: number;
  autoSql?: string;
  fields?: number;
}): Uint8Array {
  const w = new W();
  const [KEY, CHR1] = [4, [99, 104, 114, 49]];
  for (let i = 0; i < 64; i++) w.u8(0);
  const sqlAt = o.autoSql ? w.at : 0;
  if (o.autoSql) w.bytes([...new TextEncoder().encode(o.autoSql), 0]);
  const TREE = w.at;
  const leaves = o.leaves ?? [];
  w.u32(0x78ca8c91);
  w.u32(2);
  w.u32(KEY);
  w.u32(8);
  w.u64(leaves.length ? leaves.flat().length : 1);
  w.u64(0);
  for (let d = 0; d < (o.treeDepth ?? 0); d++) {
    const next = w.at + 4 + 2 * (KEY + 8);
    w.u8(0);
    w.u8(0);
    w.u16(2);
    for (let k = 0; k < 2; k++) {
      w.bytes(CHR1);
      w.u64(next);
    }
  }
  if (leaves.length) {
    let next = w.at + 4 + leaves.length * (KEY + 8);
    w.u8(0);
    w.u8(0);
    w.u16(leaves.length);
    for (const names of leaves) {
      w.bytes([...new TextEncoder().encode(names[0])]);
      w.u64(next);
      next += 4 + names.length * (KEY + 8);
    }
    let id = 0;
    for (const names of leaves) {
      w.u8(1);
      w.u8(0);
      w.u16(names.length);
      for (const name of names) {
        w.bytes([...new TextEncoder().encode(name)]);
        w.u32(id++);
        w.u32(1_000_000);
      }
    }
  } else {
    w.u8(1);
    w.u8(0);
    w.u16(1);
    w.bytes(CHR1);
    w.u32(0);
    w.u32(1_000_000);
  }
  const dataAt = w.at;
  w.u64(1);
  const blockAt = w.at;
  w.bytes(o.compressed ? [...deflateSync(Uint8Array.from(o.block))] : o.block);
  const blockSize = w.at - blockAt;
  const indexAt = w.at;
  for (const v of [0x2468ace0, o.rtreeCycle ? 2 : 1]) w.u32(v);
  w.u64(1);
  for (const v of [0, 0, 0, 1_000_000]) w.u32(v);
  w.u64(indexAt);
  w.u32(1);
  w.u32(0);
  const root = w.at;
  w.u8(o.rtreeCycle ? 0 : 1);
  w.u8(0);
  w.u16(o.rtreeCycle ? 2 : 1);
  for (let k = 0; k < (o.rtreeCycle ? 2 : 1); k++) {
    for (const v of [0, 0, 0, 1_000_000]) w.u32(v);
    w.u64(o.rtreeCycle ? root : blockAt);
    if (!o.rtreeCycle) w.u64(blockSize);
  }
  for (let i = 0; i < 64; i++) w.u8(0); // slack for fixed-size node reads
  const b = new Uint8Array(w.b);
  const dv = new DataView(b.buffer);
  dv.setUint32(0, 0x8789f2eb, true);
  dv.setUint16(4, 4, true);
  dv.setBigUint64(8, BigInt(TREE), true);
  dv.setBigUint64(16, BigInt(dataAt), true);
  dv.setBigUint64(24, BigInt(indexAt), true);
  dv.setUint16(32, o.fields ?? 4, true);
  dv.setUint16(34, Math.min(o.fields ?? 4, 4), true);
  dv.setBigUint64(36, BigInt(sqlAt), true);
  dv.setUint32(52, o.compressed ?? 0, true);
  return b;
}
const record = (start: number, end: number, name: string, nul = true) => {
  const w = new W();
  w.u32(0);
  w.u32(start);
  w.u32(end);
  w.bytes([...new TextEncoder().encode(name)]);
  if (nul) w.u8(0);
  return w.b;
};
const reader = (b: Uint8Array) => async (o: number, n: number) => b.subarray(o, o + n);
const signal = new AbortController().signal;

describe("BigBed (GI-05, GI-06)", () => {
  it("reads a good one", async () => {
    const p = bigBedProvider(reader(bigBed({ block: [...record(10, 20, "a"), ...record(30, 40, "b")] })));
    const r = await p.fetch(req("chr1", 0, 1000), signal);
    expect(r.kind === "items" && Array.from(r.items.id as ArrayLike<unknown>)).toEqual(["a@10-20", "b@30-40"]);
  });

  it("a block whose last record has no NUL ends with it", async () => {
    const p = bigBedProvider(reader(bigBed({ block: [...record(10, 20, "a"), ...record(30, 40, "b", false)] })));
    const r = await p.fetch(req("chr1", 0, 1000), signal);
    expect(r.kind === "items" && Array.from(r.items.id as ArrayLike<unknown>)).toEqual(["a@10-20", "b@30-40"]);
  });

  itWithin("a chromosome tree of shared children is walked once per node, and no deeper than 32", async (within) => {
    const ok = await within(1000, () =>
      bigBedProvider(reader(bigBed({ block: record(10, 20, "a"), treeDepth: 30 }))).fetch(req(null), signal),
    );
    expect(ok.kind === "manifest" && ok.manifest.systems.map((s) => s.key)).toEqual(["chr1"]);
    const deep = bigBedProvider(reader(bigBed({ block: record(10, 20, "a"), treeDepth: 40 })));
    await expect(deep.fetch(req(null), signal)).rejects.toThrow(/deeper than 32/);
  });

  it("lists the chromosomes in the tree's order across leaves (LA-03: the walk reversed each node's children)", async () => {
    const leaves = [["chrA", "chrB"], ["chrC", "chrD"], ["chrE"]];
    const p = bigBedProvider(reader(bigBed({ block: record(10, 20, "a"), leaves })));
    const m = await p.fetch(req(null), signal);
    expect(m.kind === "manifest" && m.manifest.systems.map((s) => s.key)).toEqual([
      "chrA",
      "chrB",
      "chrC",
      "chrD",
      "chrE",
    ]);
  });

  it("a compressed block that inflates past the header's uncompressBufSize is refused (GI-07)", async () => {
    const block = Array.from({ length: 40 }, (_, k) => record(10 * k, 10 * k + 5, `item${k}`)).flat();
    const fits = bigBedProvider(reader(bigBed({ block, compressed: block.length })));
    const r = await fits.fetch(req("chr1", 0, 1000), signal);
    expect(r.kind === "items" && (r.items.id as ArrayLike<unknown>).length).toBe(40);
    const past = bigBedProvider(reader(bigBed({ block, compressed: block.length - 1 })));
    await expect(past.fetch(req("chr1", 0, 1000), signal)).rejects.toThrow(/inflates past \d+ bytes/);
  });

  it("an autoSql field named __proto__, or like a column the reader makes, is a column of its own (LA-39)", async () => {
    const autoSql = `table t
"t"
(
string chrom; "c"
uint chromStart; "s"
uint chromEnd; "e"
string name; "n"
string start; "a field named like the item's start"
string __proto__; "a field named like the prototype"
)`;
    const p = bigBedProvider(reader(bigBed({ block: record(10, 20, "a\t7\tx"), autoSql, fields: 6 })));
    const r = await p.fetch(req("chr1", 0, 1000), signal);
    if (r.kind !== "items") throw new Error(r.kind);
    const t = r.items as unknown as Record<string, ArrayLike<unknown>>;
    const PROTO = "__proto__"; // read as data: a literal `.__proto__` reads the prototype
    expect([Object.keys(t), Array.from(t.start), Array.from(t.field5), Array.from(t[PROTO])]).toEqual([
      ["start", "stop", "id", "label", "name", "field5", "__proto__"],
      [10],
      ["7"],
      ["x"],
    ]);
  });

  it("an R-tree whose root points at itself visits it once", async () => {
    const p = bigBedProvider(reader(bigBed({ block: record(10, 20, "a"), rtreeCycle: true })));
    const r = await p.fetch(req("chr1", 0, 1000), signal);
    expect(r.kind === "items" && Array.from(r.items.id as ArrayLike<unknown>)).toEqual([]);
  });
});

describe("BGZF (GI-07)", () => {
  /** One BGZF member holding `raw`, with its ISIZE `isize`. */
  const member = (raw: Uint8Array, isize = raw.length) => {
    const body = deflateRawSync(raw, { level: 9 });
    const size = 18 + body.length + 8;
    const block = new Uint8Array(size);
    const dv = new DataView(block.buffer);
    block.set([0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff]);
    dv.setUint16(10, 6, true);
    block[12] = 66;
    block[13] = 67;
    dv.setUint16(14, 2, true);
    dv.setUint16(16, size - 1, true);
    block.set(body, 18);
    dv.setUint32(18 + body.length, crc32(raw), true);
    dv.setUint32(22 + body.length, isize, true);
    return block;
  };
  const open = (index: Uint8Array) =>
    vcfProvider(async () => new Uint8Array(0), index).fetch(req(null), new AbortController().signal);

  it("a block past the format's 64 KB is refused before it is held", async () => {
    const before = process.memoryUsage().rss;
    await expect(open(member(new Uint8Array(8 << 20)))).rejects.toThrow(/inflates past 65536 bytes/);
    expect(process.memoryUsage().rss - before).toBeLessThan(64 << 20);
  });

  it("an index fetched by a function gets the open's signal, aborted when its request leaves (LA-36)", async () => {
    let given: AbortSignal | undefined;
    const index = (s?: AbortSignal) => {
      given = s;
      return new Promise<Uint8Array>(() => {});
    };
    const a = new AbortController();
    const opening = vcfProvider(async () => new Uint8Array(0), index).fetch(req(null), a.signal);
    await new Promise((t) => setTimeout(t, 1));
    a.abort();
    await expect(opening).rejects.toThrow(/aborted/);
    expect(given?.aborted).toBe(true);
  });

  // The gzip stream checks the trailer itself and refuses data after a member (two members under one BC size), so in
  // Node no block reaches our own ISIZE check, a second line (LA-34).
  it("a block that is not its ISIZE is refused (by the gzip stream's trailer check, or ours)", async () => {
    await expect(open(member(new Uint8Array(100), 50))).rejects.toThrow();
  });
});

describe("lengths and counts from a header (GI-11, GI-12)", () => {
  it("a BigBed whose header puts its data a terabyte past the tree is refused before the read", async () => {
    const b = bigBed({ block: record(10, 20, "a") });
    new DataView(b.buffer).setBigUint64(16, 2n ** 40n, true);
    const asked: number[] = [];
    const p = bigBedProvider(async (o, n) => {
      asked.push(n);
      return b.subarray(o, o + n);
    });
    await expect(p.fetch(req(null), signal)).rejects.toThrow(/chromosome tree is \d+ bytes: past its bound/);
    expect(Math.max(...asked)).toBeLessThan(1 << 20);
  });

  /** A 2bit index of `n` one-letter sequences `s0`, `s1`, …, each record a size word of 1. */
  const twoBit = (n: number) => {
    const w = new W();
    for (const v of [0x1a412743, 0, n, 0]) w.u32(v);
    const names = Array.from({ length: n }, (_, i) => `s${i}`);
    const indexEnd = 16 + names.reduce((t, x) => t + 5 + x.length, 0);
    names.forEach((name, i) => {
      w.u8(name.length);
      w.bytes([...new TextEncoder().encode(name)]);
      w.u32(indexEnd + 16 * i);
    });
    for (let i = 0; i < n; i++) for (const v of [1, 0, 0, 0]) w.u32(v);
    return new Uint8Array(w.b);
  };

  itWithin(
    "a 2bit's sizes are read at most 16 at a time, and a count past its index is refused where the index ends",
    async (within) => {
      const b = twoBit(2000);
      let inFlight = 0;
      let most = 0;
      const p = twoBitProvider(async (o, n) => {
        inFlight++;
        most = Math.max(most, inFlight);
        await new Promise((ok) => setTimeout(ok, 0));
        inFlight--;
        return b.subarray(o, o + n);
      });
      const m = await p.fetch(req(null), signal);
      expect(m.kind === "manifest" && m.manifest.systems.length).toBe(2000);
      expect(most).toBeLessThanOrEqual(16);
      const huge = twoBit(1);
      new DataView(huge.buffer).setUint32(8, 4_000_000_000, true);
      await within(1000, () =>
        expect(twoBitProvider(reader(huge)).fetch(req(null), signal)).rejects.toThrow(
          /the index ends at sequence \d+ of the 4000000000 its header counts/,
        ),
      );
    },
  );

  it("a 2bit's bound counts the sequences `only` keeps, not the index's (LA-33: a draft assembly of 10M scaffolds)", async () => {
    const b = twoBit(2000);
    const few = twoBitProvider(reader(b), { only: (n) => n === "s7" || n === "s1999", maxSequences: 1000 });
    const m = await few.fetch(req(null), signal);
    expect(m.kind === "manifest" && m.manifest.systems.map((x) => x.key)).toEqual(["s7", "s1999"]);
    await expect(twoBitProvider(reader(b), { maxSequences: 1000 }).fetch(req(null), signal)).rejects.toThrow(
      /over 1000 sequences kept, the bound: name the ones to read/,
    );
  });

  it("a BigBed's autoSql and R-tree node are bounded (GI-11, GI-05)", async () => {
    const sql = bigBed({ block: record(10, 20, "a"), autoSql: "table t" });
    new DataView(sql.buffer).setBigUint64(44, 2n ** 30n, true);
    await expect(bigBedProvider(reader(sql)).fetch(req(null), signal)).rejects.toThrow(
      /the header's autoSql is \d+ bytes: past its bound of 1048576/,
    );
    const wide = bigBed({ block: record(10, 20, "a") });
    const dv = new DataView(wide.buffer);
    dv.setUint32(Number(dv.getBigUint64(24, true)) + 4, 70_000, true);
    await expect(bigBedProvider(reader(wide)).fetch(req("chr1", 0, 1000), signal)).rejects.toThrow(
      /block size 70000: past 65536 items/,
    );
  });
});

describe("BigBed autoSql (CB-03)", () => {
  it("names the fields after an enum or a set, its values over lines, in order", () => {
    const sql = `table bigGenePred
"bigGenePred gene models"
   (
   string chrom;       "Reference sequence chromosome or scaffold"
   uint   chromStart;  "Start position in chromosome"
   uint   chromEnd;    "End position in chromosome"
   string name;        "Name or ID of item"
   uint score;         "Score (0-1000)"
   enum(cmpl, incmpl,
        unk, none) cdsStartStat; "Status of CDS start annotation"
   set(a, b) flags;    "Flags"
   int[blockCount] blockSizes; "Comma separated list of block sizes"
   string geneName;    "Primary identifier for gene"
   )`;
    expect(sqlFields(sql)).toEqual(["name", "score", "cdsStartStat", "flags", "blockSizes", "geneName"]);
  });
});
