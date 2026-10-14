/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { decodeMsgPack, encodeMsgPack } from "@molstar/seq-core/model";
import { checkQuery } from "@molstar/seq-core/query";
// Untrusted input (the review gate's gap, GI-*): a hostile or broken file of a few bytes to a few hundred KB yields
// an issue, never a hang, a crash or a memory blow-up.
import { describe, expect, it } from "vitest";
import { readDocument, readSvsx, toJSON, toMsgPack, validate, zipFiles } from "../../src/format/index.js";
import { KernelFormats } from "../../src/format/load/formats.js";
import { itWithin } from "../model/timing.js";
import { CASES } from "./cases.js";
import { kernel } from "./fixtures.js";

const tsv = KernelFormats.find((f) => f.name === "tsv");
const read = async (text: string) => {
  const said: string[] = [];
  const ctx = { report: (i: { code: string; message: string }) => said.push(`${i.code} ${i.message}`) } as never;
  const r = await tsv?.read(new TextEncoder().encode(text), {} as never, ctx);
  return { r, said };
};

describe("TSV and CSV (GI-04)", () => {
  itWithin("a wide header over many empty lines allocates by the text's length, not header × lines", async (within) => {
    const header = Array.from({ length: 20_000 }, (_, k) => `c${k}`).join("\t");
    const text = `${header}${"\n".repeat(200_000)}`;
    const before = process.memoryUsage().heapUsed;
    const { r, said } = await within(2000, () => read(text));
    expect(process.memoryUsage().heapUsed - before).toBeLessThan(200 << 20);
    // Blank lines are no rows: a table of 20,000 empty columns.
    expect([said, (r as { length?: number }).length]).toEqual([[], 0]);
  });

  itWithin("finds a repeated name in a wide header in linear time", async (within) => {
    const names = Array.from({ length: 50_000 }, (_, k) => `c${k}`);
    const { said } = await within(1000, () => read(`${[...names, "c7"].join("\t")}\n`));
    expect(said).toEqual([expect.stringMatching(/the header names "c7" twice/)]);
  });
});

/** The offsets of a zip's structures: the end record, each central entry and each local header. */
function layout(z: Uint8Array) {
  const v = new DataView(z.buffer, z.byteOffset, z.byteLength);
  let eocd = z.length - 22;
  while (v.getUint32(eocd, true) !== 0x06054b50) eocd--;
  const central: number[] = [];
  for (let k = 0, o = v.getUint32(eocd + 16, true); k < v.getUint16(eocd + 10, true); k++) {
    central.push(o);
    o += 46 + v.getUint16(o + 28, true) + v.getUint16(o + 30, true) + v.getUint16(o + 32, true);
  }
  const local = central.map((c) => v.getUint32(c + 42, true));
  return { v, eocd, central, local };
}
const enc = (t: string) => new TextEncoder().encode(t);
const archive = (extra: Record<string, Uint8Array> = {}) =>
  zipFiles({ "index.svsj": enc(toJSON(CASES.case1())), ...extra });
const said = async (z: Uint8Array, opts = {}) => {
  const r = await readSvsx(z, opts);
  return r.issues.map((i) => `${i.code} ${i.message}`);
};

describe("archives (GI-01, GI-02, GI-08, GI-09)", () => {
  it("reads a good one", async () => {
    expect(await said(await archive())).toEqual([]);
  });

  itWithin("a truncated archive, with no end record, is invalid at once", async (within) => {
    const z = await archive();
    expect(await within(1000, () => said(z.slice(0, z.length - 22)))).toEqual([
      expect.stringMatching(/^svs\.format\.invalid .*no end of central directory/),
    ]);
  });

  it("refuses past its bounds before inflating: entries, and their declared sizes", async () => {
    const big = new Uint8Array(4 << 20); // 4 MB of zeros: a few KB deflated
    const z = await archive({ "big.bin": big });
    expect(z.length).toBeLessThan(64 << 10);
    expect(await said(z, { maxBytes: 1 << 20 })).toEqual([expect.stringMatching(/over the bound of 1048576 bytes/)]);
    expect(await said(z, { maxEntries: 1 })).toEqual([expect.stringMatching(/2 entries, over the bound of 1/)]);
  });

  it("refuses a local entry listed twice", async () => {
    const z = await archive({ "a.txt": enc("a") });
    const { v, central, local } = layout(z);
    v.setUint32(central[1] + 42, local[0], true);
    expect(await said(z)).toEqual([expect.stringMatching(/lists a local entry another lists too/)]);
  });

  it("refuses an entry that is not its declared size or fails its CRC", async () => {
    const stored = await archive({ "a.png": enc("0123456789") }); // .png is stored, not deflated
    const longer = stored.slice();
    const l = layout(longer);
    l.v.setUint32(l.central[1] + 24, 20, true);
    expect(await said(longer)).toEqual([expect.stringMatching(/runs past|holds 10 bytes, declared 20|fails/)]);
    const flipped = stored.slice();
    const f = layout(flipped);
    const name = f.v.getUint16(f.local[1] + 26, true);
    flipped[f.local[1] + 30 + name] ^= 1;
    expect(await said(flipped)).toEqual([expect.stringMatching(/"a\.png" fails its CRC-32/)]);
    const deflated = await archive({ "big.bin": new Uint8Array(1 << 16) });
    const d = layout(deflated);
    d.v.setUint32(d.central[1] + 24, 1000, true);
    expect(await said(deflated)).toEqual([expect.stringMatching(/inflates past its declared 1000 bytes/)]);
  });

  it("reads sizes from the central directory: a local header with zeros (a data descriptor, as Finder writes)", async () => {
    const z = await archive();
    const { v, local } = layout(z);
    v.setUint16(local[0] + 6, v.getUint16(local[0] + 6, true) | 8, true);
    for (const at of [14, 18, 22]) v.setUint32(local[0] + at, 0, true);
    expect(await said(z)).toEqual([]);
  });
});

describe("MessagePack (GI-03, GI-10)", () => {
  /** `{a: <a header of length 0xffffffff, with no payload>}`: 8 bytes. */
  const header = (type: number) => Uint8Array.of(0x81, 0xa1, 0x61, type, 0xff, 0xff, 0xff, 0xff);

  itWithin("a length header with no payload is refused at once, for every form", (within) => {
    for (const type of [0xdd, 0xdf, 0xc6, 0xdb, 0xc9]) {
      within(100, () => expect(() => decodeMsgPack(header(type)), type.toString(16)).toThrow(/truncated/));
    }
  });

  itWithin(
    "nested array headers that each claim every byte left are refused by the bytes read, not the lengths claimed",
    (within) => {
      // 2,000 nested `array32` headers, each claiming 1,000,000 items, then 1 MB of zeros: allocating each claimed array
      // before reading it ran the heap out of memory (the lane A review's LA-01).
      const depth = 2000;
      const claim = 1_000_000;
      const b = new Uint8Array(3 + 5 * depth + claim);
      b.set([0x81, 0xa1, 0x61]);
      const v = new DataView(b.buffer);
      for (let k = 0; k < depth; k++) {
        b[3 + 5 * k] = 0xdd;
        v.setUint32(3 + 5 * k + 1, claim);
      }
      const r = within(2000, () => readDocument(b));
      expect([r.doc, r.issues.map((i) => i.code)]).toEqual([null, ["svs.format.invalid"]]);
      expect(r.issues[0].message).toMatch(/truncated/);
    },
  );

  it("a cut document is svs.format.invalid, not a list of params cut off", () => {
    const whole = toMsgPack(CASES.case1());
    expect(readDocument(whole).issues).toEqual([]);
    for (const at of [0.5, 0.9, 0.99]) {
      const r = readDocument(whole.slice(0, Math.floor(whole.length * at)));
      expect([r.doc, r.issues.map((i) => i.code)]).toEqual([null, ["svs.format.invalid"]]);
      expect(r.issues[0].message).toMatch(/truncated/);
    }
  });

  it("a __proto__ key is an own key, as JSON reads it: the two forms read alike", () => {
    const v = decodeMsgPack(encodeMsgPack(JSON.parse('{"__proto__": {"hidden": 1}, "a": 2}'))) as Record<
      string,
      unknown
    >;
    expect([Object.keys(v), Object.getPrototypeOf(v) === Object.prototype, "hidden" in v]).toEqual([
      ["__proto__", "a"],
      true,
      false,
    ]);
    const doc = JSON.parse(toJSON(CASES.case1()));
    const zinc = doc.root.children[0].children.find((c: { ref?: string }) => c.ref === "zinc");
    zinc.params = JSON.parse(`{"__proto__": ${JSON.stringify(zinc.params)}}`);
    const fromJson = readDocument(new TextEncoder().encode(JSON.stringify(doc)));
    const fromMsgPack = readDocument(encodeMsgPack(doc));
    const reg = kernel();
    const codes = (d: unknown) => validate(d, reg).map((i) => i.code);
    expect(codes(fromMsgPack.doc)).toEqual(codes(fromJson.doc));
    expect(codes(fromJson.doc).length).toBeGreaterThan(0);
  });
});

describe("depth (CM-03)", () => {
  /** `expand` nested `n` deep around a sequence ref. */
  const deep = (n: number): unknown => {
    let e: unknown = "p53";
    for (let k = 0; k < n; k++) e = { head: "seq.modifier.expand", args: { 0: e, 1: 1 } };
    return e;
  };
  const withQuery = (expr: unknown) => {
    const doc = CASES.case1();
    const data = (doc.root.children ?? []).find((c) => c.kind === "data") as { children: unknown[] };
    data.children.push({
      kind: "computed",
      ref: "q",
      params: { provider: "query", query: { seqql: "1.0", expr, text: "(x)" } },
    });
    return doc;
  };

  it("a stored query 20,000 deep is svs.query.invalid, not a stack overflow; 64 deep is read", () => {
    const reg = kernel();
    const at = (n: number) => validate(withQuery(deep(n)), reg).map((i) => `${i.code} ${i.message.slice(0, 16)}`);
    // So deep the document's own bound refuses it first; past 64, the query's.
    expect(at(20_000)).toEqual(["svs.format.invalid the document nes"]);
    expect(at(65)).toEqual(["svs.query.invalid seqql.read.depth"]);
    expect(at(64)).toEqual(["svs.query.stale text is not the "]);
  });

  it("checkQuery refuses one called directly", () => {
    const r = checkQuery(deep(5000) as never);
    expect([r.expr, r.issues.map((i) => i.code)]).toEqual([null, ["seqql.read.depth"]]);
  });

  it("a document nested 20,000 deep is svs.format.invalid at the path where it goes too deep, in either form", () => {
    let custom: unknown = 1;
    for (let k = 0; k < 20_000; k++) custom = { x: custom };
    const doc = CASES.case1() as unknown as { root: { custom?: unknown } };
    doc.root.custom = { acme: custom };
    const reg = kernel();
    const issues = validate(doc, reg);
    expect(issues.map((i) => i.code)).toEqual(["svs.format.invalid"]);
    expect(issues[0].path).toMatch(/^\/root\/custom\/acme(\/x)+$/);
    // JSON.parse takes 20,000 levels (JSON.stringify would not: the text is written by hand).
    const text = JSON.stringify(CASES.case1()).replace(
      '"root":{',
      `"root":{"custom":{"acme":${'{"x":'.repeat(20_000)}1${"}".repeat(20_000)}},`,
    );
    const read = readDocument(new TextEncoder().encode(text));
    expect([read.doc, read.issues.map((i) => i.code)]).toEqual([null, ["svs.format.invalid"]]);
    // 2,000 deep in MessagePack, which decodes, then the same bound
    let shallow: unknown = 1;
    for (let k = 0; k < 2000; k++) shallow = { x: shallow };
    doc.root.custom = { acme: shallow };
    expect(readDocument(encodeMsgPack(doc)).issues.map((i) => i.code)).toEqual(["svs.format.invalid"]);
  });
});
