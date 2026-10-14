/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from "vitest";
import { MSAMSAAlignment, MSASeqAlignment } from "../../src/algorithms/biomsa/align/align.js";
import { AlignerOptionsError, AlignerSizeError, biomsa, MAX_CELLS } from "../../src/algorithms/index.js";

const many = (seqs: string[], options = {}, signal?: AbortSignal) =>
  (biomsa.alignMany as NonNullable<typeof biomsa.alignMany>)(seqs, options, signal);

const AA = "ACDEFGHIKLMNPQRSTVWY";
function family(n: number, length: number): string[] {
  let seed = 3;
  const rnd = () => {
    seed = (Math.imul(seed, 1103515245) + 12345) & 0x7fffffff;
    return seed / 0x7fffffff;
  };
  const base = Array.from({ length }, () => AA[Math.floor(rnd() * 20)]).join("");
  return Array.from({ length: n }, () => [...base].map((c) => (rnd() < 0.2 ? AA[Math.floor(rnd() * 20)] : c)).join(""));
}

describe("biomsa alignMany", () => {
  it("aligns identical-but-for-an-end sequences into equal-length rows in input order", async () => {
    const rows = await many(["MKVLAGW", "MKVLAGW", "mkvlagw", "VLAGW"]);
    expect(rows).toEqual(["MKVLAGW", "MKVLAGW", "MKVLAGW", "--VLAGW"]);
    expect(await many(["ACGTTGCAAC", "ACGTTGCAAC", "ACGTTG"])).toEqual(["ACGTTGCAAC", "ACGTTGCAAC", "ACGTTG----"]);
  });
  it("aligns two sequences as align does", async () => {
    const [a, b] = ["MKVLAGW", "MKVAGW"];
    const rows = (await many([a, b])) as string[];
    expect(rows.map((r) => r.replace(/-/g, ""))).toEqual([a, b]);
    expect(rows[0].length).toBe(rows[1].length);
  });
  it("is deterministic, and the rows keep every residue of their sequence", async () => {
    const seqs = family(12, 60).map((s, i) => s.slice(i % 5));
    const [x, y] = [await many(seqs), await many(seqs)];
    expect(x).toEqual(y);
    expect(x?.map((r) => r.replace(/-/g, ""))).toEqual(seqs);
    expect(new Set(x?.map((r) => r.length)).size).toBe(1);
  });
  it("is null for sequences it cannot align, and an options error for a bad request", async () => {
    expect(await many(["MKV", "M-V", "MKV"])).toBeNull();
    expect(await many(["MKVLAGW", "ACGTACGTAC", "MKVLAGW"])).toBeNull(); // protein and DNA together
    await expect(many(["MKV"])).rejects.toBeInstanceOf(AlignerOptionsError);
    await expect(many(["MKV", "MKV"], { gap_open: -1 })).rejects.toBeInstanceOf(AlignerOptionsError);
  });
  it("rejects with the signal's reason before starting, and while a long run waits between its steps", async () => {
    await expect(many(["MKV", "MKV"], {}, AbortSignal.abort(new Error("newer")))).rejects.toThrow("newer");
    const c = new AbortController();
    const running = many(family(120, 300), {}, c.signal);
    setTimeout(() => c.abort(new Error("newer request")), 1);
    await expect(running).rejects.toThrow("newer request");
  });
  it("runs two multiple alignments at once without mixing their sequences", async () => {
    const [a, b] = [family(40, 120), family(40, 130).map((s) => s.toLowerCase())];
    const [ra, rb, ref] = await Promise.all([many(a), many(b), many(a)]);
    expect(ra).toEqual(ref);
    expect(ra?.map((r) => r.replace(/-/g, ""))).toEqual(a);
    expect(rb?.map((r) => r.replace(/-/g, ""))).toEqual(b.map((s) => s.toUpperCase()));
  });
});

describe("the profile alignments", () => {
  // B15: a profile column left unpaired at the start was rewarded; four sequences of one family, the last with a
  // different start, aligned with two gap columns before.
  const ROWS = [
    "CTSMTGNQYHQKFFDPCDPNTENEMVDNANLNAVCMHHAT",
    "CCSATLNQYHTKFYDGCDPNTENEMVDNNTLNAVAMGHAF",
    "CCSATLNQYHTKDGGPCDPNTENEMTDNPNLNAVCMGHAF",
    "RFFIKWQQYHTKFYAPCDPNTENEMVDNNNIFAVCMGHAF",
  ];
  it("pairs a sequence with a profile's first columns (no gap columns of nothing)", async () => {
    const rows = (await many(ROWS)) as string[];
    expect(rows).toEqual(ROWS);
  });
  it("refuses a profile and a sequence, or two profiles, past MAX_CELLS before allocating", () => {
    const node = (length: number) => ({ profile: { length }, seq: { rawSeq: { length } } }) as never;
    const side = Math.ceil(Math.sqrt(MAX_CELLS)) + 1;
    expect(() => MSASeqAlignment(node(side), node(side), {} as never)).toThrow(AlignerSizeError);
    expect(() => MSAMSAAlignment(node(side), node(side), {} as never)).toThrow(AlignerSizeError);
    expect(() => MSASeqAlignment(node(side), node(side), {} as never)).toThrow(/limit is 67108864/);
  });
});
