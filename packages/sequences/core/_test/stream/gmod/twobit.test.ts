/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// 2bit against the recorded slices of hg38.2bit (_test/cases.ts). The expected letters are UCSC's REST API's
// (getData/sequence, which reads the same assembly with its soft mask), recorded in fixtures/expected.json.

import { readFileSync } from "node:fs";
import { StreamNotServed } from "@molstar/seq-core/format";
import { describe, expect, it } from "vitest";
import { CASES, makers, req } from "./cases.js";
import { sliceReader } from "./slices.js";

const expected = JSON.parse(readFileSync(new URL("./fixtures/expected.json", import.meta.url), "utf8"))
  .twobit as Record<string, string>;
const signal = new AbortController().signal;
const provider = () => makers.twobit(sliceReader("twobit"));

describe("2bit", () => {
  it("lists the sequences and their sizes", async () => {
    const p = await provider().fetch(req(null), signal);
    expect(p.kind === "manifest" && p.manifest).toEqual({
      systems: [
        { key: "chr1", extent: [0, 248_956_422], numbering: null, label: "chr1" },
        { key: "chrM", extent: [0, 16_569], numbering: null, label: "chrM" },
      ],
      levels: [],
      alphabet: "dna",
    });
  });

  for (const [name, r] of CASES.twobit.slice(1))
    it(`${name}: the letters, case included, are the API's`, async () => {
      const p = await provider().fetch(r, signal);
      expect(p.kind === "letters" && p.letters).toBe(expected[name]);
    });

  /**
   * A 2bit reader whose reads, once `hold` is set, wait for `release` (or their signal's abort); it honors a signal
   * already aborted, and counts the reads it holds.
   */
  const gated = () => {
    let release = () => {};
    const gate = new Promise<void>((ok) => {
      release = ok;
    });
    const base = sliceReader("twobit");
    const g = { hold: false, held: 0, release: () => release() };
    const read = async (o: number, l: number, s?: AbortSignal) => {
      if (s?.aborted) throw s.reason;
      if (g.hold) {
        g.held++;
        await Promise.race([gate, new Promise((_, no) => s?.addEventListener("abort", () => no(s.reason)))]);
      }
      return base(o, l, s);
    };
    return { g, p: makers.twobit(read) };
  };
  /** Until `ok()`, a macrotask at a time; then a few more, for a request started meanwhile to join the shared read. */
  const until = async (ok: () => boolean) => {
    for (let k = 0; k < 200 && !ok(); k++) await new Promise((t) => setTimeout(t, 1));
    expect(ok()).toBe(true);
    for (let k = 0; k < 5; k++) await new Promise((t) => setTimeout(t, 1));
  };

  it("a request's abort fails no other request on the same sequence (CB-01)", async () => {
    const { g, p } = gated();
    await p.fetch(req(null), signal);
    g.hold = true;
    const [a, b] = [new AbortController(), new AbortController()];
    const first = p.fetch(req("chrM", 0, 40), a.signal);
    // The first request holds the record's read; the second joins it (LA-16: the abort came before either had).
    await until(() => g.held > 0);
    const second = p.fetch(req("chrM", 0, 4096), b.signal);
    await until(() => true);
    a.abort();
    g.release();
    await expect(first).rejects.toThrow(/aborted/);
    const got = await second;
    expect(got.kind === "letters" && got.letters.slice(0, 10)).toBe("GATCACAGGT");
  });

  it("a manifest request's abort fails no data request waiting on the same open (CB-01, `once`)", async () => {
    const { g, p } = gated();
    g.hold = true;
    const [a, b] = [new AbortController(), new AbortController()];
    const manifest = p.fetch(req(null), a.signal);
    await until(() => g.held > 0);
    const data = p.fetch(req("chrM", 0, 40), b.signal);
    await until(() => true);
    a.abort();
    g.release();
    await expect(manifest).rejects.toThrow(/aborted/);
    const got = await data;
    expect(got.kind === "letters" && got.letters.slice(0, 10)).toBe("GATCACAGGT");
  });

  it("reads the N blocks and the soft mask", async () => {
    const p = provider();
    const at = await p.fetch(req("chr1", 9_990, 10_030), signal);
    // 10 of the telomere gap's 10,000 N, then the first repeat, lower case.
    expect(at.kind === "letters" && at.letters).toBe("NNNNNNNNNNtaaccctaaccctaaccctaaccctaaccc");
    const m = await p.fetch(req("chrM", 0, 40), signal);
    expect(m.kind === "letters" && m.letters).toBe("GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCT");
  });

  it("cuts a region at the end of the sequence", async () => {
    const p = await provider().fetch(req("chrM", 16_000, 17_000), signal);
    expect(p.kind === "letters" && p.letters.length).toBe(569);
  });

  it("finds a region of chr1 without reading the mask list (331,893 blocks)", async () => {
    const read = sliceReader("twobit");
    await makers.twobit(read).fetch(req("chr1", 100_000_000, 100_004_096), signal);
    expect(read.served).toBeLessThan(40);
  });

  it("fails clearly for what it does not serve", async () => {
    const p = provider();
    await expect(p.fetch(req("chr1", 0, 10, 0), signal)).rejects.toBeInstanceOf(StreamNotServed);
    await expect(p.fetch(req("chr2", 0, 10), signal)).rejects.toThrow("2bit: no sequence chr2");
    await expect(p.fetch(req("chr1", 300_000_000, 300_000_010), signal)).rejects.toThrow("outside 0..248956422");
  });

  it("fails naming the range when a region is not in the recorded slices", async () => {
    await expect(provider().fetch(req("chr1", 150_000_000, 150_004_096), signal)).rejects.toThrow(
      /twobit: bytes \d+\.\.\d+ are not in the recorded slices/,
    );
  });
});
