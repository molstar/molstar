/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";
import { type MyersMode, myers, myersAlign } from "../../src/algorithms/index.js";

// 1BRS from the viewer spike's fixture (as in aligners.test.ts): entity 1 is barnase's mature chain, P00648 48–157.
const spike = JSON.parse(readFileSync(new URL("../fixtures/1brs.json", import.meta.url), "utf8"));
const P00648 = spike.barnase.sequence as string;
const ENTITY_1 =
  "AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIR";
const DNA = "ACGT";
const AMINO = "ACDEFGHIKLMNPQRSTVWY";
const MODES: MyersMode[] = ["global", "prefix", "infix"];

type Equal = (x: string, y: string) => boolean;
const identical: Equal = (x, y) => x === y;

/** A deterministic generator of integers below n (mulberry32: a linear congruential generator's low bits repeat). */
function generator(seed: number) {
  let s = seed;
  return (n: number) => {
    s = (s + 0x6d2b79f5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) % n;
  };
}

function randomSeq(next: (n: number) => number, length: number, letters: string): string {
  let out = "";
  for (let i = 0; i < length; i++) out += letters[next(letters.length)];
  return out;
}

/** A copy with about `rate` percent of positions substituted, deleted or followed by an insertion. */
function mutate(next: (n: number) => number, seq: string, letters: string, rate: number): string {
  let out = "";
  for (const x of seq) {
    const r = next(300);
    if (r < rate) continue;
    out += r < 2 * rate ? letters[next(letters.length)] : x;
    if (r >= 300 - rate) out += letters[next(letters.length)];
  }
  return out;
}

/** The textbook dynamic program: `prefix` and `infix` leave gaps in `b` after `a` free, `infix` also before it. */
function naiveDistance(a: string, b: string, mode: MyersMode, equal: Equal): number {
  let prev = new Int32Array(b.length + 1);
  let cur = new Int32Array(b.length + 1);
  for (let j = 0; j <= b.length; j++) prev[j] = mode === "infix" ? 0 : j;
  for (let i = 1; i <= a.length; i++) {
    cur[0] = i;
    for (let j = 1; j <= b.length; j++) {
      const diagonal = prev[j - 1] + (equal(a[i - 1], b[j - 1]) ? 0 : 1);
      cur[j] = Math.min(diagonal, prev[j] + 1, cur[j - 1] + 1);
    }
    [prev, cur] = [cur, prev];
  }
  return mode === "global" ? prev[b.length] : Math.min(...prev);
}

/**
 * Replays a CIGAR over both sequences: the letters each side consumes (I consumes a, D consumes b) and the path's
 * cost, the leading D's free in `infix` and the trailing ones in `prefix` and `infix`.
 */
function replay(cigar: string, a: string, b: string, mode: MyersMode, equal: Equal) {
  expect(cigar).toMatch(/^(\d+[MID])+$/);
  const ops = [...cigar.matchAll(/(\d+)([MID])/g)].flatMap(([, n, op]) => Array<string>(Number(n)).fill(op));
  const firstAligned = ops.findIndex((op) => op !== "D");
  const lastAligned = ops.findLastIndex((op) => op !== "D");
  let i = 0;
  let j = 0;
  let cost = 0;
  ops.forEach((op, k) => {
    const free = (mode === "infix" && k < firstAligned) || (mode !== "global" && k > lastAligned);
    if (op === "M") cost += equal(a[i], b[j]) ? 0 : 1;
    else if (!free) cost += 1;
    if (op !== "D") i++;
    if (op !== "I") j++;
  });
  return { a: i, b: j, cost };
}

/** Checks the distance against the naive one and the path against the distance; returns the distance. */
function check(a: string, b: string, mode: MyersMode, pairs: string[] = [], equal: Equal = identical): number {
  const expected = naiveDistance(a, b, mode, equal);
  const result = myersAlign(a, b, mode, null, pairs);
  expect(result?.distance, `${mode} ${a.length}×${b.length}`).toBe(expected);
  const r = replay(result?.cigar as string, a, b, mode, equal);
  expect(r).toEqual({ a: a.length, b: b.length, cost: expected });
  return expected;
}

/** Whether edlib's traceback would store more than 1 MB, so the path comes from Hirschberg's split. */
const hirschberg = (queryLength: number, targetLength: number) =>
  12 * Math.ceil(queryLength / 32) * targetLength + 8 * targetLength >= 1024 * 1024;

describe("myers: edit distances and paths against a naive dynamic program", () => {
  it.each(MODES)("%s, short sequences (one word)", (mode) => {
    const next = generator(1);
    for (let t = 0; t < 200; t++) {
      const letters = t % 2 ? DNA : AMINO;
      const a = randomSeq(next, 1 + next(30), letters);
      // Half related (a mutated copy), half unrelated (large distances, every k doubling).
      const b = t % 4 < 2 ? mutate(next, a, letters, 20) || "A" : randomSeq(next, 1 + next(30), letters);
      check(a, b, mode);
    }
  });

  it.each(MODES)("%s, longer than a word", (mode) => {
    const next = generator(2);
    for (let t = 0; t < 60; t++) {
      const letters = t % 2 ? DNA : AMINO;
      const a = randomSeq(next, 33 + next(300), letters);
      const b =
        t % 3 === 0
          ? randomSeq(next, 33 + next(300), letters)
          : mode === "global"
            ? mutate(next, a, letters, 15)
            : randomSeq(next, next(50), letters) + mutate(next, a, letters, 15) + randomSeq(next, next(50), letters);
      check(a, b, mode);
    }
  });

  it.each(MODES)("%s, longer than Hirschberg's threshold", (mode) => {
    const next = generator(3);
    const a = randomSeq(next, 2100, DNA);
    const b =
      mode === "global"
        ? mutate(next, a, DNA, 20)
        : mode === "prefix"
          ? mutate(next, a, DNA, 10) + randomSeq(next, 900, DNA)
          : randomSeq(next, 1500, DNA) + mutate(next, a, DNA, 10) + randomSeq(next, 1200, DNA);
    expect(hirschberg(a.length, 2000)).toBe(true);
    check(a, b, mode);
    const c = randomSeq(next, 1800, AMINO);
    check(c, mutate(next, c, AMINO, 40), mode);
  });

  it("aligns a query that sits before the target entirely when that is as good as any other", () => {
    expect(check("CCTT", "AAGG", "infix")).toBe(4);
    expect(check("CCTT", "AAGG", "prefix")).toBe(4);
  });
});

describe("myers aligner", () => {
  it("finds barnase's mature chain inside its UniProt precursor", async () => {
    expect(await myers.align(ENTITY_1, P00648, { mode: "infix" })).toBe("47D110M");
    expect(myersAlign(ENTITY_1, P00648, "infix")).toMatchObject({ distance: 0, start: 47, end: 156 });
    // A construct with a tag and two point mutations, found by its distance alone.
    const construct = `MHHHHHH${ENTITY_1.slice(0, 30)}W${ENTITY_1.slice(31, 80)}K${ENTITY_1.slice(81)}`;
    expect(check(construct, P00648, "infix")).toBeLessThanOrEqual(9);
  });

  it("aligns end to end, against the start, or anywhere inside", async () => {
    expect(check("CGTA", "ACGTAGG", "global")).toBe(3);
    expect(check("CGTA", "ACGTAGG", "prefix")).toBe(1);
    expect(check("CGTA", "ACGTAGG", "infix")).toBe(0);
    expect(await myers.align("CGTA", "ACGTAGG", { mode: "infix" })).toBe("1D4M2D");
    expect(await myers.align("acgt", "ACGTTTT", { mode: "prefix" })).toBe("4M3D");
  });

  it("returns null beyond max_distance", async () => {
    expect(await myers.align(P00648, ENTITY_1, { max_distance: 46 })).toBeNull();
    expect(myersAlign(P00648, ENTITY_1, "global", 47)?.distance).toBe(47);
    const construct = `${ENTITY_1.slice(0, 30)}W${ENTITY_1.slice(31)}`;
    expect(await myers.align(construct, P00648, { mode: "infix", max_distance: 0 })).toBeNull();
    expect(await myers.align(construct, P00648, { mode: "infix", max_distance: 1 })).toBe("47D110M");
    expect(await myers.align("ACGTACGT", "TTTTTTTT", { mode: "prefix", max_distance: 3 })).toBeNull();
    for (const mode of MODES) {
      const next = generator(4);
      const a = randomSeq(next, 500, DNA);
      const b = mutate(next, a, DNA, 30);
      const d = naiveDistance(a, b, mode, identical);
      expect(myersAlign(a, b, mode, d - 1)).toBeNull();
      expect(myersAlign(a, b, mode, d)?.distance).toBe(d);
      expect(myersAlign(a, b, mode, d + 100)?.distance).toBe(d);
    }
  });

  it("matches letters named in equalities, both ways and case-insensitive", async () => {
    expect(await myers.align("ACGTNNAC", "ACGTGTAC", { max_distance: 0 })).toBeNull();
    expect(await myers.align("ACGTNNAC", "ACGTGTAC", { max_distance: 0, equalities: "NG NT" })).toBe("8M");
    expect(await myers.align("ACGTGTAC", "acgtnnac", { max_distance: 0, equalities: "ng,tn" })).toBe("8M");
    // Pairs naming letters neither sequence has change nothing.
    expect(await myers.align("ACGT", "ACGA", { equalities: "XY QW" })).toBe("4M");
    const N: Equal = (x, y) => x === y || x === "N" || y === "N";
    const next = generator(5);
    for (const mode of MODES) {
      const a = randomSeq(next, 200, "ACGTN");
      check(a, mutate(next, a, "ACGTN", 20), mode, ["NA", "NC", "NG", "NT"], N);
    }
    await expect(myers.align("ACGT", "ACGT", { equalities: "N" })).rejects.toThrow(/"equalities" must be a string/);
    await expect(myers.align("ACGT", "ACGT", { equalities: 3 })).rejects.toThrow(/"equalities" must be a string/);
  });

  it("refuses empty input and rejects invalid options", async () => {
    expect(await myers.align("", "ACGT")).toBeNull();
    expect(await myers.align("ACGT", "")).toBeNull();
    expect(() => myersAlign("", "ACGT")).toThrow(RangeError);
    await expect(myers.align("ACGT", "ACGT", { mode: "local" })).rejects.toThrow(/"mode" must be one of/);
    await expect(myers.align("ACGT", "ACGT", { max_distance: -1 })).rejects.toThrow(/at least 0/);
    await expect(myers.align("ACGT", "ACGT", { max_distance: 1.5 })).rejects.toThrow(/an integer/);
  });
});
