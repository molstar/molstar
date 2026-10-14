/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { beforeAll, describe, expect, it } from "vitest";
import {
  type ApproximateEnds,
  approximateMatcher,
  classMatcher,
  compileMotif,
  exactMatcher,
  type Matcher,
  regexMatcher,
  textIndex,
  toCodes,
} from "../../src/algorithms/search/index.ts";

/** A deterministic generator of integers below n (mulberry32, as in myers.test.ts). */
function generator(seed: number) {
  let s = seed;
  return (n: number) => {
    s = (s + 0x6d2b79f5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) % n;
  };
}
type Next = (n: number) => number;

const DNA = toCodes("ACGT");
const AMINO = toCodes("ACDEFGHIKLMNPQRSTVWY");
const BINARY = toCodes("AB");
const EXTREME = Uint8Array.of(0, 1, 127, 128, 254, 255);
const ALPHABETS = [DNA, AMINO, BINARY, EXTREME];

type Hit = readonly [start: number, stop: number, score: number];

function run(matcher: Matcher, seq: Uint8Array, from: number, to: number, cyclic: boolean): Hit[] {
  const hits: Hit[] = [];
  const count = matcher.find(seq, from, to, cyclic, (a, b, c) => hits.push([a, b, c]));
  expect(count).toBe(hits.length);
  return hits;
}

const randomSeq = (next: Next, alphabet: Uint8Array, length: number): Uint8Array =>
  Uint8Array.from({ length }, () => alphabet[next(alphabet.length)]);

/** A window of the text, and whether it is read as a circle. */
function randomWindow(next: Next, n: number): [number, number, boolean] {
  const from = next(3) === 0 ? 0 : next(n + 1);
  const to = next(3) === 0 ? n : from + next(n - from + 1);
  return [from, to, next(2) === 0];
}

/** What a hit [start, start + length) reads as in the array's indices: a stop past the end wraps. */
const wrapped = (n: number, start: number, length: number): Hit => [
  start,
  start + length > n ? start + length - n : start + length,
  0,
];

/**
 * The matches of a position test by brute force: every start in the window whose `length` codes (read as a circle when
 * cyclic, a lap at most) pass `test`, ordered by start.
 */
function bruteFixed(
  seq: Uint8Array,
  from: number,
  to: number,
  cyclic: boolean,
  length: number,
  test: (i: number, code: number) => boolean,
): Hit[] {
  const n = seq.length;
  const out: Hit[] = [];
  if (length === 0) return out;
  for (let s = from; s < to; s++) {
    if (cyclic ? length > n : s + length > to) continue;
    let ok = true;
    for (let i = 0; i < length && ok; i++) ok = test(i, seq[(s + i) % n]);
    if (ok) out.push(wrapped(n, s, length));
  }
  return out;
}

describe("exactMatcher", () => {
  it("equals the brute force on random texts, windows and alphabets, linear and cyclic", () => {
    const next = generator(1);
    for (let trial = 0; trial < 3000; trial++) {
      const alphabet = ALPHABETS[trial % ALPHABETS.length];
      const seq = randomSeq(next, alphabet, next(40));
      const pattern =
        seq.length > 1 && next(2) === 0
          ? seq.slice(next(seq.length), seq.length).slice(0, 1 + next(6))
          : randomSeq(next, alphabet, next(6));
      const [from, to, cyclic] = randomWindow(next, seq.length);
      const expected = bruteFixed(seq, from, to, cyclic, pattern.length, (i, c) => c === pattern[i]);
      expect(run(exactMatcher(pattern), seq, from, to, cyclic)).toEqual(expected);
    }
  });

  it("finds overlapping hits, a hit at the very end, and hits across the end of a circle", () => {
    const aaa = toCodes("AAA");
    expect(run(exactMatcher(toCodes("AA")), aaa, 0, 3, false)).toEqual([
      [0, 2, 0],
      [1, 3, 0],
    ]);
    expect(run(exactMatcher(toCodes("GT")), toCodes("ACGT"), 0, 4, false)).toEqual([[2, 4, 0]]);
    // The circle AACGTA reads "TAA" from 4, across the end, to 1.
    expect(run(exactMatcher(toCodes("TAA")), toCodes("AACGTA"), 0, 6, true)).toEqual([[4, 1, 0]]);
    expect(run(exactMatcher(toCodes("TAA")), toCodes("AACGTA"), 0, 6, false)).toEqual([]);
    // Both ways: the same circle read from the other side of the origin.
    expect(run(exactMatcher(toCodes("AAC")), toCodes("ACGTAA"), 0, 6, true)).toEqual([[5, 2, 0]]);
    // As long as the circle: start == stop.
    expect(run(exactMatcher(toCodes("CGA")), toCodes("ACG"), 0, 3, true)).toEqual([[1, 1, 0]]);
  });

  it("has no hits for an empty pattern, a pattern longer than the text, or an empty window", () => {
    const seq = toCodes("ACGT");
    expect(run(exactMatcher(new Uint8Array(0)), seq, 0, 4, true)).toEqual([]);
    expect(run(exactMatcher(toCodes("ACGTA")), seq, 0, 4, false)).toEqual([]);
    expect(run(exactMatcher(toCodes("ACGTA")), seq, 0, 4, true)).toEqual([]);
    expect(run(exactMatcher(toCodes("A")), seq, 2, 2, false)).toEqual([]);
    expect(run(exactMatcher(toCodes("A")), new Uint8Array(0), 0, 0, true)).toEqual([]);
  });

  it("copies its pattern, and rejects a window outside the array", () => {
    const pattern = toCodes("AC");
    const matcher = exactMatcher(pattern);
    pattern[0] = 71;
    expect(run(matcher, toCodes("ACGAC"), 0, 5, false)).toEqual([
      [0, 2, 0],
      [3, 5, 0],
    ]);
    expect(() => matcher.find(toCodes("ACG"), 0, 4, false, () => {})).toThrow(RangeError);
    expect(() => matcher.find(toCodes("ACG"), 2, 1, false, () => {})).toThrow(RangeError);
    expect(() => matcher.find(toCodes("ACG"), 0.5, 2, false, () => {})).toThrow(RangeError);
  });
});

/** A random set of codes of the alphabet: one code, or a few. */
const randomSet = (next: Next, alphabet: Uint8Array): Uint8Array =>
  Uint8Array.from(new Set(Array.from({ length: 1 + next(3) }, () => alphabet[next(alphabet.length)])));

describe("classMatcher", () => {
  it("equals the brute force with sets per position, patterns up to 64, linear and cyclic", () => {
    const next = generator(2);
    for (let trial = 0; trial < 3000; trial++) {
      const alphabet = ALPHABETS[trial % ALPHABETS.length];
      const m = trial % 7 === 0 ? 30 + next(35) : next(7);
      const pattern = Array.from({ length: m }, () => randomSet(next, alphabet));
      const seq = randomSeq(next, alphabet, next(90));
      // Plant a match, so the long patterns hit.
      if (seq.length > m && m > 0 && next(2) === 0) {
        const at = next(seq.length - m + 1);
        for (let i = 0; i < m; i++) seq[at + i] = pattern[i][next(pattern[i].length)];
      }
      const [from, to, cyclic] = randomWindow(next, seq.length);
      const expected = bruteFixed(seq, from, to, cyclic, m, (i, c) => pattern[i].includes(c));
      expect(run(classMatcher(pattern), seq, from, to, cyclic)).toEqual(expected);
    }
  });

  it("compiles IUPAC letters from the caller's classes, and finds the degenerate motif", () => {
    const iupac = { R: "AG", Y: "CT", N: "ACGT" };
    const motif = compileMotif("GRNCY", iupac);
    expect(motif.map((set) => Array.from(set))).toEqual([[71], [65, 71], [65, 67, 71, 84], [67], [67, 84]]);
    // GAACC and GGTCT match; GCACC does not (R is A or G).
    expect(run(classMatcher(motif), toCodes("GAACCTGGTCTGCACC"), 0, 16, false)).toEqual([
      [0, 5, 0],
      [6, 11, 0],
    ]);
    // A letter without a class stands for itself, a class does not include its own letter, and it is case-sensitive.
    expect(compileMotif("rR", { R: "AG" }).map((set) => Array.from(set))).toEqual([[114], [65, 71]]);
    expect(() => compileMotif("é€", {})).toThrow(RangeError);
    expect(compileMotif("", iupac)).toEqual([]);
  });

  it("limits a pattern to 64 positions, and an empty pattern has no hits", () => {
    expect(() => classMatcher(Array.from({ length: 65 }, () => Uint8Array.of(65)))).toThrow(RangeError);
    const max = classMatcher(Array.from({ length: 64 }, () => Uint8Array.of(65)));
    expect(run(max, new Uint8Array(70).fill(65), 0, 70, false).length).toBe(7);
    expect(run(classMatcher([]), toCodes("ACGT"), 0, 4, true)).toEqual([]);
  });
});

describe("regexMatcher", () => {
  const expressions = [/A[CG]+T/, /(AC)+/, /G.A/, /^A/, /T$/, /A{2,3}/, /(?<=G)C/, /.{0}/, /[AC]|GG/];
  it("equals a sticky test at every start, linear and cyclic, overlapping", () => {
    const next = generator(3);
    for (let trial = 0; trial < 1500; trial++) {
      const re = expressions[trial % expressions.length];
      const seq = randomSeq(next, DNA, next(40));
      const [from, to, cyclic] = randomWindow(next, seq.length);
      const maxLength = 8;
      const n = seq.length;
      // The text a search reads: the window, and maxLength - 1 codes more across the end of a circle (a lap at most).
      const readable = cyclic ? Math.min(to - from + maxLength - 1, to - from + n) : to - from;
      const text = String.fromCharCode(...Array.from({ length: readable }, (_, i) => seq[(from + i) % n]));
      const sticky = new RegExp(re.source, "y");
      const found: [number, number][] = [];
      for (let s = from; s < to; s++) {
        sticky.lastIndex = s - from;
        const m = sticky.exec(text);
        if (m && m[0].length > 0 && m[0].length <= n) found.push([s, m[0].length]);
      }
      // In order of their ends (the Matcher contract), then of their starts.
      found.sort((a, b) => a[0] + a[1] - (b[0] + b[1]) || a[0] - b[0]);
      const expected: Hit[] = found.map(([s, l]) => wrapped(n, s, l));
      expect(run(regexMatcher(re, maxLength), seq, from, to, cyclic)).toEqual(expected);
    }
  });

  it("reads anchors against the window, counts one hit per start, and keeps the caller's expression unchanged", () => {
    const re = /A/y;
    expect(run(regexMatcher(re), toCodes("CAA"), 1, 3, false)).toEqual([
      [1, 2, 0],
      [2, 3, 0],
    ]);
    expect(re.lastIndex).toBe(0);
    expect(run(regexMatcher(/^A/), toCodes("CAA"), 1, 3, false)).toEqual([[1, 2, 0]]);
    expect(run(regexMatcher(/AAC/), toCodes("ACGTAA"), 0, 6, true)).toEqual([[5, 2, 0]]);
    // A cyclic search reads maxLength - 1 codes past the window, so a maxLength below 1 would cut the window short.
    expect(() => regexMatcher(/A/, 0)).toThrow(RangeError);
  });

  it("visits hits in order of their ends, as every matcher does (CB-04)", () => {
    expect(run(regexMatcher(/A.{3}|X/), toCodes("AXXXY"), 0, 5, false)).toEqual([
      [1, 2, 0],
      [2, 3, 0],
      [0, 4, 0],
      [3, 4, 0],
    ]);
  });

  it("visits many overlapping hits of mixed lengths in order of their ends, then starts (LA-09's heap)", () => {
    let x = 7;
    const rand = () => {
      x = (x * 1103515245 + 12345) % 2147483648;
      return x / 2147483648;
    };
    const text = Array.from({ length: 4000 }, () => "ACGT"[Math.floor(rand() * 4)]).join("");
    for (const re of [/A[CGT]{0,12}T/, /A.{0,40}?G|C.{5}|T/, /(?:AC|G){1,9}/]) {
      const hits = run(regexMatcher(re), toCodes(text), 0, text.length, false);
      const sorted = [...hits].sort((a, b) => a[1] - b[1] || a[0] - b[0]);
      expect(hits.length).toBeGreaterThan(100);
      expect(hits).toEqual(sorted);
    }
  });

  it("is exact over 3 million codes in chunks", () => {
    const seq = new Uint8Array(3_000_000).fill(67);
    seq[2_999_999] = 71;
    seq[8190] = 71;
    expect(run(regexMatcher(/G/), seq, 0, seq.length, false).map((h) => h[0])).toEqual([8190, 2_999_999]);
  });
});

/** Edit distance of two code arrays by the textbook DP, `equal` telling when two codes match. */
function editDistance(a: ArrayLike<number>, b: ArrayLike<number>, equal: (x: number, y: number) => boolean): number {
  let prev = Array.from({ length: b.length + 1 }, (_, j) => j);
  for (let i = 1; i <= a.length; i++) {
    const cur = [i];
    for (let j = 1; j <= b.length; j++) {
      cur[j] = Math.min(prev[j - 1] + (equal(a[i - 1], b[j - 1]) ? 0 : 1), prev[j] + 1, cur[j - 1] + 1);
    }
    prev = cur;
  }
  return prev[b.length];
}

/**
 * The approximate matcher's specification by brute force: over the text (the window, or a circle read on past `to` by
 * a lap), the cost of an end is the least edit distance of the pattern to a stretch ending there. With `"all"` an end
 * is a hit when its cost is at most k; with the local-best policies also when its cost is 0, or when it is the last
 * (`"local_best_last"`) or the first (`"local_best_first"`) end of a run of ends of that cost with a higher cost before
 * the run (the text's start: m) and after it (or the text's end). The hit is the shortest stretch of that cost, at most
 * a lap long; it counts when it starts in the window.
 */
function bruteApproximate(
  seq: Uint8Array,
  pattern: Uint8Array,
  k: number,
  equalPairs: readonly number[],
  from: number,
  to: number,
  cyclic: boolean,
  ends: ApproximateEnds,
): Hit[] {
  const n = seq.length;
  const equal = (x: number, y: number) => {
    if (x === y) return true;
    for (let i = 0; i < equalPairs.length; i += 2) {
      if ((equalPairs[i] === x && equalPairs[i + 1] === y) || (equalPairs[i] === y && equalPairs[i + 1] === x))
        return true;
    }
    return false;
  };
  const text = Array.from({ length: (cyclic ? to + n : to) - from }, (_, i) => seq[(from + i) % n]);
  const m = pattern.length;
  const out: Hit[] = [];
  if (m === 0) return out;
  // The infix DP: the first row is free, so each column's last cell is the cost of the best stretch ending there.
  const cost: number[] = [];
  let column = Array.from({ length: m + 1 }, (_, i) => i);
  for (const c of text) {
    const next = [0];
    for (let i = 1; i <= m; i++) {
      next[i] = Math.min(column[i - 1] + (equal(pattern[i - 1], c) ? 0 : 1), column[i] + 1, next[i - 1] + 1);
    }
    cost.push(next[m]);
    column = next;
  }
  for (let j = 0; j < text.length; j++) {
    if (cost[j] > k) continue;
    if (ends !== "all" && cost[j] > 0) {
      // The run of ends of this cost around j: [first, after).
      let first = j;
      while (first > 0 && cost[first - 1] === cost[j]) first--;
      let after = j + 1;
      while (after < text.length && cost[after] === cost[j]) after++;
      const before = first > 0 ? cost[first - 1] : m;
      if (before < cost[j] || (after < text.length && cost[after] < cost[j])) continue;
      if (j !== (ends === "local_best_first" ? first : after - 1)) continue;
    }
    let s = j;
    while (s >= 0 && editDistance(pattern, text.slice(s, j + 1), equal) !== cost[j]) s--;
    if (s < 0) throw new Error("no stretch has the end's cost");
    if (j + 1 - s > n || from + s >= to) continue;
    const stop = from + j + 1;
    out.push([from + s, stop > n ? stop - n : stop, cost[j]]);
  }
  return out;
}

/** A copy of the pattern with some edits. */
function mutate(next: Next, seq: Uint8Array, alphabet: Uint8Array, edits: number): Uint8Array {
  const out = Array.from(seq);
  for (let e = 0; e < edits && out.length > 0; e++) {
    const at = next(out.length);
    const kind = next(3);
    if (kind === 0) out[at] = alphabet[next(alphabet.length)];
    else if (kind === 1) out.splice(at, 1);
    else out.splice(at, 0, alphabet[next(alphabet.length)]);
  }
  return Uint8Array.from(out);
}

describe("approximateMatcher", () => {
  it("equals the dynamic-programming oracle on planted and random patterns, up to 3 edits, linear and cyclic", () => {
    const next = generator(4);
    let hits = 0;
    let withEdits = 0;
    let crossing = 0;
    let allHits = 0;
    let firstDiffers = 0;
    for (let trial = 0; trial < 1500; trial++) {
      const alphabet = ALPHABETS[trial % ALPHABETS.length];
      const seq = randomSeq(next, alphabet, next(70));
      const m = 1 + next(trial % 5 === 0 ? 12 : 8);
      const planted = seq.length > m && next(3) > 0;
      const pattern = planted
        ? mutate(next, seq.slice(next(seq.length - m + 1)).slice(0, m), alphabet, next(4))
        : randomSeq(next, alphabet, m);
      if (pattern.length === 0) continue;
      const k = Math.min(next(4), pattern.length - 1);
      const [from, to, cyclic] = randomWindow(next, seq.length);
      const best = bruteApproximate(seq, pattern, k, [], from, to, cyclic, "local_best_last");
      expect(run(approximateMatcher(pattern, k), seq, from, to, cyclic), `trial ${trial}`).toEqual(best);
      const first = bruteApproximate(seq, pattern, k, [], from, to, cyclic, "local_best_first");
      expect(
        run(approximateMatcher(pattern, k, { ends: "local_best_first" }), seq, from, to, cyclic),
        `trial ${trial}`,
      ).toEqual(first);
      const all = bruteApproximate(seq, pattern, k, [], from, to, cyclic, "all");
      expect(run(approximateMatcher(pattern, k, { ends: "all" }), seq, from, to, cyclic), `trial ${trial}`).toEqual(
        all,
      );
      hits += best.length;
      withEdits += best.filter((h) => h[2] > 0).length;
      crossing += best.filter((h) => h[0] > h[1]).length;
      allHits += all.length;
      firstDiffers += best.some((h, i) => first[i]?.[1] !== h[1]) ? 1 : 0;
    }
    // The comparison was not vacuous, and the policies differ.
    expect(hits).toBeGreaterThan(2000);
    expect(withEdits).toBeGreaterThan(1000);
    expect(crossing).toBeGreaterThan(30);
    expect(allHits).toBeGreaterThan(hits * 1.5);
    expect(firstDiffers).toBeGreaterThan(100);
  });

  /** `trials` seeded patterns of 31 to 70 letters, each against the oracle at each `ends`. */
  function longPatterns(trials: number) {
    const next = generator(5);
    for (let trial = 0; trial < trials; trial++) {
      const alphabet = trial % 2 === 0 ? DNA : AMINO;
      const seq = randomSeq(next, alphabet, 80 + next(80));
      // A third at the blocks' edges (a full word: the score at the sign bit), the rest anywhere over a word.
      const m = trial % 3 === 0 ? [31, 32, 33, 63, 64, 65][next(6)] : 33 + next(38);
      const pattern = mutate(next, seq.slice(next(seq.length - m)).slice(0, m), alphabet, next(4));
      // Mostly the few edits of a search, sometimes many.
      const k = Math.min(trial % 4 === 0 ? next(12) : next(4), pattern.length - 1);
      const [from, to, cyclic] = randomWindow(next, seq.length);
      for (const ends of ["local_best_last", "local_best_first", "all"] as const) {
        const expected = bruteApproximate(seq, pattern, k, [], from, to, cyclic, ends);
        expect(run(approximateMatcher(pattern, k, { ends }), seq, from, to, cyclic), `trial ${trial} ${ends}`).toEqual(
          expected,
        );
      }
    }
  }

  it("equals the oracle for 31 to 70 letters (one to three blocks, full words too), 40 patterns", () =>
    longPatterns(40));

  it(
    "equals the oracle for 31 to 70 letters (one to three blocks, full words too)",
    { tags: ["slow"], timeout: 60000 },
    () => longPatterns(120),
  );

  it("counts a pattern's letters that match others (`equal`, both ways)", () => {
    const next = generator(6);
    const N = 78;
    const pairs = [N, 65, N, 67, N, 71, N, 84];
    for (let trial = 0; trial < 400; trial++) {
      const seq = randomSeq(next, Uint8Array.of(65, 67, 71, 84, N), 10 + next(50));
      const pattern = mutate(next, seq.slice(next(8)).slice(0, 4 + next(8)), Uint8Array.of(65, 67, 71, 84, N), next(3));
      if (pattern.length < 2) continue;
      const k = Math.min(next(3), pattern.length - 1);
      const [from, to, cyclic] = randomWindow(next, seq.length);
      const ends = (["local_best_last", "local_best_first", "all"] as const)[trial % 3];
      const expected = bruteApproximate(seq, pattern, k, pairs, from, to, cyclic, ends);
      expect(
        run(approximateMatcher(pattern, k, { equal: Uint8Array.from(pairs), ends }), seq, from, to, cyclic),
        `trial ${trial}`,
      ).toEqual(expected);
    }
  });

  it("with no edits equals the exact matcher, hits overlapping and across the end of a circle", () => {
    const next = generator(7);
    for (let trial = 0; trial < 500; trial++) {
      const seq = randomSeq(next, trial % 2 === 0 ? BINARY : DNA, 1 + next(50));
      const pattern = randomSeq(next, trial % 2 === 0 ? BINARY : DNA, 1 + next(4));
      const [from, to, cyclic] = randomWindow(next, seq.length);
      expect(run(approximateMatcher(pattern, 0), seq, from, to, cyclic)).toEqual(
        run(exactMatcher(pattern), seq, from, to, cyclic),
      );
    }
  });

  it("reports one hit for a site and its edits, with the stretch and the cost", () => {
    // ACGT with one extra G: the ends from ACG to ACGGT cost 1; the last is the hit (by default), the insertion.
    expect(run(approximateMatcher(toCodes("ACGT"), 1), toCodes("TTACGGTTT"), 0, 9, false)).toEqual([[2, 7, 1]]);
    // An exact site beats the neighbors that cost 1.
    expect(run(approximateMatcher(toCodes("ACGT"), 1), toCodes("TTACGTTT"), 0, 8, false)).toEqual([[2, 6, 0]]);
    // Across the end of a circle: GT|AC with one substitution.
    expect(run(approximateMatcher(toCodes("GTTC"), 1), toCodes("ACCCCCGT"), 0, 8, true)).toEqual([[6, 2, 1]]);
    expect(run(approximateMatcher(toCodes("GTTC"), 1), toCodes("ACCCCCGT"), 0, 8, false)).toEqual([]);
    // A run of ends of cost 1 that falls to an exact site is not a minimum: one hit, not AAA with a deletion as well.
    expect(run(approximateMatcher(toCodes("AAAC"), 1), toCodes("ACAAAACA"), 0, 8, false)).toEqual([[3, 7, 0]]);
  });

  it('with ends: "all" reports every end of cost at most k, each as its shortest stretch', () => {
    // The ends of AAA from 2 and from 3 (a deletion: shorter than AAAA, a substitution), AAAC exact, AAACA (an insertion).
    expect(run(approximateMatcher(toCodes("AAAC"), 1, { ends: "all" }), toCodes("ACAAAACA"), 0, 8, false)).toEqual([
      [2, 5, 1],
      [3, 6, 1],
      [3, 7, 0],
      [3, 8, 1],
    ]);
    // A run of equal cost: local-best gives its last end or its first, all every end.
    const run1 = (ends: ApproximateEnds) =>
      run(approximateMatcher(toCodes("ACGT"), 1, { ends }), toCodes("TTACGGTTT"), 0, 9, false);
    expect(run1("local_best_last")).toEqual([[2, 7, 1]]);
    expect(run1("local_best_first")).toEqual([[2, 5, 1]]);
    expect(run1("all").map((h) => h[1])).toEqual([5, 6, 7]);
  });

  it("takes a run's last end or its first: a mismatch at the pattern's end as a substitution or a deletion", () => {
    const at = (ends: ApproximateEnds, pattern: string, text: string) =>
      run(approximateMatcher(toCodes(pattern), 1, { ends }), toCodes(text), 0, text.length, false);
    // The story's cases: RHSWV over RHSVVV, ATGGCC over ATGGCT.
    expect(at("local_best_last", "RHSWV", "TFRHSVVVPY")).toEqual([[2, 7, 1]]);
    expect(at("local_best_first", "RHSWV", "TFRHSVVVPY")).toEqual([[2, 6, 1]]);
    expect(at("local_best_last", "ATGGCC", "CCATGGCTAA")).toEqual([[2, 8, 1]]);
    expect(at("local_best_first", "ATGGCC", "CCATGGCTAA")).toEqual([[2, 7, 1]]);
    // A deletion of the pattern's last letter: the first end reads it so, the last a substitution after it.
    expect(at("local_best_first", "ACGTA", "GGACGTCC")).toEqual([[2, 6, 1]]);
    expect(at("local_best_last", "ACGTA", "GGACGTCC")).toEqual([[2, 7, 1]]);
  });

  it("has no hits for an empty pattern, a pattern longer than the text, or an empty window, and checks its options", () => {
    expect(run(approximateMatcher(new Uint8Array(0), 2), toCodes("ACGT"), 0, 4, true)).toEqual([]);
    expect(run(approximateMatcher(toCodes("ACGTACGT"), 1), toCodes("ACGT"), 0, 4, false)).toEqual([]);
    expect(run(approximateMatcher(toCodes("A"), 0), toCodes("ACGT"), 2, 2, false)).toEqual([]);
    expect(approximateMatcher(toCodes("ACGT"), 3).length).toBe(7);
    expect(() => approximateMatcher(toCodes("ACGT"), 4)).toThrow(RangeError);
    expect(() => approximateMatcher(toCodes("ACGT"), -1)).toThrow(RangeError);
    expect(() => approximateMatcher(toCodes("ACGT"), 1.5)).toThrow(RangeError);
    expect(() => approximateMatcher(toCodes("ACGT"), 1, { equal: Uint8Array.of(65) })).toThrow(RangeError);
    expect(() => approximateMatcher(toCodes("ACGT"), 1, { ends: "local-best" as ApproximateEnds })).toThrow(RangeError);
  });
});

describe("textIndex", () => {
  /** The reference: a plain array of [key, id], sorted by key then id. */
  class Reference {
    pairs: [string, number][] = [];
    add(id: number, text: string) {
      const key = text.normalize("NFKC").toLowerCase();
      // GM-08: a text added again for its id is held once.
      if (!this.pairs.some(([k, i]) => k === key && i === id)) this.pairs.push([key, id]);
    }
    remove(ids: Set<number>) {
      this.pairs = this.pairs.filter(([, id]) => !ids.has(id));
    }
    sorted() {
      return [...this.pairs].sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : a[1] - b[1]));
    }
  }
  const collect = (fn: (visit: (id: number) => void) => number): number[] => {
    const out: number[] = [];
    const count = fn((id) => out.push(id));
    expect(count).toBe(out.length);
    return out;
  };

  it("interleaves add, remove and queries against a sorted array", () => {
    const next = generator(8);
    const words = [
      "",
      "a",
      "ab",
      "abc",
      "b",
      "Ab",
      "ABC",
      "été",
      "Été",
      "Ångström",
      "ångström",
      "ﬁsh",
      "ｆｉｓｈ",
      "z",
    ];
    const index = textIndex();
    const ref = new Reference();
    let nextId = 1;
    for (let step = 0; step < 400; step++) {
      const op = next(10);
      if (op < 5) {
        const count = 1 + next(5);
        const ids = Int32Array.from({ length: count + 2 }, () => nextId++ % 40);
        const texts = Array.from({ length: count + 2 }, () => words[next(words.length)] + (next(3) === 0 ? "x" : ""));
        index.add(ids, texts, 1, count);
        for (let i = 1; i < 1 + count; i++) ref.add(ids[i], texts[i]);
      } else if (op < 7) {
        const ids = Int32Array.from({ length: 4 }, () => next(40));
        const count = 1 + next(3);
        index.remove(ids, count);
        ref.remove(new Set(ids.subarray(0, count)));
      } else {
        const q = words[next(words.length)];
        const key = q.normalize("NFKC").toLowerCase();
        const sorted = ref.sorted();
        expect(collect((v) => index.exact(q, v))).toEqual(sorted.filter(([k]) => k === key).map(([, id]) => id));
        expect(collect((v) => index.prefix(q, v))).toEqual(
          sorted.filter(([k]) => k.startsWith(key)).map(([, id]) => id),
        );
        const re = /[bz]$/;
        expect(collect((v) => index.matches(re, v))).toEqual(sorted.filter(([k]) => re.test(k)).map(([, id]) => id));
      }
      expect(index.size).toBe(ref.pairs.length);
    }
  });

  it("normalizes keys and queries: compatibility forms, case and accents of precomposed letters", () => {
    const index = textIndex();
    index.add(Int32Array.of(1, 2, 3, 4), ["Ångström", "ＡＢＣ", "ﬁsh", "Zeta"], 0, 4);
    expect(collect((v) => index.prefix("ångs", v))).toEqual([1]);
    expect(collect((v) => index.prefix("Ångs", v))).toEqual([1]);
    expect(collect((v) => index.exact("abc", v))).toEqual([2]);
    expect(collect((v) => index.prefix("FI", v))).toEqual([3]);
    // In order of the keys' code units: "ångström" comes last.
    expect(collect((v) => index.prefix("", v))).toEqual([2, 3, 4, 1]);
    expect(collect((v) => index.exact("ab", v))).toEqual([]);
    expect(collect((v) => index.matches(/e.a$/g, v))).toEqual([4]);
  });

  it("takes the caller's normalization, and keeps several texts under one id", () => {
    const index = textIndex({ normalize: (s) => s });
    index.add(Int32Array.of(7, 7, 8), ["Gene", "gene", "Gene"], 0, 3);
    expect(collect((v) => index.exact("Gene", v))).toEqual([7, 8]);
    expect(collect((v) => index.exact("gene", v))).toEqual([7]);
    index.remove(Int32Array.of(7), 1);
    expect(collect((v) => index.prefix("", v))).toEqual([8]);
    index.remove(Int32Array.of(8), 0);
    expect(index.size).toBe(1);
  });

  it("GM-08: a text added again for its id (a region fetched again) is held once", () => {
    const index = textIndex();
    index.add(Int32Array.of(1, 2, 1), ["BRCA1", "TP53", "brca1"], 0, 3);
    expect(index.size).toBe(2);
    index.add(Int32Array.of(1), ["BRCA1"], 0, 1);
    expect(collect((v) => index.exact("brca1", v))).toEqual([1]);
    expect(index.size).toBe(2);
    index.remove(Int32Array.of(1), 1);
    expect(collect((v) => index.prefix("", v))).toEqual([2]);
  });
});

// Wall-clock bounds over 3 million nucleotides: the slow tier (the oracles above check the matchers' results).
describe("large inputs", { tags: ["slow"] }, () => {
  const N = 3_000_000;
  let seq: Uint8Array = new Uint8Array(0);
  let pattern = seq;
  beforeAll(() => {
    seq = randomSeq(generator(9), DNA, N);
    pattern = seq.slice(1_234_567, 1_234_587);
  });
  const timed = (matcher: Matcher): { hits: Hit[]; ms: number } => {
    const hits: Hit[] = [];
    const t = performance.now();
    matcher.find(seq, 0, N, false, (a, b, c) => {
      if (hits.length < 100) hits.push([a, b, c]);
    });
    return { hits, ms: performance.now() - t };
  };

  it("finds a 20-mer in 3 million nucleotides, exactly, by class and with 2 edits, in a bounded time", () => {
    const exact = timed(exactMatcher(pattern));
    expect(exact.hits).toContainEqual([1_234_567, 1_234_587, 0]);
    const sets = Array.from(pattern, (c) => Uint8Array.of(c));
    const cls = timed(classMatcher(sets));
    expect(cls.hits).toEqual(exact.hits);
    const approx = timed(approximateMatcher(pattern, 2));
    expect(approx.hits).toContainEqual([1_234_567, 1_234_587, 0]);
    expect(approx.hits.every((h) => h[2] <= 2)).toBe(true);
    // A laptop does these in tens of milliseconds; the bound is for a busy machine.
    expect(exact.ms).toBeLessThan(1000);
    expect(cls.ms).toBeLessThan(1000);
    expect(approx.ms).toBeLessThan(2000);
  });

  it("searches a 60-letter pattern with 3 edits (two blocks) in 3 million nucleotides in a bounded time", () => {
    const long = seq.slice(2_000_000, 2_000_060);
    const approx = timed(approximateMatcher(long, 3));
    expect(approx.hits).toContainEqual([2_000_000, 2_000_060, 0]);
    expect(approx.ms).toBeLessThan(4000);
  });
});
