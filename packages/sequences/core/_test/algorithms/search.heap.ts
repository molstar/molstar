/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for the search matchers (SR1; seqspec:libraries#search "Tests"): a search over a text with hundreds of
// hits, repeated until a million hits have been visited, on RB0b's stepped guard: runs of `n` and `4n` steps, the
// longer growing at most 64 KB beyond the shorter, and no GC in either. `regexMatcher` allocates a match array per hit
// and is not guarded.

import { describe, expect, it } from "vitest";
import {
  approximateMatcher,
  classMatcher,
  compileMotif,
  exactMatcher,
  type Matcher,
  toCodes,
} from "../../src/algorithms/search/index.js";
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from "../model/heap.js";

/** Steps in the shorter run; the longer runs four times as many. */
const N = 2000;
/** A text of 2,000 letters in which every pattern below has 250 to 500 hits: 250,000 to 500,000 over the short run. */
const seq = toCodes("ACGTACGAAC".repeat(200));

let sum = 0;
const visit = (start: number, stop: number, score: number) => {
  sum = (sum + start + stop + score) & 0xfffff;
};

function guard(name: string, matcher: Matcher, cyclic: boolean) {
  it(name, () => {
    let hits = 0;
    const step = () => {
      hits += matcher.find(seq, 0, seq.length, cyclic, visit);
    };
    const r = heapScaling(N, step);
    expect(hits, "the guard visits hits").toBeGreaterThan(N * 200);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // The template's control (CB-05): the same step with a leak is caught.
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });
}

describe("search heap guards", () => {
  guard("exact, linear and across the end of a circle", exactMatcher(toCodes("AC")), true);
  guard("class (a degenerate motif), linear", classMatcher(compileMotif("ARG", { R: "CGT" })), false);
  guard(
    "class across the end of a circle, 40 positions (two words)",
    classMatcher(compileMotif("ACGTACGAAC".repeat(4), {})),
    true,
  );
  guard("approximate, one word, 1 edit", approximateMatcher(toCodes("ACGAA"), 1), false);
  guard("approximate, across the end of a circle", approximateMatcher(toCodes("GTACGAAC"), 2), true);
  guard("approximate, several words, 2 edits", approximateMatcher(toCodes("ACGTACGAAC".repeat(5)), 2), false);
  guard("approximate, every end (ends: all), 2 edits", approximateMatcher(toCodes("ACGAAC"), 2, { ends: "all" }), true);
});
