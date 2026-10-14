/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Marks (IN1; seqspec:interaction#channels, the marks; CP1a): a MarkSet's runs per block, sorted and disjoint, merged
// as they are added (in order, out of order, adjacent in either direction), deduplicated with no cap, clipped to the
// block's piece; the edges of `has`, `any` and `first`; row rectangles joining a band; a section's marks as the overlay
// reads them, and the hover's primary as the first lit position the layout shows. Indices are block coordinate system
// indices unless a title says shown.

import { describe, expect, it } from "vitest";
import type { ShownIndex } from "../src/axis.js";
import { MarkBit } from "../src/index.js";

import { MarkSet, RowMarks, SectionMarks, type ShownRuns } from "../src/marks.js";
import { axisOf, rangedDoc, resolved } from "./axis-fixtures.js";

/** Block b's runs as pairs. */
const runsOf = (m: MarkSet, b: number) =>
  Array.from({ length: m.to(b) - m.from(b) }, (_, r) => [m.start[m.from(b) + r], m.stop[m.from(b) + r]]);
/** One block filled with `runs` in the order given, closed to [lo, hi). */
function filled(runs: number[][], lo = 0, hi = 1000) {
  const m = new MarkSet();
  m.reset(1);
  m.begin(0);
  for (const [s, t] of runs) m.add(s, t);
  m.close(lo, hi);
  return m;
}

describe("a MarkSet's runs", () => {
  it("in-order adds merge when they overlap or touch, and stay apart otherwise", () => {
    expect(
      runsOf(
        filled([
          [1, 3],
          [3, 5],
          [4, 6],
          [8, 9],
        ]),
        0,
      ),
    ).toEqual([
      [1, 6],
      [8, 9],
    ]);
  });

  it("descending adjacency merges into the last run: 5, then 4, then 3", () => {
    const m = filled([
      [5, 6],
      [4, 5],
      [3, 4],
    ]);
    expect(runsOf(m, 0)).toEqual([[3, 6]]);
    expect(m.count).toBe(1);
  });

  it("out-of-order adds are sorted at close, then merged", () => {
    expect(
      runsOf(
        filled([
          [50, 60],
          [10, 20],
          [30, 40],
          [15, 35],
          [70, 71],
          [0, 1],
        ]),
        0,
      ),
    ).toEqual([
      [0, 1],
      [10, 40],
      [50, 60],
      [70, 71],
    ]);
  });

  it("a merge that extends the last run back over an earlier one is merged at close", () => {
    expect(
      runsOf(
        filled([
          [0, 2],
          [5, 6],
          [1, 5],
        ]),
        0,
      ),
    ).toEqual([[0, 6]]);
  });

  it("duplicates collapse and nothing caps them: 100,000 positions in a shuffled order are one run", () => {
    const order = Array.from({ length: 100_000 }, (_, i) => (i * 7919) % 100_000);
    const m = new MarkSet();
    m.reset(1);
    m.begin(0);
    for (const i of order) m.add(i, i + 1);
    for (const i of order) m.add(i, i + 1);
    m.close(0, 1 << 30);
    expect(runsOf(m, 0)).toEqual([[0, 100_000]]);
  });

  it("10,000 disjoint positions in reverse are kept, sorted", () => {
    const m = new MarkSet();
    m.reset(1);
    m.begin(0);
    for (let i = 9999; i >= 0; i--) m.add(2 * i, 2 * i + 1);
    m.close(0, 1 << 30);
    expect(m.to(0) - m.from(0)).toBe(10_000);
    expect([m.start[0], m.start[9999]]).toEqual([0, 19_998]);
    for (let r = 1; r < 10_000; r++) expect(m.start[r]).toBeGreaterThan(m.stop[r - 1]);
  });

  it("clipping to the block's piece trims runs and drops those outside", () => {
    expect(
      runsOf(
        filled(
          [
            [0, 4],
            [6, 8],
            [9, 12],
            [20, 30],
          ],
          3,
          10,
        ),
        0,
      ),
    ).toEqual([
      [3, 4],
      [6, 8],
      [9, 10],
    ]);
    expect(runsOf(filled([[0, 3]], 3, 10), 0)).toEqual([]);
  });

  it("empty runs add nothing; blocks are filled in order, a skipped block is empty", () => {
    const m = new MarkSet();
    m.reset(4);
    m.begin(0);
    m.add(5, 5);
    m.close(0, 10);
    m.begin(2);
    m.add(1, 2);
    m.close(0, 10);
    expect([runsOf(m, 0), runsOf(m, 1), runsOf(m, 2), runsOf(m, 3)]).toEqual([[], [], [[1, 2]], []]);
    expect([m.first(1), m.first(2), m.first(3), m.first(-1)]).toEqual([-1, 1, -1, -1]);
  });

  it("reused: a reset keeps the buffers, the next fill sees none of the last", () => {
    const m = filled([
      [0, 1],
      [2, 3],
      [4, 5],
    ]);
    const buffers = [m.start, m.stop, m.at];
    m.reset(1);
    m.begin(0);
    m.add(7, 8);
    m.close(0, 10);
    expect(runsOf(m, 0)).toEqual([[7, 8]]);
    expect([m.start, m.stop, m.at]).toEqual(buffers);
  });

  it("carets: appended in any order, sorted, distinct and clipped per block by close (seqspec:interaction#channels)", () => {
    const m = new MarkSet();
    const carets = (b: number) => Array.from(m.caretAt.subarray(m.caretFrom(b), m.caretTo(b)));
    m.reset(4);
    m.begin(0);
    for (const at of [9, 3, 3, 12, 9, 0, 5, 5]) m.boundary(at);
    m.close(0, 10);
    // Block 1 skipped; block 2 in order, one past its range.
    m.begin(2);
    for (const at of [1, 2, 2, 4]) m.boundary(at);
    m.close(0, 3);
    expect([carets(0), carets(1), carets(2), carets(3)]).toEqual([[0, 3, 5, 9], [], [1, 2], []]);
    expect(Array.from(m.caretBlock.subarray(0, m.carets))).toEqual([0, 0, 0, 0, 2, 2]);
    m.reset(4);
    expect([m.carets, m.caretTo(0)]).toEqual([0, 0]);
  });
});

describe("has, any and first at their edges", () => {
  const m = filled([
    [2, 4],
    [7, 8],
    [10, 15],
  ]);

  it("has: a run's start in, its stop out, between runs out, in order and back", () => {
    const at = (i: number) => m.has(0, i);
    expect([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 99].map(at)).toEqual([
      false,
      false,
      true,
      true,
      false,
      false,
      false,
      true,
      false,
      false,
      true,
      true,
      false,
      false,
    ]);
    expect([at(14), at(2), at(7), at(3)]).toEqual([true, true, true, true]);
    expect([m.has(1, 3), m.has(-1, 3)]).toEqual([false, false]);
  });

  it("any: [s, t) overlapping a run, touching one, between two, empty", () => {
    expect(m.any(0, 0, 3)).toBe(true);
    expect(m.any(0, 0, 2)).toBe(false);
    expect(m.any(0, 4, 7)).toBe(false);
    expect(m.any(0, 4, 8)).toBe(true);
    expect(m.any(0, 8, 10)).toBe(false);
    expect(m.any(0, 14, 100)).toBe(true);
    expect(m.any(0, 15, 100)).toBe(false);
    expect(m.any(0, 3, 3)).toBe(false);
    expect(m.any(1, 0, 100)).toBe(false);
  });

  it("first: the lowest lit index", () => {
    expect(m.first(0)).toBe(2);
  });
});

describe("row rectangles", () => {
  it("consecutive rows with the same columns join into one band; a different row or columns starts another", () => {
    const r = new RowMarks();
    for (let row = 3; row < 9; row++) r.add(MarkBit.selection, 0, row, 10, 20);
    r.add(MarkBit.selection, 0, 9, 10, 21);
    r.add(MarkBit.highlight, 0, 10, 10, 21);
    r.add(MarkBit.highlight, 0, 12, 10, 21);
    expect(r.count).toBe(4);
    expect([r.rowStart[0], r.rowStop[0], r.start[0], r.stop[0]]).toEqual([3, 9, 10, 20]);
    r.clear(MarkBit.selection);
    expect([r.count, r.bit[0], r.rowStart[0], r.rowStart[1]]).toEqual([2, MarkBit.highlight, 10, 12]);
  });
});

describe("a section's marks as the overlay reads them", () => {
  const fill = (m: MarkSet, runs: number[][]) => {
    m.reset(2);
    m.begin(1);
    for (const [s, t] of runs) m.add(s, t);
    m.close(0, 100);
  };

  it("bitsAt ORs each channel's bit; any takes a mask; runs gives one channel's set", () => {
    const s = new SectionMarks();
    fill(s.hover, [[5, 6]]);
    fill(s.selection, [[0, 10]]);
    fill(s.accent, [[20, 30]]);
    for (const m of [s.highlight, s.band]) m.reset(2);
    expect(s.bitsAt(1, 5)).toBe(MarkBit.hover | MarkBit.selection);
    expect(s.bitsAt(1, 25)).toBe(MarkBit.accent);
    expect(s.bitsAt(0, 5)).toBe(0);
    expect(s.any(1, 11, 21, MarkBit.selection | MarkBit.hover)).toBe(false);
    expect(s.any(1, 11, 21, MarkBit.accent)).toBe(true);
    expect(s.runs(MarkBit.accent)).toBe(s.accent);
    expect(s.runs(MarkBit.band)).toBe(s.band);
  });

  it("the primary is the first lit position the layout shows; none shown, no tick", async () => {
    // [5, 15) filtered to [7, 10) and [12, 14) (view f): 5 and 6 are not shown.
    const ax = axisOf(await resolved(rangedDoc(), "f"));
    const s = new SectionMarks();
    s.hover.reset(1);
    s.hover.begin(0);
    s.hover.add(5, 8);
    s.hover.add(13, 14);
    s.hover.close(5, 15);
    s.updatePrimary(ax, 1);
    expect(s.primary(0)).toBe(7);
    s.hover.reset(1);
    s.hover.begin(0);
    s.hover.add(5, 7);
    s.hover.add(10, 12);
    s.hover.close(5, 15);
    s.updatePrimary(ax, 1);
    expect(s.primary(0)).toBe(-1);
    expect(s.primary(3)).toBe(-1);
  });

  it("on a reversed piece the primary is the lit position shown first", () => {
    // A stand-in axis over [0, 10) laid out in reverse: shown k is index 9 - k.
    const reversed: ShownRuns = {
      runToShown: (_b, start, stop, out) => {
        [out[0], out[1]] = [10 - stop, 10 - start];
        return 1;
      },
      indexOf: (_b, k: ShownIndex) => 9 - k,
    };
    const s = new SectionMarks();
    s.hover.reset(1);
    s.hover.begin(0);
    s.hover.add(2, 3);
    s.hover.add(6, 8);
    s.hover.close(0, 10);
    s.updatePrimary(reversed, 1);
    expect(s.primary(0)).toBe(7);
  });
});
