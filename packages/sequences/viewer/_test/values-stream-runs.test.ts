/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// VT-06's draw side: a streamed values column read run by run. A device column across a cyclic block's origin (or a
// filter's gap) asks the stream each run's positions, not the hull of its two ends (the whole system but the column).

import { expect, it } from "vitest";
import { streamColumn } from "../src/draw/api.js";

it("a column across a cyclic origin asks its two runs, and merges their statistics", () => {
  const L = 5_000_000;
  // Shown k from the origin at 2,500,000: block index (k + L − 2.5 M) mod L.
  const indexOf = (_b: number, k: number) => (k + L - 2_500_000) % L;
  const asked: [number, number][] = [];
  const v = {
    axis: { indexOf },
    cell: 1 / 5000,
    span: [1],
    dpr: 1,
    shownToRuns: (_b: number, k0: number, k1: number, out: Int32Array) => {
      // [k0, 2.5 M) and [2.5 M, k1) in shown order: [L − 2,500, L) then [0, 2,500).
      out.set([indexOf(0, k0), L, 0, indexOf(0, k1 - 1) + 1]);
      return 2;
    },
    stream: {
      values: (_b: number, i0: number, i1: number, out: Float64Array) => {
        asked.push([i0, i1]);
        // Each run its own range of values: the left one 1..3, the right one 5..9.
        const left = i0 > 0;
        out[0] = left ? 1 : 5;
        out[1] = left ? 3 : 9;
        out[2] = left ? 20 : 70;
        out[3] = 10;
        return true;
      },
    },
  };
  const out = new Float64Array(5);
  streamColumn(v as never, 0, 2_497_500, 2_502_500, out);
  expect(asked).toEqual([
    [L - 2_500, L],
    [0, 2_500],
  ]);
  expect(Array.from(out.subarray(0, 4))).toEqual([1, 9, 90, 20]);
});
