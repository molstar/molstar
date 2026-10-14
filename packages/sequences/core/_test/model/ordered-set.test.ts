/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from "vitest";
import { Interval, OrderedSet, SortedArray, SortedRanges } from "../../src/model/index.ts";

/** Seeded and deterministic (an LCG). */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

/** A random set as an OrderedSet and as a naive sorted list: contiguous, sparse, dense with holes, or empty. */
function randomSet(rnd: () => number): [OrderedSet, number[]] {
  const shape = Math.floor(rnd() * 4);
  const lo = Math.floor(rnd() * 300);
  if (shape === 0) {
    const stop = lo + Math.floor(rnd() * 200);
    return [OrderedSet.ofRange(lo, stop), Array.from({ length: stop - lo }, (_, i) => lo + i)];
  }
  if (shape === 3 && rnd() < 0.3) return [OrderedSet.Empty, []];
  const n = Math.floor(rnd() * (shape === 1 ? 20 : 400));
  const p = shape === 1 ? 0.05 : 0.9;
  const xs: number[] = [];
  for (let x = lo; xs.length < n; x++) if (rnd() < p) xs.push(x);
  return [OrderedSet.ofSortedArray(Int32Array.from(xs)), xs];
}

const toList = (s: OrderedSet) => {
  const out: number[] = [];
  OrderedSet.forEach(
    s,
    (x, _i, ctx) => {
      ctx.push(x);
    },
    out,
  );
  return out;
};

describe("OrderedSet against a naive model", () => {
  it("answers membership, positions and lower bounds on random sets", () => {
    const rnd = random(42);
    // Disagreements collected and compared once: an expect per query was most of the run.
    const wrong: string[] = [];
    for (let trial = 0; trial < 300; trial++) {
      const [s, xs] = randomSet(rnd);
      expect(OrderedSet.size(s)).toBe(xs.length);
      expect(toList(s)).toEqual(xs);
      for (let k = 0; k < 40; k++) {
        const q = Math.floor(rnd() * 800) - 50;
        const i = xs.indexOf(q);
        const lb = xs.findIndex((x) => x >= q);
        const got = [OrderedSet.has(s, q), OrderedSet.indexOf(s, q), OrderedSet.lowerBound(s, q)];
        const want = [i >= 0, i, lb < 0 ? xs.length : lb];
        if (got.join() !== want.join()) wrong.push(`#${trial} ${q}: has, indexOf, lowerBound ${got} (want ${want})`);
      }
      for (let i = 0; i < xs.length; i += 7)
        if (OrderedSet.getAt(s, i) !== xs[i])
          wrong.push(`#${trial} getAt ${i}: ${OrderedSet.getAt(s, i)} (want ${xs[i]})`);
    }
    expect(wrong).toEqual([]);
  });

  it("unions and intersects as sets do, in both forms", () => {
    const rnd = random(7);
    for (let trial = 0; trial < 500; trial++) {
      const [a, xa] = randomSet(rnd);
      const [b, xb] = randomSet(rnd);
      const union = [...new Set([...xa, ...xb])].sort((x, y) => x - y);
      const inter = xa.filter((x) => xb.includes(x));
      const u = OrderedSet.union(a, b);
      const n = OrderedSet.intersect(a, b);
      expect(toList(u)).toEqual(union);
      expect(toList(n)).toEqual(inter);
      // Contiguous results are intervals.
      const contiguous = (xs: number[]) => xs.length > 0 && xs[xs.length - 1] - xs[0] + 1 === xs.length;
      if (contiguous(union)) expect(Interval.is(u)).toBe(true);
      if (contiguous(inter)) expect(Interval.is(n)).toBe(true);
      expect(OrderedSet.areEqual(u, OrderedSet.ofSortedArray(Int32Array.from(union)))).toBe(true);
      expect(OrderedSet.areEqual(a, b)).toBe(xa.join() === xb.join());
    }
  });

  it("makes an interval of a contiguous array and keeps a sparse one", () => {
    expect(Interval.is(OrderedSet.ofSortedArray(Int32Array.of(3, 4, 5)))).toBe(true);
    expect(OrderedSet.ofSortedArray(Int32Array.of(3, 5))).toEqual(Int32Array.of(3, 5));
    expect(OrderedSet.size(OrderedSet.ofRange(5, 2))).toBe(0);
    expect(SortedArray.ofUnsortedArray(Int32Array.of(5, 1, 5, 3))).toEqual(Int32Array.of(1, 3, 5));
  });
});

describe("the empty set", () => {
  // Every way of making an empty set is the one `Empty`, so equality is transitive: an empty ofRange(5,5) used to
  // differ from Empty (raw words) while both equalled an empty sorted array.
  it("is one value however it was made", () => {
    const empties: OrderedSet[] = [
      OrderedSet.Empty,
      Interval.Empty,
      OrderedSet.ofRange(5, 5),
      OrderedSet.ofRange(7, 3),
      Interval.ofRange(-4, -4),
      OrderedSet.ofSortedArray(new Int32Array(0)),
      OrderedSet.intersect(OrderedSet.ofRange(0, 3), OrderedSet.ofRange(5, 8)),
      OrderedSet.intersect(OrderedSet.ofRange(0, 3), Int32Array.of(7, 9)),
      OrderedSet.union(OrderedSet.ofRange(2, 2), OrderedSet.ofRange(9, 9)),
      new Int32Array(0),
    ];
    for (const a of empties) {
      expect(OrderedSet.size(a)).toBe(0);
      for (const b of empties) expect(OrderedSet.areEqual(a, b)).toBe(true);
    }
    for (const a of empties.filter((e) => typeof e === "number")) expect(a).toBe(OrderedSet.Empty);
    expect(Interval.areEqual(Interval.ofRange(5, 5), Interval.ofRange(9, 1))).toBe(true);
  });

  it("is not equal to a set with a member", () => {
    expect(OrderedSet.areEqual(OrderedSet.ofRange(5, 5), OrderedSet.ofRange(5, 6))).toBe(false);
    expect(OrderedSet.areEqual(OrderedSet.Empty, Int32Array.of(0))).toBe(false);
  });
});

describe("SortedRanges", () => {
  it("answers membership over half-open pairs", () => {
    const r = Int32Array.of(2, 5, 8, 9, 20, 30);
    const naive = (x: number) => (x >= 2 && x < 5) || x === 8 || (x >= 20 && x < 30);
    for (let x = -1; x < 32; x++) expect(SortedRanges.has(r, x)).toBe(naive(x));
    expect(SortedRanges.size(r)).toBe(14);
    expect(SortedRanges.count(r)).toBe(3);
  });
});
