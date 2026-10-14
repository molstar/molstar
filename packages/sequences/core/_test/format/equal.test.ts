/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `equal`: the canonical JSON's equality, compared in place (the equality sweep's one helper).
import { describe, expect, it } from "vitest";
import { canonical, equal } from "../../src/format/index.js";
import { itWithin } from "../model/timing.js";

/** A seeded generator (mulberry32), so a failure repeats. */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (s + 0x6d2b79f5) >>> 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const NUMBERS = [0, -0, 1, 0.1, 1.5, Math.fround(0.1), 2 ** 31, Number.NaN, Number.POSITIVE_INFINITY];

/** A small JSON-like value, typed arrays, NaN, undefined and shared objects included. */
function value(r: () => number, depth: number): unknown {
  const k = Math.floor(r() * (depth > 2 ? 6 : 10));
  const pick = <T>(xs: readonly T[]) => xs[Math.floor(r() * xs.length)] as T;
  if (k === 0) return pick([null, undefined]);
  if (k === 1) return pick(NUMBERS);
  if (k === 2) return pick(["", "a", "b", "0"]);
  if (k === 3) return r() < 0.5;
  if (k === 4)
    return new (pick([Float32Array, Float64Array, Int32Array]))(
      Array.from({ length: Math.floor(r() * 3) }, () => pick(NUMBERS)),
    );
  if (k === 5) return Array.from({ length: Math.floor(r() * 3) }, () => pick(NUMBERS));
  if (k < 8) return Array.from({ length: Math.floor(r() * 3) }, () => value(r, depth + 1));
  const o: Record<string, unknown> = {};
  for (const key of ["x", "y", "z"]) if (r() < 0.5) o[key] = value(r, depth + 1);
  return o;
}

/** The same value rebuilt: other objects, keys in another order, typed arrays as lists of their JSON numbers. */
function rebuilt(x: unknown, r: () => number): unknown {
  if (x === null || typeof x !== "object") return x;
  if (ArrayBuffer.isView(x)) return r() < 0.5 ? x : JSON.parse(canonical(x));
  if (Array.isArray(x)) return x.map((v) => rebuilt(v, r));
  const keys = Object.keys(x).reverse();
  return Object.fromEntries(keys.map((k) => [k, rebuilt((x as Record<string, unknown>)[k], r)]));
}

describe("equal", () => {
  it("agrees with canonical text on 20,000 random pairs, equal and not", () => {
    const r = random(7);
    let same = 0;
    const wrong: string[] = [];
    for (let i = 0; i < 20000; i++) {
      const a = value(r, 0);
      const b = r() < 0.5 ? rebuilt(a, r) : value(r, 0);
      const want = canonical(a) === canonical(b);
      if (want) same++;
      if (equal(a, b) !== want || equal(b, a) !== want) wrong.push(`${canonical(a)} vs ${canonical(b)}: not ${want}`);
    }
    expect(wrong).toEqual([]);
    expect(same).toBeGreaterThan(5000);
  });

  it("keeps canonical's equivalences", () => {
    expect(equal(new Float32Array([0.1, Number.NaN]), [0.1, null])).toBe(true);
    expect(equal(new Float32Array([0.1]), [Math.fround(0.1)])).toBe(false);
    expect(equal(new Float64Array([1, 2]), new Int32Array([1, 2]))).toBe(true);
    expect(equal(Number.NaN, null)).toBe(true);
    expect(equal(undefined, Number.POSITIVE_INFINITY)).toBe(true);
    expect(equal(0, -0)).toBe(true);
    expect(equal({ a: 1, b: [2] }, { b: [2], a: 1 })).toBe(true);
    expect(equal({ a: undefined }, {})).toBe(false);
    expect(equal([1], { 0: 1 })).toBe(false);
    expect(equal("1", 1)).toBe(false);
  });

  itWithin("compares a 9M-cell grid in place, sharing nothing, in a few tens of ms", (within) => {
    const a = new Float32Array(3000 * 3000).map((_, i) => i % 97);
    const b = a.slice();
    within(500, () => {
      expect(equal({ grid: a }, { grid: b })).toBe(true);
      b[b.length - 1] = -1;
      expect(equal({ grid: a }, { grid: b })).toBe(false);
    });
  });
});
