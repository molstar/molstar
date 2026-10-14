/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `float32Number` (common.ts): the shortest decimal that reads back as a float32, from 7 digits first, against the
// search from 1 digit up on a million random float32s, every float32 in a few dense ranges and the edge values (the slow
// tier; the fast tier a sample of each). Every float32 was checked once the same way.

import { describe, expect, it } from "vitest";
import { float32Number } from "../../src/model/index.js";

/** The search from 1 digit up: the nearest decimal of 1, 2, … 9 digits that reads back. */
function slow(x: number): number {
  for (let p = 1; ; p++) {
    const y = Number(x.toPrecision(p));
    if (p === 9 || Math.fround(y) === x) return y;
  }
}

/** The three powers of two where an 8-digit decimal other than the nearest reads back (the search gives 9 digits). */
const POWERS: [number, number][] = [
  [2 ** -96, 1.2621775e-29],
  [2 ** 87, 1.5474251e26],
  [2 ** 90, 1.2379401e27],
];
const expected = (x: number) => POWERS.find(([p]) => p === Math.abs(x))?.[1] ?? slow(x);

const f = new Float32Array(1);
const bits = new Uint32Array(f.buffer);
/** The float32s whose bit patterns are `from` to `to` (inclusive) that differ from the expected; NaN skipped. */
function differences(from: number, to: number): string[] {
  const out: string[] = [];
  for (let b = from; b <= to; b++) {
    bits[0] = b;
    const x = f[0];
    if (!Number.isNaN(x) && float32Number(x) !== Math.sign(x) * expected(Math.abs(x)))
      out.push(`${b.toString(16)} ${x}`);
  }
  return out;
}

/** `count` xorshift float32s from `seed` that differ from the expected. */
function randomDifferences(count: number, seed = 2463534242): string[] {
  let s = seed;
  const out: string[] = [];
  for (let i = 0; i < count; i++) {
    s ^= s << 13;
    s ^= s >>> 17;
    s ^= s << 5;
    out.push(...differences(s >>> 0, s >>> 0));
  }
  return out;
}

const at = (x: number) => {
  f[0] = x;
  return bits[0];
};

describe("float32Number", () => {
  it("equals the search from 1 digit up on 50,000 random float32s", () => {
    expect(randomDifferences(50_000)).toEqual([]);
  });

  it("equals the search from 1 digit up on a million random float32s", { tags: ["slow"] }, () => {
    expect(randomDifferences(1_000_000)).toEqual([]);
  });

  it("and on float32s in dense ranges, a sample: subnormals, around 1, past 2^43, about every power of two", () => {
    expect(differences(0, 10_000)).toEqual([]);
    expect(differences(at(0.9999), at(1.0001))).toEqual([]);
    expect(differences(at(2 ** 43), at(2 ** 43) + 20_000)).toEqual([]);
    expect(differences(at(1e13) - 20_000, at(1e13))).toEqual([]);
    for (let e = -126; e <= 127; e++) expect(differences(at(2 ** e) - 20, at(2 ** e) + 20)).toEqual([]);
  });

  it("and on every float32 in dense ranges: subnormals, around 1, and where float32 steps exceed 7 digits", {
    tags: ["slow"],
  }, () => {
    expect(differences(0, 100_000)).toEqual([]);
    expect(differences(at(0.999), at(1.001))).toEqual([]);
    // [2^43, 1e13): float32's step (2^20) exceeds the 7-digit one (10^6), so a shorter decimal may read back too.
    expect(differences(at(2 ** 43), at(1e13))).toEqual([]);
    for (let e = -126; e <= 127; e++) expect(differences(at(2 ** e) - 200, at(2 ** e) + 200)).toEqual([]);
  });

  it("edge values: zeros, not finite, subnormals, the largest, integers, the three powers of two", () => {
    const cases: [number, number][] = [
      [0, 0],
      [-0, 0],
      [Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY],
      [Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY],
      [Math.fround(1e-45), 1e-45],
      [Math.fround(-1.17549435e-38), -1.1754944e-38],
      [Math.fround(3.4028234663852886e38), 3.4028235e38],
      [16777216, 16777216],
      [Math.fround(16777217), 16777216],
      [Math.fround(0.15), 0.15],
      [Math.fround(1 / 3), 0.33333334],
      [-1, -1],
      [123456, 123456],
      ...POWERS,
      [-(2 ** 87), -1.5474251e26],
    ];
    expect(cases.map(([x]) => float32Number(x))).toEqual(cases.map(([, y]) => y));
    expect(Object.is(float32Number(-0), 0)).toBe(true);
    expect(float32Number(Number.NaN)).toBeNaN();
    for (const [x, y] of POWERS)
      expect([Math.fround(y) === x, String(y).length < String(slow(x)).length]).toEqual([true, true]);
  });
});
