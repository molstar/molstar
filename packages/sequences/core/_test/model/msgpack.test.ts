/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// MessagePack (seqspec:format#values; msgpack/VENDORED.md): the typed-array extension decodes to identical numbers in
// every kind and payload length (fixext and ext forms), a decoder without it sees `{ type, data }`, and numbers JSON
// cannot hold follow JSON's rules.

import { describe, expect, it } from "vitest";
import {
  columnNumbers,
  decodeMsgPack,
  decodeTypedArrays,
  encodeMsgPack,
  encodeTypedArrays,
  float32Number,
} from "../../src/model/index.ts";
import { TypedArrayExt } from "../../src/model/msgpack/index.ts";

const roundTrip = (v: unknown) => decodeMsgPack(encodeMsgPack(v, encodeTypedArrays), decodeTypedArrays);
/** Every element the same number, -0 and NaN included. */
const same = (a: ArrayLike<number>, b: ArrayLike<number>) =>
  a.length === b.length && Array.from(a).every((x, i) => Object.is(x, b[i]));

describe("the typed-array extension", () => {
  const edges = {
    int8: Int8Array.of(-128, -1, 0, 127),
    uint8: Uint8Array.of(0, 1, 255),
    int16: Int16Array.of(-32768, 0, 32767),
    uint16: Uint16Array.of(0, 65535),
    int32: Int32Array.of(-(2 ** 31), -1, 0, 2 ** 31 - 1),
    uint32: Uint32Array.of(0, 2 ** 32 - 1),
    float32: Float32Array.of(0.15, -0, Number.NaN, 3.4028234663852886e38, 1e-45),
    float64: Float64Array.of(0.15, 0.3, 0.1 + 0.2, -0, Number.NaN, 1e-310, Number.MAX_VALUE, -Infinity, 2 ** 53 + 2),
  };

  it("each kind comes back as that kind with the same numbers", () => {
    for (const [name, a] of Object.entries(edges)) {
      const back = roundTrip({ a }) as { a: typeof a };
      expect(back.a.constructor, name).toBe(a.constructor);
      expect(same(back.a, a), name).toBe(true);
    }
    expect([...TypedArrayExt.keys()]).toEqual([0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18]);
  });

  it("every payload length: fixext 1, 2, 4, 8, 16 and ext 8, 16, 32; aligned when the payload is not", () => {
    for (const n of [0, 1, 2, 3, 4, 5, 31, 32, 8191, 8192, 10_000]) {
      const u8 = Uint8Array.from({ length: n }, (_, i) => i * 7);
      const f64 = Float64Array.from({ length: n }, (_, i) => i / 10);
      // A string first puts each payload at an odd offset.
      const back = roundTrip(["x", u8, "yy", f64]) as [string, Uint8Array, string, Float64Array];
      expect([same(back[1], u8), same(back[3], f64)], `n = ${n}`).toEqual([true, true]);
    }
    // The header forms by payload size: fixext 8 (one float64), fixext 16, ext 16, ext 32.
    const head = (n: number) => encodeMsgPack(new Float64Array(n), encodeTypedArrays)[0];
    expect([head(1), head(2), head(100), head(10_000)]).toEqual([0xd7, 0xd8, 0xc8, 0xc9]);
  });

  it("a Float32Array's numbers as its shortest decimals; a subarray only its elements", () => {
    const f32 = Float32Array.of(0.15, 1.5, 0.1 + 0.2, 16777217, 3.4028234663852886e38, Number.NaN);
    expect(columnNumbers(f32, null)).toEqual([0.15, 1.5, 0.3, 16777216, 3.4028235e38, null]);
    for (const x of Array.from({ length: 1000 }, (_, i) => Math.fround(Math.sin(i) * 10 ** ((i % 40) - 20)))) {
      expect(Math.fround(float32Number(x))).toBe(x);
    }
    const big = Float64Array.from({ length: 100 }, (_, i) => i);
    expect(Array.from(roundTrip(big.subarray(10, 13)) as Float64Array)).toEqual([10, 11, 12]);
  });

  it("a decoder without it sees { type, data }, a copy of the little-endian elements", () => {
    const plain = decodeMsgPack(encodeMsgPack(Int32Array.of(1, -2), encodeTypedArrays)) as {
      type: number;
      data: Uint8Array;
    };
    expect(plain.type).toBe(0x15);
    expect(Array.from(plain.data)).toEqual([1, 0, 0, 0, 0xfe, 0xff, 0xff, 0xff]);
    expect(() => decodeMsgPack(Uint8Array.of(0xc7, 3, 0x18, 1, 2, 3), decodeTypedArrays)).toThrow(/3 bytes/);
  });
});

describe("JSON's rules", () => {
  it("numbers: integers past 32 bits exactly, not finite as nil; undefined keys left out, undefined items nil", () => {
    const v = { big: 2 ** 40, neg: -(2 ** 40) - 1, max: Number.MAX_SAFE_INTEGER, nan: Number.NaN, gone: undefined };
    expect(roundTrip(v)).toEqual(JSON.parse(JSON.stringify(v)));
    expect(Object.keys(roundTrip(v) as object)).toEqual(["big", "neg", "max", "nan"]);
    expect(roundTrip([1, undefined, Number.POSITIVE_INFINITY])).toEqual([1, null, null]);
  });

  it("strings, nested maps and arrays, booleans and the int and float forms", () => {
    const v = { s: "Abé☃🧬", n: [0, -1, 127, 128, -33, 65536, -70000, 0.5], m: { t: true, f: false, z: null } };
    expect(roundTrip(v)).toEqual(v);
  });
});
