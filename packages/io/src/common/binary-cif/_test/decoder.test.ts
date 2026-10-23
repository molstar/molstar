/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ArrayEncoding as E } from '../array-encoder.js';
import { decode, decodeSteps, decodeStringDictionary } from '../decoder.js';
import { Encoding, type EncodedData } from '../encoding.js';

/** mulberry32: a seeded PRNG, so a failure reproduces */
function rng(seed: number) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const IntTypes = [
  [Encoding.IntDataType.Int8, Int8Array],
  [Encoding.IntDataType.Int16, Int16Array],
  [Encoding.IntDataType.Int32, Int32Array],
  [Encoding.IntDataType.Uint8, Uint8Array],
  [Encoding.IntDataType.Uint16, Uint16Array],
  [Encoding.IntDataType.Uint32, Uint32Array],
] as const;

/** Same type, length and bits (`Object.is`: -0 and NaN count). */
function expectIdentical(actual: ArrayLike<number>, expected: ArrayLike<number>) {
  expect(actual.constructor).toBe(expected.constructor);
  expect(actual.length).toBe(expected.length);
  for (let i = 0; i < expected.length; i++) {
    if (!Object.is(actual[i], expected[i])) expect([i, actual[i]]).toEqual([i, expected[i]]);
  }
}

/** A random chain of the fused form, with random (often inconsistent) sizes and values near the packing limits. */
function randomChain(r: () => number): EncodedData {
  const pick = <T>(xs: readonly T[]) => xs[Math.floor(r() * xs.length)];
  const n = pick([0, 1, 2, 3, 7, 50, 300]);
  const [type, ctor] = pick(IntTypes);
  const raw = new ctor(n);
  // the step decoder writes every repeat of a run, so run lengths stay small
  const hasRunLength = r() < 0.5;
  const limits = hasRunLength
    ? [0x7f, -0x80, 0xff]
    : [0x7f, -0x80, 0xff, 0x7fff, -0x8000, 0xffff, 0x7fffffff, -0x80000000, 0xffffffff];
  const scales = hasRunLength ? [4, 30] : [4, 300, 1e5, 1e10];
  for (let i = 0; i < n; i++) raw[i] = r() < 0.3 ? pick(limits) : Math.floor((r() - 0.4) * pick(scales));
  // a negative value in an unsigned array would be a run of billions
  if (hasRunLength && type >= Encoding.IntDataType.Uint8) for (let i = 0; i < n; i++) raw[i] &= 0xff;

  const chain: Encoding[] = [{ kind: 'ByteArray', type }];
  const size = () => pick([n, n + 1, Math.max(0, n - 1), 2 * n, Math.floor(n / 2)]);
  const intType = () => (r() < 0.85 ? Encoding.IntDataType.Int32 : pick(IntTypes)[0]);
  if (r() < 0.7) {
    chain.unshift({ kind: 'IntegerPacking', byteCount: pick([1, 2]), isUnsigned: r() < 0.5, srcSize: size() });
  }
  if (hasRunLength) chain.unshift({ kind: 'RunLength', srcType: intType(), srcSize: pick([size(), 3 * n]) });
  if (r() < 0.6) chain.unshift({ kind: 'Delta', origin: pick([0, 1, -5, 1e9]), srcType: intType() });
  const float = r() < 0.5 ? Encoding.FloatDataType.Float32 : Encoding.FloatDataType.Float64;
  const s = r();
  if (s < 0.35) chain.unshift({ kind: 'FixedPoint', factor: pick([1, 10, 1000, 3, -100]), srcType: float });
  else if (s < 0.55) chain.unshift({ kind: 'IntervalQuantization', min: -1.5, max: 7, numSteps: 13, srcType: float });

  // as msgpack gives it: a copy that owns its buffer
  return { encoding: chain, data: new Uint8Array(raw.buffer).slice() };
}

describe('binary-cif fused decode', () => {
  it('is bit-identical to the step-by-step decode on random chains', () => {
    const r = rng(1);
    for (let k = 0; k < 4000; k++) {
      const data = randomChain(r);
      let expected: ArrayLike<number> | undefined;
      try {
        expected = decodeSteps(data);
      } catch {
        // a chain the step decoder rejects (a negative size) must be rejected as well
        expect(() => decode(data)).toThrow();
        continue;
      }
      expectIdentical(decode(data), expected);
    }
  });

  it('maps StringArray indices as the step-by-step decode does', () => {
    const r = rng(3);
    const strings = ['ALA', 'GLY', 'C', 'CA', 'HOH'];
    const offsets = Int32Array.from({ length: strings.length + 1 }, (_, i) => strings.slice(0, i).join('').length);
    const offsetEncoded = E.by(E.delta).and(E.integerPacking).encode(offsets);
    for (let k = 0; k < 2000; k++) {
      const indices = randomChain(r);
      if (indices.encoding[0].kind === 'FixedPoint' || indices.encoding[0].kind === 'IntervalQuantization') continue;
      const data: EncodedData = {
        encoding: [
          {
            kind: 'StringArray',
            dataEncoding: indices.encoding,
            stringData: strings.join(''),
            offsetEncoding: offsetEncoded.encoding,
            offsets: offsetEncoded.data,
          },
        ],
        data: indices.data,
      };
      let expected: ArrayLike<number> | undefined;
      try {
        expected = decodeSteps(data);
      } catch {
        expect(() => decode(data)).toThrow();
        continue;
      }
      expectIdentical(decode(data), expected);
      // as a dictionary: the same strings, or none when an index is out of range
      const dictionary = decodeStringDictionary(data);
      if (Array.from(expected).some((v) => v === void 0)) {
        expect(dictionary).toBe(undefined);
      } else {
        expect(dictionary!.indices).toBeInstanceOf(Uint8Array);
        expect(Array.from(dictionary!.indices, (i) => dictionary!.values[i])).toEqual(Array.from(expected));
      }
    }
  });

  it('keeps dictionary indices in the smallest type that holds them', () => {
    for (const [count, type] of [
      [255, Uint8Array],
      [256, Uint16Array],
      [65535, Uint16Array],
      [70000, Int32Array],
    ] as const) {
      const strings = Array.from({ length: count }, (_, i) => `s${i}`);
      const offsets = new Int32Array(count + 1);
      for (let i = 0; i < count; i++) offsets[i + 1] = offsets[i] + strings[i].length;
      const offsetEncoded = E.by(E.delta).and(E.integerPacking).encode(offsets);
      const raw = Int32Array.from({ length: 10 }, (_, i) => (i * 7919) % count);
      const dataEncoded = E.by(E.byteArray).encode(raw);
      const data: EncodedData = {
        encoding: [
          {
            kind: 'StringArray',
            dataEncoding: dataEncoded.encoding,
            stringData: strings.join(''),
            offsetEncoding: offsetEncoded.encoding,
            offsets: offsetEncoded.data,
          },
        ],
        data: dataEncoded.data,
      };
      const dictionary = decodeStringDictionary(data)!;
      expect(dictionary.indices).toBeInstanceOf(type);
      expect(Array.from(dictionary.indices, (i) => dictionary.values[i])).toEqual(Array.from(raw, (i) => strings[i]));
    }
  });

  it('round-trips the encoders of real columns', () => {
    const r = rng(2);
    const coords = Float32Array.from({ length: 1000 }, () => Math.round((r() - 0.5) * 2e5) / 1000);
    const ids = Int32Array.from({ length: 1000 }, (_, i) => 1 + i + (i > 500 ? 40000 : 0));
    const seq = Int32Array.from({ length: 1000 }, (_, i) => 1 + Math.floor(i / 8));
    const cases: [ArrayLike<number>, ReturnType<typeof E.by>][] = [
      [coords, E.by(E.fixedPoint(1000)).and(E.delta).and(E.integerPacking)],
      [coords, E.by(E.fixedPoint(100)).and(E.runLength).and(E.integerPacking)],
      [ids, E.by(E.delta).and(E.runLength).and(E.integerPacking)],
      [seq, E.by(E.runLength).and(E.integerPacking)],
      [seq, E.by(E.delta).and(E.integerPacking)],
    ];
    for (const [values, encoder] of cases) {
      const encoded = encoder.encode(values);
      const decoded = decode(encoded);
      expectIdentical(decoded, decodeSteps(encoded));
      if (values instanceof Int32Array) expect(Array.from(decoded)).toEqual(Array.from(values));
    }
  });
});
