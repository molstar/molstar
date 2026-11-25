/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { itWithin, ratioWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import {
  hashMurmur128o,
  murmurHash3_32,
  murmurHash3_32_fromBytes,
  murmurHash3_128,
  murmurHash3_128_fromBytes,
  PCG,
} from '../hash-functions.js';

/** The 128-bit hex (h1 to h4, each as 8 hex digits) as the reference implementation's 16 output bytes. */
function bytesOf128(hex: string): Uint8Array {
  const out = new Uint8Array(16);
  const view = new DataView(out.buffer);
  for (let i = 0; i < 4; i++) view.setUint32(4 * i, parseInt(hex.slice(8 * i, 8 * i + 8), 16), true);
  return out;
}

/**
 * SMHasher's verification value: hash keys {}, {0}, {0, 1}, …, {0, …, 254} with seed 256 - length, hash the
 * concatenated outputs with seed 0, and read the first 4 bytes little-endian.
 */
function verification(width: 4 | 16, hash: (key: Uint8Array, seed: number) => Uint8Array): number {
  const key = new Uint8Array(256);
  const hashes = new Uint8Array(256 * width);
  for (let i = 0; i < 256; i++) {
    key[i] = i;
    hashes.set(hash(key.subarray(0, i), 256 - i), i * width);
  }
  return new DataView(hash(hashes, 0).buffer).getUint32(0, true);
}

describe('MurmurHash3_x86_32', () => {
  it("matches SMHasher's verification value", () => {
    const hash = (key: Uint8Array, seed: number) => {
      const out = new Uint8Array(4);
      new DataView(out.buffer).setUint32(0, murmurHash3_32_fromBytes(key, seed), true);
      return out;
    };
    expect(verification(4, hash)).toBe(0xb0f57ee3);
  });

  it('matches published vectors', () => {
    expect(murmurHash3_32('', 0)).toBe(0);
    expect(murmurHash3_32('', 1)).toBe(0x514e28b7);
    expect(murmurHash3_32('', 0xffffffff)).toBe(0x81f16f39);
    expect(murmurHash3_32('\0\0\0\0', 0)).toBe(0x2362f9de);
    expect(murmurHash3_32('aaaa', 0x9747b28c)).toBe(0x5a97808a);
    expect(murmurHash3_32('Hello, world!', 0x9747b28c)).toBe(0x24884cba);
    expect(murmurHash3_32('The quick brown fox jumps over the lazy dog', 0x9747b28c)).toBe(0x2fa826cd);
  });
});

describe('MurmurHash3_x86_128', () => {
  it("matches SMHasher's verification value", () => {
    expect(verification(16, (key, seed) => bytesOf128(murmurHash3_128_fromBytes(key, seed)))).toBe(0xb3ece62a);
  });

  // From an independent implementation of the reference (smhasher's MurmurHash3.cpp), over the UTF-8 bytes.
  it('matches reference vectors', () => {
    expect(murmurHash3_128('', 42)).toBe('af6d2cb695c80cba95c80cba95c80cba');
    expect(murmurHash3_128('a', 42)).toBe('517b4f5294460fc694460fc694460fc6');
    expect(murmurHash3_128('Hello, world!', 42)).toBe('205186c50f5fb050233ccf1317f38f48');
    expect(murmurHash3_128('The quick brown fox jumps over the lazy dog', 42)).toBe('b0c69c19b1fd95c704c746bdb64fcfec');
    expect(murmurHash3_128('{"label":"β-sheet 🧬"}', 42)).toBe('cf79ad150d8bbd7a28a36457a30dca03');
  });

  it('every tail length agrees between the string and the bytes form', () => {
    const encoder = new TextEncoder();
    for (let n = 0; n <= 40; n++) {
      const text = 'xyzŁ€🧬'.repeat(8).slice(0, n);
      expect(murmurHash3_128(text, 7)).toBe(murmurHash3_128_fromBytes(encoder.encode(text), 7));
      expect(murmurHash3_32(text, 7)).toBe(murmurHash3_32_fromBytes(encoder.encode(text), 7));
    }
  });
});

describe('strings are hashed as UTF-8', () => {
  it('characters above U+00FF do not collide with their low byte', () => {
    // 'ł' is U+0142, whose low byte is 'B': hashing code units' low bytes made these equal.
    expect(murmurHash3_128('Ałpha helix, 16c', 42)).toBe('ed186557d3d9c6f0944f4263a5b4c9ea');
    expect(murmurHash3_128('ABpha helix, 16c', 42)).toBe('23349d609686e845aa7abc00ebd13af1');
    expect(murmurHash3_32('Ał', 0)).not.toBe(murmurHash3_32('AB', 0));
    expect(hashMurmur128o({ label: 'ł' })).not.toBe(hashMurmur128o({ label: 'B' }));
  });

  it('a long string grows the scratch buffer and still hashes right', () => {
    const long = 'ł'.repeat(5000);
    expect(murmurHash3_128(long, 1)).toBe(murmurHash3_128_fromBytes(new TextEncoder().encode(long), 1));
    expect(murmurHash3_128('short', 1)).toBe(murmurHash3_128_fromBytes(new TextEncoder().encode('short'), 1));
  });
});

/** The PCG before PCG32 (32 bits of state, the 64-bit output step): the speed baseline only. */
class OldPCG {
  private state: number;
  constructor(seed = 26699) {
    this.state = seed >>> 0;
  }
  int(): number {
    const oldstate = this.state;
    this.state = (Math.imul(this.state, 1664525) + 1013904223) >>> 0;
    const xorshifted = ((oldstate >>> 18) ^ oldstate) >>> 5;
    const rot = oldstate >>> 27;
    return ((xorshifted >>> rot) | (xorshifted << (32 - rot))) >>> 0;
  }
}

describe('PCG', () => {
  const Draws = 1_000_000;

  it("matches pcg32's reference outputs (pcg32_srandom(seed, 54))", () => {
    // pcg32-demo's first round for seed 42, stream 54 (https://www.pcg-random.org/, pcg-c-basic)
    const draw = (seed: number) => {
      const rng = new PCG(seed);
      return Array.from({ length: 6 }, () => rng.int());
    };
    expect(draw(42)).toEqual([0xa15c02b7, 0x7b47f409, 0xba1d3330, 0x83d2f293, 0xbfa4784b, 0xcbed606e]);
    // the same reference (64-bit arithmetic) for the default seed, 0 and the largest seed
    expect(draw(26699)).toEqual([0x3ac5b7c3, 0x9d02ce30, 0xa488c533, 0x2300c4fa, 0x4a00946f, 0x805f2ce4]);
    expect(draw(0)).toEqual([0x47c28b93, 0xb98f6a27, 0x7d3dcb1e, 0xf0761116, 0x9cc33f5b, 0xbe0e744d]);
    expect(draw(0xffffffff)).toEqual([0x1836f28a, 0x41720992, 0x67039735, 0xe39c0b28, 0xdd934c44, 0x04b356a9]);
    const first = new PCG().int();
    expect(first).toBe(0x3ac5b7c3);
  });

  it('equals the published algorithm in 64-bit arithmetic over many seeds and long runs', () => {
    // pcg32_srandom_r and pcg32_random_r of pcg-c-basic, written with BigInt as an independent oracle: the two-halves
    // arithmetic's carries are exercised deep into each sequence, not only in the first draws
    const Mask = (1n << 64n) - 1n;
    const Mul = 6364136223846793005n;
    const oracle = (seed: number) => {
      const inc = (54n << 1n) | 1n;
      let state = 0n;
      const step = () => {
        state = (state * Mul + inc) & Mask;
      };
      step();
      state = (state + BigInt(seed >>> 0)) & Mask;
      step();
      return () => {
        const old = state;
        step();
        const xorshifted = Number((((old >> 18n) ^ old) >> 27n) & 0xffffffffn);
        const rot = Number(old >> 59n);
        return ((xorshifted >>> rot) | (xorshifted << (-rot & 31))) >>> 0;
      };
    };
    for (const seed of [0, 1, 7, 42, 26699, 0x7fffffff, 0x80000000, 0xdeadbeef, 0xffffffff]) {
      const rng = new PCG(seed);
      const ref = oracle(seed);
      for (let i = 0; i < 2000; i++) {
        const a = rng.int();
        const b = ref();
        if (a !== b) expect({ seed, draw: i, value: a }).toEqual({ seed, draw: i, value: b });
      }
    }
    // the carry from the low half into the high half after adding the increment (about one step in 4·10^7, so no
    // sequence above reaches it): a state whose low half times the multiplier ends within the increment of 2^32
    const MulLo = Mul & 0xffffffffn;
    let inv = 1n; // MulLo⁻¹ mod 2^32 by Newton's iteration
    for (let i = 0; i < 5; i++) inv = (inv * (2n - MulLo * inv)) & 0xffffffffn;
    for (const hi of [0n, 0x12345678n, 0xffffffffn]) {
      const lo = (0xffffffffn * inv) & 0xffffffffn; // lo · MulLo ≡ 2^32 − 1, so adding the increment carries
      const start = (hi << 32n) | lo;
      const rng = new PCG(0) as unknown as { hi: number; lo: number; int(): number };
      rng.hi = Number(hi);
      rng.lo = Number(lo);
      let state = start;
      const inc = (54n << 1n) | 1n;
      for (let i = 0; i < 3; i++) {
        const old = state;
        state = (state * Mul + inc) & Mask;
        const xorshifted = Number((((old >> 18n) ^ old) >> 27n) & 0xffffffffn);
        const rot = Number(old >> 59n);
        expect(rng.int()).toBe(((xorshifted >>> rot) | (xorshifted << (-rot & 31))) >>> 0);
      }
    }
  });

  it('floats fall into eighths evenly, within 1% over 1M draws', () => {
    const rng = new PCG(1);
    const counts = new Float64Array(8);
    let min = 1;
    let max = 0;
    for (let i = 0; i < Draws; i++) {
      const f = rng.float();
      if (f < min) min = f;
      if (f > max) max = f;
      counts[Math.floor(f * 8)]++;
    }
    expect(min).toBeGreaterThanOrEqual(0);
    expect(max).toBeLessThan(1);
    for (const c of counts) expect(Math.abs(c - Draws / 8) / (Draws / 8)).toBeLessThan(0.01);
  });

  it('consecutive draws are uncorrelated (|r| < 0.01 over 1M pairs)', () => {
    const rng = new PCG(2);
    let prev = rng.float();
    let sx = 0;
    let sy = 0;
    let sxx = 0;
    let syy = 0;
    let sxy = 0;
    for (let i = 0; i < Draws; i++) {
      const next = rng.float();
      sx += prev;
      sy += next;
      sxx += prev * prev;
      syy += next * next;
      sxy += prev * next;
      prev = next;
    }
    const n = Draws;
    const r = (n * sxy - sx * sy) / Math.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
    expect(Math.abs(r)).toBeLessThan(0.01);
  });

  itWithin('draws no slower than twice the old generator', (within) => {
    const run = (rng: { int(): number }) => () => {
      let acc = 0;
      for (let i = 0; i < Draws; i++) acc ^= rng.int();
      return acc;
    };
    ratioWithin(within, 2, run(new OldPCG(3)), run(new PCG(3)));
  });
});
