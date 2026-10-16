/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  hashMurmur128o,
  murmurHash3_32,
  murmurHash3_32_fromBytes,
  murmurHash3_128,
  murmurHash3_128_fromBytes,
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
