/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { utf8ByteCount, utf8Read, utf8Write } from '../utf8.js';

function write(s: string) {
  const data = new Uint8Array(utf8ByteCount(s));
  utf8Write(data, 0, s);
  return data;
}

describe('utf8', () => {
  it('writes ASCII and BMP text as TextEncoder does', () => {
    for (const s of ['', 'abc', 'é', '☃', 'Abé☃ ࠀ￿']) {
      expect(Array.from(write(s)), s).toEqual(Array.from(new TextEncoder().encode(s)));
      expect(utf8ByteCount(s)).toBe(new TextEncoder().encode(s).length);
    }
  });

  it('writes a surrogate pair as one 4-byte sequence', () => {
    for (const s of ['🧬', 'a🧬b😀', '\u{10ffff}']) {
      expect(utf8ByteCount(s), s).toBe(new TextEncoder().encode(s).length);
      expect(Array.from(write(s)), s).toEqual(Array.from(new TextEncoder().encode(s)));
      expect(utf8Read(write(s), 0, utf8ByteCount(s))).toBe(s);
    }
  });

  it('writes a lone surrogate as 3 bytes', () => {
    expect(utf8ByteCount('\ud800')).toBe(3);
  });
});
