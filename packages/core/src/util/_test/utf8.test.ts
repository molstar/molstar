/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { asciiIndexOf, asciiStartsWith, utf8ByteCount, utf8Read, utf8ReadValue, utf8Write } from '../utf8.js';

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

  it('reads short ASCII as TextDecoder does, and falls back for the rest', () => {
    const decoder = new TextDecoder();
    for (const s of [
      '',
      'a',
      'CA',
      'HOH',
      'ATOM',
      "C1'",
      'abcdefgh',
      'pdbx_PDB_model_num',
      'x'.repeat(24),
      'x'.repeat(25),
    ]) {
      const data = new TextEncoder().encode(` ${s} `);
      expect(utf8Read(data, 1, data.length - 2)).toBe(s);
    }
    for (const s of ['é', 'aé', '☃b', '🧬', '\ufeffdata_']) {
      const data = new TextEncoder().encode(s);
      expect(utf8Read(data, 0, data.length)).toBe(decoder.decode(data));
    }
    // a range past the end reads what is there, as the decoder of a subarray does
    const data = new TextEncoder().encode('abc');
    expect(utf8Read(data, 1, 5)).toBe('bc');
  });

  it('reads a value with a leading U+FEFF as it is', () => {
    const bom = String.fromCharCode(0xfeff);
    const data = new TextEncoder().encode(`${bom}abc`);
    expect(utf8Read(data)).toBe('abc');
    expect(utf8ReadValue(data)).toBe(`${bom}abc`);
    expect(utf8ReadValue(new TextEncoder().encode(' é '), 1, 2)).toBe('é');
  });
});

describe('asciiIndexOf / asciiStartsWith', () => {
  const texts = ['', 'a', 'data_x\n_atom_site.id 1\n', 'ééé _refln.pdbx_FWT ü', '__a_ab_abc', 'abcabcabd'];
  const needles = ['', 'a', '_', '_refln.pdbx_FWT', 'abd', 'abc', 'x\n_', 'zz'];

  it('equals String.indexOf on ASCII text, as byte indices', () => {
    for (const t of texts) {
      const bytes = write(t);
      for (const n of needles) {
        for (const from of [0, 1, 3]) {
          const i = t.indexOf(n, from);
          const expected = i < 0 ? -1 : utf8ByteCount(t.substring(0, i));
          if (!/[^\x00-\x7f]/.test(t)) expect(asciiIndexOf(bytes, n, from)).toBe(expected);
          else if (from === 0) expect(asciiIndexOf(bytes, n) >= 0).toBe(i >= 0);
        }
        expect(asciiStartsWith(bytes, n)).toBe(t.startsWith(n));
      }
    }
  });
});
