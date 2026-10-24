/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Sniffable } from '../sniffable.js';

describe('Sniffable', () => {
  const text = 'data_x\n_atom_site.id 1\nééé _refln.pdbx_FWT\n';
  const forms = { string: Sniffable.of(text), bytes: Sniffable.of(new TextEncoder().encode(text)) };

  it('searches a string and its bytes alike', () => {
    for (const s of Object.values(forms)) {
      expect(s.sniffText('_refln.pdbx_FWT')).toBe(true);
      expect(s.sniffText('_volume_data_3d_info.')).toBe(false);
      expect(s.sniffText('data_', 'start')).toBe(true);
      expect(s.sniffText('_atom', 'start')).toBe(false);
      expect(s.sniffBytes('data_')).toBe(true);
      expect(s.sniffBytes([0x5f, 0x61], 7)).toBe(true);
      expect(Array.from(s.header.subarray(0, 5))).toEqual([100, 97, 116, 97, 95]);
    }
    expect(forms.bytes.heldBytes).toBeInstanceOf(Uint8Array);
    expect(forms.string.heldBytes).toBeUndefined();
  });

  it('keeps a header of at most HeaderBytes', () => {
    const big = new Uint8Array(Sniffable.HeaderBytes + 10);
    expect(Sniffable.of(big).header.length).toBe(Sniffable.HeaderBytes);
    expect(Sniffable.of('x'.repeat(Sniffable.HeaderBytes + 10)).header.length).toBe(Sniffable.HeaderBytes);
  });

  it('is told by its sniffs, and `of` gives back a Sniffable as it is', () => {
    expect(Sniffable.is(forms.bytes)).toBe(true);
    expect(Sniffable.is(new Uint8Array(1))).toBe(false);
    expect(Sniffable.is('x')).toBe(false);
    expect(Sniffable.of(forms.bytes)).toBe(forms.bytes);
  });
});
