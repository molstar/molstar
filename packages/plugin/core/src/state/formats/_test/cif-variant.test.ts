/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { getFileNameInfo } from '@molstar/core/util/file-info';
import { Sniffable } from '@molstar/core/util/sniffable';
import { ChunkedBigString, type StringLike } from '@molstar/core/util/string-like';
import { guessCifVariant } from '../provider.js';

/** The former text sniff of a `.cif`: one search per needle */
function reference(text: string) {
  const has = (needle: string) => text.includes(needle);
  if (text.startsWith('data_SERVER\n#\n_density_server_result')) return 'dscif';
  if (text.startsWith('data_SERVER\n#\ndata_SEGMENTATION_DATA')) return 'segcif';
  if (has('_volume_data_3d_info.')) {
    if (has('_volume_data_3d.')) return 'dscif';
    if (has('_segmentation_data_3d.')) return 'segcif';
  }
  if (has('_refln.pdbx_FWT') || has('_refln.pdbx_DELFWT')) return 'sfcif';
  if (has('atom_site_fract_x') || has('atom_site.fract_x')) return 'coreCif';
  return -1;
}

const Fragments = [
  'data_x\n',
  '_atom_site.id 1\n',
  '_atom_site.Cartn_x 1.0\n',
  '_atom_site_fract_x 0.5\n',
  '_atom_site.fract_x 0.5\n',
  'atom_site',
  'atom_site.fract_',
  '_refln.pdbx_',
  '_refln.pdbx_FWT 1\n',
  '_refln.pdbx_DELFWT 1\n',
  '_refln.pdbx_DEL',
  '_refln.pdbx_FW',
  '_refln.pdbx_refln.pdbx_FWT',
  '_volume_data_3d_info.name a\n',
  '_volume_data_3d.values 1\n',
  '_segmentation_data_3d.v 1\n',
  'ATOM 1 C CA ALA A 1 \n',
  'é😀 ',
];

/** Deterministic LCG in [0, 1) */
function rng(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

describe('guessCifVariant', () => {
  const info = getFileNameInfo('x.cif');

  it('tells the same variant as one search per needle, for bytes, responses, strings and big strings (perf-audit-3 01/P13)', () => {
    const next = rng(7);
    const enc = new TextEncoder();
    const texts = ['', 'data_SERVER\n#\n_density_server_result.x 1\n', 'data_SERVER\n#\ndata_SEGMENTATION_DATA\n'];
    for (let t = 0; t < 400; t++) {
      let text = '';
      const n = 1 + Math.floor(next() * 6);
      for (let i = 0; i < n; i++) text += Fragments[Math.floor(next() * Fragments.length)];
      texts.push(text);
    }
    const seen = new Set<unknown>();
    for (const text of texts) {
      const expected = reference(text);
      seen.add(expected);
      const bytes = enc.encode(text);
      const forms: (Uint8Array | StringLike | Sniffable)[] = [
        bytes,
        Sniffable.of(enc.encode(text)),
        text,
        // chunks of 64 characters: needles cross chunk ends
        ChunkedBigString.fromString(text, 6),
      ];
      for (const data of forms) expect(guessCifVariant(info, data), text).toBe(expected);
    }
    // every variant is told by some text
    expect([...seen].sort()).toEqual([-1, 'coreCif', 'dscif', 'segcif', 'sfcif'].sort());
  });

  it('reads a plain mmCIF three times, not five', () => {
    const bytes = new TextEncoder().encode(`data_x\n${'_atom_site.Cartn_x 1.0\nATOM 1 C CA ALA A 1 \n'.repeat(50)}`);
    const indexOf = vi.spyOn(Uint8Array.prototype, 'indexOf');
    try {
      expect(guessCifVariant(info, bytes)).toBe(-1);
      // a search that finds nothing ends with a native scan to the end
      expect(indexOf.mock.results.filter((r) => r.value === -1).length).toBe(3);
    } finally {
      indexOf.mockRestore();
    }
  });
});
