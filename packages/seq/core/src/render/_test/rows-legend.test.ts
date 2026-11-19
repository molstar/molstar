/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { viewLegends } from '../index.js';
import { settle } from './helpers.js';

// The rows layer's letters legend lists the letters of its alphabet (an RNA's U, never its T, a DNA's T), the letters
// drawn when they are cells, the alphabet's own when they are image tiles.

function spec(alphabet: 'rna' | 'dna' | 'protein', letters: string[]) {
  const b = createSeqSpec();
  const aln = b.data.alignment('a', {
    alphabet,
    rows: { label: letters.map((_, i) => `r${i}`), letters },
  });
  b.views.view('v', (v) => {
    const s = v.section('s', { of: aln });
    s.track('rows').layer('rows', { type: 'rows', data: aln });
  });
  return b.build();
}

async function labels(alphabet: 'rna' | 'dna' | 'protein', letters: string[], width = 640) {
  const st = await settle(spec(alphabet, letters));
  const f = await st.frame({ size: { width, height: 200 }, dpr: 1, settled: true });
  const legend = viewLegends(f.layout).legends.find((l) => l.layers.some((k) => k.endsWith('/rows')));
  return legend && legend.entries.kind === 'swatches' ? legend.entries.items.map((i) => i.label).sort() : [];
}

describe('the rows layer’s legend', () => {
  it('an RNA alignment: U, not T', async () => {
    expect(await labels('rna', ['ACGUACGU', 'ACGUACGU'])).toEqual(['A', 'C', 'G', 'U']);
  });

  it('a DNA alignment: T, not U', async () => {
    expect(await labels('dna', ['ACGTACGT', 'ACGTACGT'])).toEqual(['A', 'C', 'G', 'T']);
  });

  it('zoomed out to image tiles an RNA alignment still says U (the alphabet’s own)', async () => {
    const row = 'ACGU'.repeat(400);
    expect(
      await labels(
        'rna',
        Array.from({ length: 40 }, () => row),
        300,
      ),
    ).toEqual(['A', 'C', 'G', 'U']);
  });
});
