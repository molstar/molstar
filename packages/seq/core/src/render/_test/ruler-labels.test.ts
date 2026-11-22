/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ruler's labels in bases are kept between layouts: a pan's frame over the same ticks formats none again
// (`formatGenomic` parses a format and builds its string: 16 KB a frame of p53's chromosome pan, the storybook's
// `chromosome.quick.ts`). Mocks `@molstar/core/util/number` to count its calls (the file runs isolated).

import { SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it, vi } from 'vitest';
import { type RulerLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

const calls = vi.hoisted(() => ({ n: 0 }));

vi.mock('@molstar/core/util/number', async (importOriginal) => {
  const m = await importOriginal<typeof import('@molstar/core/util/number')>();
  return {
    ...m,
    formatGenomic: (n: number, step: number) => {
      calls.n++;
      return m.formatGenomic(n, step);
    },
  };
});

function chromosomeRuler() {
  const b = createSeqSpec({ title: 'ruler' });
  const value = seqToValue({
    coverage: { systems: [{ name: 'chr', length: 83_257_441 }], shapes: [], levels: [1] },
    pairs: [],
  });
  const source = b.data.source('chr-source', { format: 'seq', value: value as never });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source, name: 'chr' });
  b.views.view('main', (v) => {
    v.section('chr', { of: chr }).track('r').layer('r', { type: 'ruler', units: 'bases' });
  });
  return b.build();
}

describe('the ruler’s labels in bases', () => {
  it('a pan over the same ticks formats none of them again', async () => {
    const s = await settle(chromosomeRuler(), [SeqRuntime]);
    const size = { width: 700, height: 60 };
    const at = (d: number) => ({ start: 7_600_000 + d, end: 7_700_000 + d });
    const first = layerLay<RulerLay>(await s.frame({ size, viewport: at(0) }), 0);
    expect(first.labels[0]).toMatch(/ Mb$/);
    expect(calls.n).toBeGreaterThan(0);
    calls.n = 0;
    for (let d = 1; d <= 8; d++) await s.frame({ size, viewport: at(d * 100) });
    expect(calls.n).toBe(0);
  });
});
