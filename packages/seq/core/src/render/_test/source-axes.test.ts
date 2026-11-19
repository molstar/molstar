/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Axes over sequences whose letters come from a source (_spec/README.md "Resolve", Systems): a grid's values are read
// once its axes' lengths are the source's, not dropped at resolve; a section of blocks over such sequences has its
// axis from their coverages, never pending once the sources are open.

import { inlineSeq, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { describe, expect, it } from 'vitest';
import type { GridLay } from '../index.js';
import { settle } from './helpers.js';

type Builder = ReturnType<typeof createSeqSpec>;

/** A `seq` source of one protein record. */
function file(b: Builder, id: string, letters: string) {
  return b.data.source(`${id}-file`, {
    format: 'seq',
    value: seqToValue(
      inlineSeq(
        { name: id, length: letters.length, alphabet: 'protein' },
        { letters: { form: 'letters', start: 0, letters } },
      ),
    ),
  });
}

describe('axes over a source’s letters', () => {
  it('a grid over them is read once the source is open: its values held, its cells laid out', async () => {
    const b = createSeqSpec();
    const a = b.data.sequence('a', { alphabet: 'protein', source: file(b, 'a', 'ACDEFGHIKLMN') });
    const values = new Float32Array(144).map((_, k) => (k % 13 === 0 ? 1 : 0));
    const grid = b.data.grid('g', { x: a, values });
    b.views.view('v', (v) =>
      v.section('plot', { of: a, layout: 'grid' }).track('t').layer('g', { type: 'grid', data: grid }),
    );
    const spec = b.build();
    // at resolve, before the source: no length, so no values yet
    expect(resolveSpec(spec, { runtimes: [SeqRuntime] }).annotations.get('/data/g')?.grid).toBeUndefined();
    const s = await settle(spec, [SeqRuntime]);
    expect(s.model.annotations.get('/data/g')?.grid).toMatchObject({ nx: 12, ny: 12 });
    expect(s.model.issues).toEqual([]);
    const f = await s.frame({ size: { width: 400, height: 400 }, settled: true });
    expect(f.layout.issues).toEqual([]);
    expect((f.layout.sections[0].tracks[0].layers[0].lay as GridLay).grid).toMatchObject({ nx: 12, ny: 12 });
  });

  it('a section of blocks over them has its axis from their coverages', async () => {
    const b = createSeqSpec();
    const a = b.data.sequence('a', { alphabet: 'protein', source: file(b, 'a', 'ACDEFGHIKLMN') });
    const c = b.data.sequence('c', { alphabet: 'protein', source: file(b, 'c', 'PQRSTVWY') });
    b.views.view('v', (v) => {
      const s = v.section('complex', {});
      s.block('a', { of: a });
      s.block('c', { of: c });
      s.track('l').layer('l', { type: 'letters' });
    });
    const s = await settle(b.build(), [SeqRuntime]);
    expect(s.view.sections[0].axis.ranges).toEqual([0, 20]);
    const f = await s.frame({ size: { width: 600, height: 200 }, settled: true });
    expect(f.layout.issues).toEqual([]);
    expect(f.layout.sections[0].axis?.total).toBe(20);
  });
});
