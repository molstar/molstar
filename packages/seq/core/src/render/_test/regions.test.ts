/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type SeqData, SeqRuntime, seqToValue } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { P53 } from './fixtures.js';
import { settle } from './helpers.js';

// Planning by region: a values layer zoomed out reads a level of the store's value pyramid (not level 1, which would
// cover more tiles than one want may), and a protein's view wants the gene's region its residues' codons cover,
// carried back through the static codon map.

const coverageOnly = (name: string, length: number, shapes: ('letters' | 'values')[]): SeqData => ({
  coverage: { systems: [{ name, length }], shapes: shapes.map((shape) => ({ system: name, shape })), levels: [1] },
  pairs: [],
});

describe('planning by region', () => {
  it('a values layer zoomed out reads a coarse level of the pyramid, within the tile bound', async () => {
    const length = 300_000;
    const b = createSeqSpec({ title: 'signal' });
    const src = b.data.source('chrom', { format: 'seq', value: seqToValue(coverageOnly('s', length, ['letters'])) });
    const s = b.data.sequence('s', { alphabet: 'dna', source: src });
    const v = b.data.values('signal', { target: s, values: Float32Array.from({ length }, (_, i) => (i % 100) / 100) });
    b.views.view('main', (view) => view.section('all', { of: s }).track('signal').layer({ type: 'values', data: v }));
    const { store, frame } = await settle(b.build(), [SeqRuntime]);
    const f = await frame({ size: { width: 640, height: 60 } });
    const levels = new Set(f.wants.map((w) => w.level));
    expect([...levels]).toEqual([256]);
    expect(f.layout.issues).toEqual([]);
    expect(store.issues.value).toEqual([]);
    const keys = store.keysOf(f.wants[0]);
    expect(keys.length).toBeGreaterThan(0);
    expect(keys.every((k) => store.get(k)?.status === 'loaded')).toBe(true);
    // zoomed in to a few hundred positions, the layer reads level 1 again
    const near = await frame({ size: { width: 640, height: 60 }, viewport: { start: 1000, end: 1300 } });
    expect(new Set(near.wants.map((w) => w.level))).toEqual(new Set([1]));
  });

  it('a protein’s view wants the gene’s region its codons cover, through the static codon map', async () => {
    const b = createSeqSpec({ title: 'p53 on its gene' });
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    const file = b.data.source('gene', { format: 'seq', value: seqToValue(coverageOnly('tp53', 1500, ['values'])) });
    const gene = b.data.sequence('tp53', { alphabet: 'dna', source: file });
    const exons = { start: [100, 480, 960], end: [400, 900, 1422] };
    b.data.map('cds', { from: p53, to: gene, method: 'codon', exons, strand: '-' });
    const signal = b.data.values('gc', { target: gene, source: file });
    b.views.view('main', (v) =>
      v
        .section('protein', { of: p53, range: [94, 292] })
        .track('gc')
        .layer({ type: 'values', data: signal }),
    );
    const { frame } = await settle(b.build(), [SeqRuntime]);
    const f = await frame({ size: { width: 640, height: 60 } });
    const view = f.wants.filter((w) => w.priority === 'view').map((w) => [w.source, w.system, w.start, w.end]);
    // residues 94 to 291 on the minus strand: the end of the second exon and the start of the third
    expect(view).toEqual([
      ['/data/gene', 'tp53', 486, 900],
      ['/data/gene', 'tp53', 960, 1140],
    ]);
  });
});
