/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A band's label passes are clipped to its box (layouts.md, Loading: the labels stay in the band): a scrolled band's
// labels, leaders and backdrops cut at its top and bottom (TP53-233 above, TP53-221 below), a pending track's
// placeholder lane's labels not below it (TP53-227).

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import type { RangesLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay, outside } from './small-helpers.js';

const SIZE = { width: 1000, height: 200 } as const;

function spec(line: string) {
  const lines = ['##gff-version 3', '##sequence-region c 1 2000'];
  for (let t = 0; t < 20; t++)
    lines.push(
      ['c', '.', 'mRNA', `${101 + t}`, `${900 + t}`, '.', '+', '.', `ID=T${t};Name=TP53-${220 + t}`].join('\t'),
      ['c', '.', 'exon', `${101 + t}`, `${300 + t}`, '.', '+', '.', `Parent=T${t}`].join('\t'),
      ['c', '.', 'exon', `${700 + t}`, `${900 + t}`, '.', '+', '.', `Parent=T${t}`].join('\t'),
    );
  const b = createSeqSpec({ title: 'many' });
  const file = b.data.source('gff', { format: 'gff3', value: lines.join('\n') as never });
  const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    v.section('c', { of: chr })
      .track('tx', { height: 60 } as never)
      .layer('tx', { type: 'ranges', data: tx, ...(line === 'gene' ? { line } : {}) } as never);
  });
  return b.build();
}

describe('a band’s labels stay in its box', () => {
  for (const line of ['blocks', 'gene'])
    for (const top of [0, 5, 100]) {
      it(`${line}: a scrolled band (top ${top}) draws no label outside it`, async () => {
        const s = await settle(spec(line), [GffRuntime], { track: { maxPx: 240 } });
        const f = await s.frame({
          size: SIZE,
          settled: true,
          viewports: { '/views/main/c': { start: 0, end: 1000, bands: { '/views/main/c/tx': { bandTop: top } } } },
        } as never);
        const lay = layerLay<RangesLay>(f, 0);
        expect(lay.scrolls).toBe(true);
        const r = outside(f.lists, 'labels', lay.box.y, lay.box.y + lay.box.height);
        expect(r.drawn).toBeGreaterThan(0);
        expect(r.bad).toBe(0);
      });
    }

  it('a pending track draws no label below its band', async () => {
    const s = await settle(spec('gene'), [GffRuntime], { track: { maxPx: 240 } });
    const f = s.renderer.frame({ size: SIZE });
    const lay = layerLay<RangesLay>(f, 0);
    expect(outside(f.lists, 'labels', lay.box.y, lay.box.y + lay.box.height).bad).toBe(0);
  });
});
