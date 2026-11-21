/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A ranges layer over a source (its whole data's lanes not known) with more items than its band's lanes
// (_spec/README.md "Render", Track heights): the lanes are centred by the most the layer has shown, the last lane the
// density row, every lane inside the band (the lanes were centred by the one lane of the first frame, before any tile
// arrived, and ran past the band's bottom).

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import type { RangesLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay } from './small-helpers.js';

const SIZE = { width: 1000, height: 120 } as const;

describe('a source’s items past the lanes', () => {
  for (const line of ['blocks'])
    it(`${line}: a source’s items past its lanes: every lane inside the band, the last the density row`, async () => {
      // 12 transcripts over one another in a band of a few lanes: the rest go to the density row
      const lines = ['##gff-version 3', '##sequence-region c 1 2000'];
      for (let t = 0; t < 12; t++)
        lines.push(
          ['c', '.', 'mRNA', `${101 + t}`, `${900 + t}`, '.', '+', '.', `ID=x${t}`].join('\t'),
          ['c', '.', 'exon', `${101 + t}`, `${300 + t}`, '.', '+', '.', `Parent=x${t}`].join('\t'),
          ['c', '.', 'exon', `${700 + t}`, `${900 + t}`, '.', '+', '.', `Parent=x${t}`].join('\t'),
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
      const s = await settle(b.build(), [GffRuntime]);
      const f = await s.frame({ size: SIZE });
      const lay = layerLay<RangesLay>(f, 0);
      expect(lay.overflow).toBeGreaterThan(0);
      expect(lay.m).toBeGreaterThan(0);
      for (let j = 0; j < lay.m; j++) {
        const y = lay.laneY + lay.mlane[j] * lay.laneHeight;
        expect(y).toBeGreaterThanOrEqual(lay.box.y - 0.01);
        expect(y + lay.laneHeight).toBeLessThanOrEqual(lay.box.y + lay.box.height + 0.01);
      }
      expect(lay.densityY + lay.densityH).toBeLessThanOrEqual(lay.box.y + lay.box.height + 0.01);
    });
});
