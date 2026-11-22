/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A source's track is as tall as its reachable region needs from the first settled frame (_spec/README.md "Render",
// Track heights): viewing a stretch of few lanes first, then panning to one of many, never grows or shrinks the track.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import { settle } from './helpers.js';

const SIZE = { width: 1000, height: 300 } as const;

function spec() {
  const lines = ['##gff-version 3', '##sequence-region c 1 20000'];
  // two overlapping transcripts at the start, six at the end
  for (const [at, count] of [
    [100, 2],
    [18_000, 6],
  ])
    for (let t = 0; t < count; t++)
      lines.push(['c', '.', 'mRNA', `${at + t}`, `${at + 500 + t}`, '.', '+', '.', `ID=t${at}-${t}`].join('\t'));
  const b = createSeqSpec({ title: 'two clusters' });
  const file = b.data.source('gff', { format: 'gff3', value: lines.join('\n') as never });
  const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    v.section('c', { of: chr })
      .track('tx', {})
      .layer('tx', { type: 'ranges', data: tx } as never);
  });
  return b.build();
}

describe('a source’s track over its reachable region', () => {
  it('is as tall as the whole region’s lanes at the first settled frame, whichever stretch is in view', async () => {
    const s = await settle(spec(), [GffRuntime]);
    const heightIn = async (start: number, end: number) => {
      const f = await s.frame({ size: SIZE, settled: true, viewports: { '/views/main/c': { start, end } } } as never);
      return f.layout.sections[0].tracks[0].box.height;
    };
    const first = await heightIn(0, 1500);
    expect(await heightIn(17_000, 19_000)).toBe(first);
    expect(await heightIn(0, 20_000)).toBe(first);
    expect(await heightIn(0, 1500)).toBe(first);
    // six lanes (a lane and its padding, 14 and 4+4 px) rather than the two in view
    expect(first).toBeGreaterThanOrEqual(6 * 14 + 8);
  });
});
