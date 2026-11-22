/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A source's track is as tall as its reachable region needs from the first settled frame (_spec/README.md "Render",
// Track heights): viewing a stretch of few lanes first, then panning to one of many, never grows or shrinks the track.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { DefaultOpenLimits } from '@molstar/seqx-core/data';
import { GffRuntime } from '../../runtimes/gff.js';
import { memoryAssets } from '../../runtimes/_test/bed-fixtures.js';
import { settle } from './helpers.js';

const SIZE = { width: 1000, height: 300 } as const;

function gffText(clusters: readonly (readonly [number, number])[]) {
  const lines = ['##gff-version 3', '##sequence-region c 1 20000'];
  // two overlapping transcripts at the start, six at the end
  for (const [at, count] of clusters)
    for (let t = 0; t < count; t++)
      lines.push(['c', '.', 'mRNA', `${at + t}`, `${at + 500 + t}`, '.', '+', '.', `ID=t${at}-${t}`].join('\t'));
  return lines.join('\n');
}

function spec(
  text = gffText([
    [100, 2],
    [18_000, 6],
  ]),
  sections = 1,
) {
  const b = createSeqSpec({ title: 'two clusters' });
  const file = b.data.source('gff', { format: 'gff3', value: text as never });
  const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    for (let k = 0; k < sections; k++)
      v.section(k ? `c${k}` : 'c', { of: chr })
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

  it('a data swap on one source starts again: the track follows the new data’s lanes, grown and shrunk', async () => {
    const s = await settle(spec(), [GffRuntime]);
    const heightIn = async (end = 20_000) => {
      const f = await s.frame({
        size: SIZE,
        settled: true,
        viewports: { '/views/main/c': { start: 0, end } },
      } as never);
      return f.layout.sections[0].tracks[0].box.height;
    };
    const six = await heightIn();
    expect(six).toBeGreaterThanOrEqual(6 * 14 + 8);
    const swap = async (text: string) => {
      const address = [...s.store.named()][0].address;
      const opened = await GffRuntime.open(
        { address, format: 'gff3', value: text as never },
        { assets: memoryAssets({}), signal: new AbortController().signal, limits: DefaultOpenLimits },
      );
      expect(s.store.reattach(address, opened)).toBe(true);
    };
    // two lanes only: a pan never shrinks the track, a new reading of the source does
    await swap(gffText([[100, 2]]));
    expect(await heightIn()).toBeLessThan(six);
    await swap(gffText([[18_000, 8]]));
    expect(await heightIn()).toBeGreaterThanOrEqual(8 * 14 + 8);
  });

  it('a lane the reach’s packing does not model (a label beside a short mark) still grows the track', async () => {
    // overlapping marks of 1 px (a name beside each, in its own room): the view packs more lanes than first-fit by start
    const lines = ['##gff-version 3', '##sequence-region c 1 20000'];
    for (let t = 0; t < 8; t++)
      lines.push(
        ['c', '.', 'mRNA', `${100 + 2 * t}`, `${102 + 2 * t}`, '.', '+', '.', `ID=x${t};Name=a-long-name-${t}`].join(
          '\t',
        ),
      );
    const s = await settle(spec(lines.join('\n')), [GffRuntime]);
    const f = await s.frame({
      size: SIZE,
      settled: true,
      viewports: { '/views/main/c': { start: 0, end: 1000 } },
    } as never);
    const lay = f.layout.sections[0].tracks[0].layers[0].lay as { laneCount: number };
    const h = f.layout.sections[0].tracks[0].box.height;
    // whatever the view packs, the track is at least as tall as its lanes
    expect(h).toBeGreaterThanOrEqual(lay.laneCount * 14 + 8);
    const fit = await s.frame({ size: SIZE, settled: true } as never);
    expect(fit.layout.sections[0].tracks[0].box.height).toBe(h);
  });

  it('two sections over one source size alike, each by its own layer (a prep is a layer’s)', async () => {
    const s = await settle(spec(undefined, 2), [GffRuntime]);
    const f = await s.frame({
      size: SIZE,
      settled: true,
      viewports: { '/views/main/c': { start: 0, end: 1500 }, '/views/main/c1': { start: 17_500, end: 19_000 } },
    } as never);
    const a = f.layout.sections[0].tracks[0].box.height;
    expect(a).toBeGreaterThanOrEqual(6 * 14 + 8);
    expect(f.layout.sections[1].tracks[0].box.height).toBe(a);
  });
});
