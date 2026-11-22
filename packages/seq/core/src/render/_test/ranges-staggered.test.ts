/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A source's track while its tiles arrive one by one, each after a delay of its own, on a fake clock (_spec/README.md
// "Render", Track heights): its height changes once per settled read, the lanes of the whole region from the first
// (a pan to the other tiles does not grow it), and a band still loading draws no label outside its box (EC-169: a
// pending tile is never drawn as empty, and its placeholder lane is the band's).

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { delayedRuntime } from '@molstar/seqx-core/data';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { GffRuntime } from '../../runtimes/gff.js';
import type { RangesLay } from '../index.js';
import { settle } from './helpers.js';
import { layerLay, outside } from './small-helpers.js';

const SIZE = { width: 1000, height: 300 } as const;
const Tiles = 6;

beforeEach(() => {
  vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout', 'Date', 'performance'] });
});
afterEach(() => {
  vi.useRealTimers();
});

/** Six tiles of 1,024 positions, tile `t` holding `t + 1` overlapping transcripts with names longer than their marks. */
function spec(line: string) {
  const lines = ['##gff-version 3', '##sequence-region c 1 8000'];
  for (let t = 0; t < Tiles; t++)
    for (let k = 0; k <= t; k++)
      lines.push(
        [
          'c',
          '.',
          'mRNA',
          `${1024 * t + 100 + k}`,
          `${1024 * t + 160 + k}`,
          '.',
          '+',
          '.',
          `ID=x${t}-${k};Name=TP53-${220 + t * 10 + k}`,
        ].join('\t'),
      );
  const b = createSeqSpec({ title: 'staggered' });
  const file = b.data.source('gff', { format: 'gff3', value: lines.join('\n') as never });
  const chr = b.data.sequence('c', { alphabet: 'dna', source: file, name: 'c' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    v.section('c', { of: chr })
      .track('tx', {})
      .layer('tx', { type: 'ranges', data: tx, ...(line === 'gene' ? { line } : {}) } as never);
  });
  return b.build();
}

/** The runtime answering each tile after a delay of its own: later tiles later. */
const Staggered = delayedRuntime(GffRuntime, {
  delayMs: (r) => 100 + 150 * Math.floor(r.start / 1024),
});

async function load(line: string, window: readonly [number, number]) {
  const s = await settle(spec(line), [Staggered]);
  const heights: number[] = [];
  const unclipped: number[] = [];
  let pendingFrames = 0;
  const frame = (viewport: { start: number; end: number }) => {
    const f = s.renderer.frame({ size: SIZE, settled: s.store.settled.value, viewport } as never);
    s.store.want(f.wants);
    const h = f.layout.sections[0].tracks[0].box.height;
    if (heights[heights.length - 1] !== h) heights.push(h);
    const lay = layerLay<RangesLay>(f, 0);
    if (lay.pending.length) {
      pendingFrames++;
      unclipped.push(outside(f.lists, 'labels', lay.box.y, lay.box.y + lay.box.height).bad);
    }
    return f;
  };
  const run = async (viewport: { start: number; end: number }, ms: number) => {
    for (let t = 0; t < ms; t += 25) {
      await vi.advanceTimersByTimeAsync(25);
      frame(viewport);
    }
  };
  frame({ start: window[0], end: window[1] });
  await run({ start: window[0], end: window[1] }, 3000);
  return { s, heights, unclipped, pendingFrames, run };
}

describe('a source’s track while its tiles arrive after delays of their own', () => {
  for (const line of ['blocks', 'gene']) {
    it(`${line}: a view of one tile first, a pan to the others: the height changes once, to the whole region’s lanes`, async () => {
      const { heights, run } = await load(line, [0, 1024]);
      // the placeholder lane, then the six lanes the last tile packs (22 and six lanes' 92 px): once
      expect(heights).toEqual([22, 92]);
      for (const [start, end] of [
        [1024, 2048],
        [3072, 4096],
        [5120, 6144],
        [0, 8000],
      ]) {
        await run({ start, end }, 1500);
        expect(heights).toEqual([22, 92]);
      }
    });

    it(`${line}: EC-169, a pending track draws its labels inside its band, never below it`, async () => {
      // the whole region in view: tiles arrive one after another under a one-lane placeholder band
      const { unclipped, pendingFrames } = await load(line, [0, 8000]);
      expect(pendingFrames).toBeGreaterThan(5);
      expect(unclipped.every((n) => n === 0)).toBe(true);
    });
  }
});
