/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An item that starts in an earlier tile is drawn when the view lies wholly past that tile's edge: a region's answer
// holds the items overlapping it, so every item in view has its marks, once.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import type { RangesLay } from '../index.js';
import { BigRuntime } from './fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 720, height: 400 };

function spec(layout?: 'text') {
  const b = createSeqSpec({ title: 'tile edge' });
  const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
  const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
  const items = b.data.ranges('items', { target: chr, source: src });
  // biome-ignore lint/suspicious/noExplicitAny: the builder's view type is not exported
  const fill = (v: any) => {
    const sec = v.section('chr', { of: chr });
    sec.track('letters').layer('letters', { type: 'letters' });
    sec.track('items').layer('items', { type: 'ranges', data: items });
  };
  if (layout) b.views.view('main', { layout }, fill);
  else b.views.view('main', fill);
  return b.build();
}

/** The synthetic source's items: every 50 positions, 70 long, `r<start>`; the ones overlapping `[a, b)`. */
function overlapping(a: number, b: number): string[] {
  const out: string[] = [];
  for (let p = Math.max(0, Math.ceil((a - 69) / 50) * 50); p < b; p += 50) if (p + 70 > a) out.push(`r${p}`);
  return out;
}

const keysOf = (l: RangesLay) => Array.from({ length: l.items }, (_, j) => String(l.keys[j]));

describe('items from an earlier tile', () => {
  it('a tracks view wholly past a tile edge draws every item overlapping it, once', async () => {
    const s = await settle(spec(), [BigRuntime]);
    for (const [a, b] of [
      [1040, 1080],
      [2100, 2110],
      [1030, 1040],
    ]) {
      const f = await s.frame({ size: Size, viewport: { start: a, end: b } });
      const keys = keysOf(f.layout.sections[0].tracks[1].layers[0].lay as RangesLay);
      expect(new Set(keys).size).toBe(keys.length);
      // the layout may hold more (whole tiles, clipped at draw)
      for (const k of overlapping(a, b)) expect(keys).toContain(k);
    }
  });

  it('a text view of lines past a tile edge draws every item overlapping each, once', async () => {
    const s = await settle(spec('text'), [BigRuntime], { text: { collapseEmpty: false } });
    const f = await s.frame({ size: Size, viewport: { start: 1040, end: 1040 } });
    const sec = f.layout.sections[0];
    const t = sec.placement as unknown as { strips: { n: number; bands: number; axes: { u0: number; u1: number }[] } };
    let checked = 0;
    for (let i = 0; i < t.strips.n; i++) {
      if (!sec.lays[i * t.strips.bands].length) continue;
      if (!sec.lays[i * t.strips.bands + 1].length) continue; // the items' band below the box
      const a = t.strips.axes[i];
      const keys = keysOf(sec.lays[i * t.strips.bands + 1][0].lay as RangesLay);
      expect(new Set(keys).size).toBe(keys.length);
      for (const k of overlapping(a.u0, a.u1)) expect(keys).toContain(k);
      checked++;
    }
    expect(checked).toBeGreaterThan(3);
  });
});
