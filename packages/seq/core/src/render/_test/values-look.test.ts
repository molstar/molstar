/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The values layer's look (_spec/README.md "Render", Layer types, `values`): bars with a gap and rounded free ends from
// 4 px a value, a profile below, 10% headroom, a hairline baseline, a missing value a gap; at the edges of the window, on
// both sides of zero, with a per-value color, and in the line encodings.

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { LightTokens, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { readPaths, readRects } from './small-helpers.js';

const Size = { width: 400, height: 56 };

/** A values layer over `values` (a section `s` of as many positions) in one track of 32 px, rendered. */
async function render(
  values: readonly number[],
  opts: {
    viewport?: { start: number; end: number };
    dpr?: number;
    size?: { width: number; height: number };
    layer?: Record<string, unknown>;
    config?: Parameters<typeof renderSvg>[1] extends infer O ? (O extends { config?: infer C } ? C : never) : never;
  } = {},
) {
  const b = createSeqSpec({ title: 'look' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'A'.repeat(values.length) });
  const v = b.data.values('v', { target: s, values: Float32Array.from(values) });
  b.views.view('main', (view) =>
    view
      .section('s', { of: s })
      .track('v')
      .layer('values', { type: 'values', data: v, ...opts.layer }),
  );
  return renderSvg(b.build(), {
    size: opts.size ?? Size,
    ...(opts.viewport ? { viewport: opts.viewport } : {}),
    ...(opts.dpr ? { dpr: opts.dpr } : {}),
    ...(opts.config ? { config: opts.config } : {}),
  });
}

/** The data pass's bars: the rounded rects of the layer (the series color at 60%). */
const barsOf = (r: Awaited<ReturnType<typeof render>>) => readRects(r.lists).filter((x) => (x.fill & 255) === 0x99);

describe('the values look', () => {
  it('a negative bar is round at its free end, below the baseline (the top corners square), a positive at its top', async () => {
    // 6 positions in 400 px: bars. The domain [-3, 3]: the baseline in the middle of the headroom's share
    const r = await render([3, -3, 3, -3, 3, -3]);
    const lines = recordDrawLists(r.lists, ['data']).split('\n');
    const bars = lines.filter((l) => l.startsWith('rect ') && / #[0-9a-f]{6}99 /.test(l));
    expect(bars.length).toBe(6);
    // a positive's bottom corners are square (cut-bottom), a negative's top corners are (cut-top)
    expect(bars.filter((l) => l.endsWith('cut-bottom')).length).toBe(3);
    expect(bars.filter((l) => l.endsWith('cut-top')).length).toBe(3);
  });

  it('a column wholly below zero draws its area once: the band from the baseline to its maximum is not doubled', async () => {
    // 2,000 positions in 400 px: columns. Every column but the first spans [-3, -1]; the domain [-3, 3] has its zero
    const r = await render(
      Array.from({ length: 2000 }, (_, i) => (i === 0 ? 3 : i % 2 ? -1 : -3)),
      { size: { width: 400, height: 56 } },
    );
    expect(readPaths(r.lists).filter((p) => p.closed).length).toBeGreaterThan(0);
    const baseline = readRects(r.lists).find((x) => x.fill === LightTokens.colors['seq-sequence'] && x.w > 100)!;
    const blocks = readRects(r.lists).filter((x) => (x.fill & 255) === 0x59 && x.x > baseline.x + 3);
    expect(blocks.length).toBeGreaterThan(0);
    // each column's maximum (-1) is about 4.8 px under the baseline: its block begins there, not at the baseline the
    // area (baseline to maximum) already covers
    for (const x of blocks) expect(x.y).toBeGreaterThan(baseline.y + 3);
  });

  it('a bar cut by the window’s edge has a positive width and radius, canvas and export, at dpr 1 and 2', async () => {
    for (const dpr of [1, 2]) {
      const r = await render(
        Array.from({ length: 200 }, (_, i) => 1 + (i % 5)),
        { viewport: { start: 100.9, end: 170.9 }, dpr, size: { width: 400, height: 56 } },
      );
      expect(r.svg).not.toMatch(/ (width|height)="-/);
      // the sliver of the position the edge cuts is drawn too (a gap wider than it dropped it): 100 to 170 are 71 bars
      expect(barsOf(r).length, `dpr ${dpr}`).toBe(71);
      for (const x of readRects(r.lists)) {
        expect(x.w, `dpr ${dpr}`).toBeGreaterThanOrEqual(0);
        expect(x.h, `dpr ${dpr}`).toBeGreaterThanOrEqual(0);
      }
      const radii = recordDrawLists(r.lists, ['data'])
        .split('\n')
        .map((l) => /squircle (-?[\d.]+)/.exec(l)?.[1])
        .filter((q): q is string => q !== undefined);
      for (const q of radii) expect(Number(q), `dpr ${dpr}`).toBeGreaterThanOrEqual(0);
    }
  });

  it('a per-value color keeps the profile a few paths, not two a column, and its top line continuous', async () => {
    // 1,200 positions in 400 px: a profile; the value picks a ramp: every column its own color
    const values = Array.from({ length: 1200 }, (_, i) => 5 + 4 * Math.sin(i / 90));
    const b = createSeqSpec({ title: 'colored' });
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'A'.repeat(values.length) });
    const v = b.data.values('v', { target: s, values: Float32Array.from(values) });
    b.views.view('main', (view) =>
      view
        .section('s', { of: s })
        .track('v')
        .layer('values', { type: 'values', data: v, color: theme.field('value', { domain: [0, 10] }) }),
    );
    const r = await renderSvg(b.build(), { size: Size });
    const paths = readPaths(r.lists);
    // far fewer than a column (400 of them, two paths each)
    expect(paths.length).toBeLessThan(80);
    // the lines (open paths) join: each run's line starts where the one before it ended
    const lines = paths.filter((p) => !p.closed);
    for (let i = 1; i < lines.length; i++) {
      const a = lines[i - 1].xy;
      expect(lines[i].xy[0]).toBeCloseTo(a[a.length - 2], 3);
    }
  });

  it('a color’s own alpha is kept: the bars multiply it (a 50% series at 60% is 30%)', async () => {
    const r = await render([1, 2, 3, 2, 1, 2], { layer: { color: '#5b7c9980' } });
    const bars = recordDrawLists(r.lists, ['data'])
      .split('\n')
      .filter((l) => l.startsWith('rect ') && l.includes('#5b7c99'));
    expect(bars.length).toBeGreaterThan(0);
    // 0x80 × 153 / 255 = 77 = 0x4d
    for (const l of bars) if (l.includes('squircle')) expect(l).toContain('#5b7c994d');
  });

  it('the line and area share the headroom, the baseline and the gap of the bars: aligned on one scale', async () => {
    const values = [1, 2, 3, Number.NaN, 2, 3, 1];
    const bars = await render(values);
    const line = await render(values, { layer: { encoding: 'line' } });
    const top = (r: Awaited<ReturnType<typeof render>>) => Math.min(...readRects(r.lists).map((x) => x.y));
    const baseline = (r: Awaited<ReturnType<typeof render>>) =>
      readRects(r.lists).filter((x) => x.fill === LightTokens.colors['seq-sequence'] && x.w > 100);
    expect(baseline(line).length).toBe(1);
    expect(baseline(line)[0].y).toBeCloseTo(baseline(bars)[0].y, 4);
    // the highest value (3) sits at the same y on the line's path as at the bar's top
    const lineTop = Math.min(...readPaths(line.lists).flatMap((p) => p.xy.filter((_, i) => i % 2 === 1)));
    expect(lineTop).toBeCloseTo(Math.min(...barsOf(bars).map((x) => x.y)), 3);
    // a missing value: a gap, no mark in the missing color
    expect(top(line)).toBeGreaterThan(0);
    expect(line.svg).not.toMatch(/#aeb4bb|seq-missing/);
    expect(readPaths(line.lists).filter((p) => !p.closed).length).toBe(2);
  });
});
