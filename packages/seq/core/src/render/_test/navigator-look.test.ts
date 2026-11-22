/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigator's look (navigators.md, "The window's look: one control row, one stroke, one readout place"): over a
// genome two rows always, whatever the window's width (a chromosome whole or a few Mb of it: chr6 and chr13 alike), the
// grips on the bottom row only; a window at the fit one quiet frame with no grips, its row a sweep; one stroke weight
// for every line of the part, the window lit by its fill and grips alone.

import { loadGenome } from '@molstar/seqx-core/resolve';
import { beforeAll, describe, expect, it } from 'vitest';
import { LightTokens, type ViewLayout } from '../index.js';
import type { NavigatorFrame, NavigatorFrameInput } from '../navigator.js';
import { IdleAlpha, NavStrokePx } from '../window.js';
import { genomeSpec, toyGenome } from './genome-fixtures.js';
import { settle } from './helpers.js';
import { regionSpec, settleNavigator, withNavigator } from './navigator-fixtures.js';

const Section = '/views/v/genome';
const Size = { width: 1200, height: 300 };
const NavSize = { width: 1200, height: 28 };

beforeAll(async () => {
  await loadGenome();
});

const spec = () => withNavigator(genomeSpec(toyGenome(25, 0)), { of: Section, axes: 'positions', content: 'axis' });

/** The navigator over the toy genome, and the view's layout of the window `[a, b)` in block `k`'s own fractions. */
async function over() {
  const s = await settleNavigator(spec());
  const blocks = s.view.sections[0].blocks!;
  const layouts = new Map<string, ViewLayout>();
  const at = async (k: number, from: number, to: number) => {
    const key = `${k}:${from}:${to}`;
    if (!layouts.has(key)) {
      const v = await settle(spec());
      const o = blocks[k].offset!;
      const span = (blocks[k] as { end?: number }).end ?? blocks[k + 1].offset! - o;
      const f = await v.frame({ size: Size, viewport: { start: o + from * span, end: o + to * span } });
      layouts.set(key, f.layout);
    }
    return layouts.get(key)!;
  };
  const frame = async (
    k: number,
    from: number,
    to: number,
    lit?: NavigatorFrameInput['lit'],
  ): Promise<NavigatorFrame> => {
    const view = await at(k, from, to);
    await s.frames({ size: Size }, { size: NavSize });
    return s.nav.frame({ size: NavSize, view, moves: { [Section]: 1 }, ...(lit ? { lit } : {}) });
  };
  return { s, frame };
}

/** The window pass's rects drawn in `seq-grip` or the window's lit ink as grips: tabs, outside the window's own box. */
function gripRects(f: NavigatorFrame) {
  const r = f.lists.passes.window.rects;
  const out: { y: number }[] = [];
  const grip = LightTokens.colors['seq-grip'];
  for (let i = 0; i < r.n; i++) if (r.fill[i] === grip) out.push({ y: r.y[i] });
  return out;
}

describe('the navigator over a genome: two rows always', () => {
  it('chr6 and chr13 alike: the zoom strip, the same part height, the grips on the bottom row only', async () => {
    const { s, frame } = await over();
    const heights: number[] = [];
    for (const k of [5, 12]) {
      // half of the chromosome: a window with grips
      const f = await frame(k, 0.25, 0.5);
      heights.push(f.lists.height);
      expect(s.nav.zoomHeight).toBeGreaterThan(0);
      const strip = s.nav.pictureLayout!.sections[1];
      const grips = gripRects(f);
      expect(grips.length).toBe(2);
      for (const g of grips) expect(g.y).toBeGreaterThanOrEqual(strip.y - 1e-6);
      // the whole chromosome: still the strip (not the window's width), the window one quiet frame
      const whole = await frame(k, 0, 1);
      heights.push(whole.lists.height);
      expect(s.nav.zoomHeight).toBeGreaterThan(0);
      expect(gripRects(whole).length).toBe(0);
    }
    expect(new Set(heights).size).toBe(1);
  });

  it('every line of the part has one width: the window, the mark, the funnel, a lit window', async () => {
    const { s, frame } = await over();
    for (const lit of [undefined, 'body' as const]) {
      const f = await frame(5, 0.25, 0.5, lit && { section: Section, part: lit, axis: 'positions' });
      const { rects, paths } = f.lists.passes.window;
      const widths = new Set<number>();
      for (let i = 0; i < rects.n; i++)
        if (rects.stroke[i] !== 0 && rects.strokeWidth[i] > 0) widths.add(rects.strokeWidth[i]);
      for (let i = 0; i < paths.n; i++) widths.add(paths.strokeWidth[i]);
      expect([...widths]).toEqual([NavStrokePx]);
    }
  });
});

describe('a window at the fit is one quiet frame (F25)', () => {
  const R = '/views/v/region';
  const Box = { width: 600, height: 24 };

  it('the whole axis in view: no grips, no wash, the whole row a sweep; zoomed in, grips and the window a target', async () => {
    const s = await settleNavigator(regionSpec({ range: [0, 1000] }));
    const grip = LightTokens.colors['seq-grip'];
    const { nav: whole } = await s.frames({ size: { width: 600, height: 200 } }, { size: Box, moves: { [R]: 1 } });
    const rects = whole.lists.passes.window.rects;
    for (let i = 0; i < rects.n; i++) expect(rects.fill[i]).not.toBe(grip);
    for (const x of [1, 300, 599]) expect(s.nav.hit(x, 12)?.part).toBe('outside');
    // an outline only: its stroke the quiet alpha of the window's ink
    const quiet = ((LightTokens.colors['seq-window'] & 0xffffff00) | Math.round(IdleAlpha * 255)) >>> 0;
    expect(Array.from({ length: rects.n }, (_, i) => rects.stroke[i]).filter((c) => c !== 0)).toEqual([quiet]);

    const z = await settleNavigator(regionSpec({ range: [0, 1000] }));
    const { nav: part } = await z.frames(
      { size: { width: 600, height: 200 }, viewport: { start: 250, end: 500 } },
      { size: Box, moves: { [R]: 1 } },
    );
    const r2 = part.lists.passes.window.rects;
    expect(Array.from({ length: r2.n }, (_, i) => r2.fill[i]).filter((c) => c === grip).length).toBe(2);
    expect(z.nav.hit(225, 12)?.part).toBe('body');
  });
});
