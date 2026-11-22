/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The zoom strip (_spec/README.md "Render", Navigators, The zoom strip): a region a sliver of its chromosome is drawn as
// a mark on the context and a zoom strip under it over the region alone. Across a cyclic origin (the shorter arc, the
// strip's viewport past the origin); a `heightPx` small or 0; the gap a finger's; several bound sections (one on the
// navigator's own axis, two far apart); a side part; a strip's axis and tracks as the context; the mark centred
// whatever the floor; no hysteresis (the strip is the span's); the mark kept
// while the bound section is pending; the tokens kept over a picture made again; a reversed section.

import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { DarkTokens, DefaultRenderConfig, NavWindowFlags, type ViewLayout } from '../index.js';
import { type NavigatorRendererOptions } from '../navigator.js';
import { IdleAlpha } from '../window.js';
import { regionSpec, settleNavigator } from './navigator-fixtures.js';

const R = '/views/v/region';
const Box = { width: 600, height: 24 };

/** The region's navigator settled at `nav` (its box), the view at `width`; the strip's section and the frame. */
async function zoomed(
  spec: Parameters<typeof regionSpec>[0] = {},
  nav: { width: number; height: number } = Box,
  options: Partial<NavigatorRendererOptions> = {},
  config: Parameters<typeof settleNavigator>[3] = {},
) {
  const s = await settleNavigator(regionSpec(spec), 0, [], config, () => options);
  const { nav: f, view } = await s.frames(
    { size: { width: nav.width, height: 200 } },
    { size: nav, moves: { [R]: NavWindowFlags.x } },
  );
  const layout = s.nav.pictureLayout!;
  return { s, f, view, ctx: layout.sections[0], strip: s.nav.zoomHeight ? layout.sections[1] : undefined };
}

/** The window pass's outline over the context's height (the mark's pieces). */
function marks(f: { lists: { passes: { window: { rects: Record<string, any> } } } }, top: number, h: number) {
  const r = f.lists.passes.window.rects;
  const out: { x: number; w: number }[] = [];
  for (let i = 0; i < r.n; i++)
    if (r.y[i] === top && r.h[i] === h && r.fill[i] === 0) out.push({ x: r.x[i], w: r.w[i] });
  return out;
}

describe('the zoom strip', () => {
  it('a region across a cyclic origin: the shorter arc, the strip past the origin, the mark in two pieces', async () => {
    const { s, f, strip } = await zoomed({ length: 60_000, cyclic: true, range: [59_000, 1000] });
    expect(s.nav.zoomHeight).toBe(16);
    // the strip's axis is the arc, 59,000 on past the origin to 1,000
    expect(strip!.axis!.u0).toBeCloseTo(59_000, 6);
    expect(strip!.axis!.u1).toBeCloseTo(61_000, 6);
    expect(s.nav.markText()).toBe('chr:59,001–1,000');
    // the region's whole view is the strip's whole width, one piece
    const w = f.windows[0];
    const box = strip!.placement.box;
    expect(w.wrap).toBe(false);
    expect(w.x0).toBeCloseTo(box.x, 3);
    expect(w.x1).toBeCloseTo(box.x + box.width, 3);
    expect(w.lo).toBeCloseTo(box.x, 3);
    expect(w.hi).toBeCloseTo(box.x + box.width, 3);
    // a drag maps across the origin in the region's own units
    expect(s.nav.sectionUnitAt(R, box.x + box.width / 2)).toBeCloseTo(1000, 0);
    expect(s.nav.sectionUnitAt(R, box.x + box.width / 4)).toBeCloseTo(500, 0);
    // the mark at both ends of the context
    const m = marks(f, 0, 24);
    expect(m.length).toBe(2);
    expect(Math.min(...m.map((r) => r.x))).toBeCloseTo(0, 3);
    expect(Math.max(...m.map((r) => r.x + r.w))).toBeCloseTo(Box.width, 3);
    expect(f.lists.passes.window.paths.n).toBe(2);
  });

  it.each([
    [0, 0],
    [4, 4],
    [16, 16],
  ])(
    '`heightPx: %i`: the strip never lays out taller (%i more); 0 is off, the window back on the context',
    async (h) => {
      const { s, f, strip } = await zoomed({}, Box, { zoomStrip: { minPx: 48, heightPx: h } });
      expect(s.nav.zoomHeight).toBe(h);
      expect(f.lists.height).toBe(Box.height + h);
      const w = f.windows[0];
      if (h === 0) {
        expect(strip).toBeUndefined();
        expect([w.ly, w.hy]).toEqual([0, Box.height]);
        expect(f.lists.passes.window.paths.n).toBe(0);
      } else {
        expect(strip!.y + strip!.height).toBeLessThanOrEqual(Box.height + h + 1e-6);
        expect(strip!.y).toBeGreaterThanOrEqual(Box.height);
        expect(w.ly).toBeCloseTo(strip!.y, 6);
      }
    },
  );

  it('the gap is the strip’s: a finger 2 px above the strip lands on its window', async () => {
    const { s, ctx, strip } = await zoomed();
    const x = Box.width / 2;
    // the region is at its fit on the strip: no target, the row a sweep (`outside`)
    expect(s.nav.hit(x, strip!.y - 2, true)?.part).toBe('outside');
    expect(s.nav.hit(x, strip!.y - 2, true)?.context).toBeUndefined();
    expect(s.nav.hit(x, ctx.y + ctx.height - 1, true)?.context).toBe('outside');
  });

  it('a bound section on the navigator’s whole axis: no strip, its window right; two far regions keep the hull', async () => {
    // the navigator on the chromosome's section's own axis, and a short sequence mapped onto it beside
    const b = createSeqSpec({ title: 'own axis' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', letters: 'ACGT'.repeat(250) });
    const short = b.data.sequence('short', { alphabet: 'dna', letters: 'ACGT'.repeat(5) });
    b.data.map('short-on-chr', { from: short, to: chr, method: 'offset', by: 300 });
    b.views.view('v', (v) => {
      const c = v.section('chr', { of: chr });
      c.track('letters').layer('letters', { type: 'letters' });
      const s = v.section('short', { of: short });
      s.track('letters').layer('letters', { type: 'letters' });
      v.navigator('chr', stripNavigator([c.address, s.address]));
    });
    const own = await settleNavigator(b.build());
    const { nav: f } = await own.frames({ size: { width: 600, height: 200 } }, { size: Box });
    expect(own.nav.zoomHeight).toBe(0);
    const axis = own.nav.pictureSection!.axis!;
    const w = f.windows[0];
    expect(w.x0).toBeCloseTo(axis.xOfUnit(0), 3);
    expect(w.x1).toBeCloseTo(axis.xOfUnit(1000), 3);
    const whole = await zoomed({ others: ['whole'] });
    expect(whole.s.nav.zoomHeight).toBe(0);
    // the hull of two far regions is wide: one strip (the spec's known case)
    const far = await zoomed({ range: [100, 110], others: [[900, 910]] });
    expect(far.s.nav.zoomHeight).toBe(0);
  });

  it('a side part gets no zoom strip', async () => {
    const { s } = await zoomed({ place: 'right', content: 'axis' });
    expect(s.nav.zoomHeight).toBe(0);
  });

  it.each(['axis', 'tracks'] as const)('a strip of its %s: the zoom strip at the box’s bottom', async (content) => {
    const box = { width: 600, height: 48 };
    const { s, f, ctx, strip } = await zoomed({ content }, box);
    expect(s.nav.zoomHeight).toBe(16);
    expect(ctx.y + ctx.height).toBeLessThanOrEqual(box.height + 1e-6);
    expect(strip!.y + strip!.height).toBeCloseTo(box.height + 16, 6);
    expect(f.lists.height).toBe(box.height + 16);
    // the connectors from the context's foot to the strip's top
    const paths = f.lists.passes.window.paths;
    expect(paths.n).toBe(2);
    for (let i = 0; i < 2; i++) {
      const at = paths.cmdStart[i];
      expect(paths.points[2 * at + 3]).toBeCloseTo(strip!.y, 6);
    }
  });

  it('the mark is centred on the truth, whatever `window.floor`; the truth inside it', async () => {
    const thumb = { navigator: { window: { ...DefaultRenderConfig.navigator.window, floor: 'thumb' as const } } };
    const { f, ctx } = await zoomed({ range: [100, 110] }, Box, {}, thumb);
    const axis = ctx.axis!;
    const [m] = marks(f, 0, 24);
    const a = axis.xOfUnit(100);
    const b = axis.xOfUnit(110);
    expect(m.w).toBeCloseTo(16, 6);
    expect(m.x + m.w / 2).toBeCloseTo((a + b) / 2, 6);
    expect(m.x).toBeLessThanOrEqual(a);
    expect(m.x + m.w).toBeGreaterThanOrEqual(b);
  });

  it('one system: the strip is the reachable span’s, the same at a width on either side of `minPx` (no hysteresis)', async () => {
    // the region is 60 of 1,000: 48 px at a width of 800
    const { s, view } = await zoomed({ range: [300, 360] }, { width: 780, height: 24 });
    const at = (width: number) => {
      s.nav.frame({ size: { width, height: 24 }, view: view.layout });
      return s.nav.zoomHeight > 0;
    };
    expect([780, 790, 799].map(at)).toEqual([true, true, true]);
    expect([801, 820, 940, 920, 810].map(at)).toEqual([false, false, false, false, false]);
    expect(at(790)).toBe(true);
  });

  it('the mark and the connectors kept while the bound section is pending (its last span)', async () => {
    const { s, view } = await zoomed();
    const pending: ViewLayout = {
      ...view.layout,
      sections: view.layout.sections.map((x) => ({ ...x, axis: undefined })),
    } as ViewLayout;
    const g = s.nav.frame({ size: Box, view: pending });
    expect(s.nav.zoomHeight).toBe(16);
    expect(g.windows[0].shown).toBe(false);
    expect(marks(g, 0, 24).length).toBe(1);
    expect(g.lists.passes.window.paths.n).toBe(2);
    expect(s.nav.markText()).toBe('chr:301–330');
  });

  it('the tokens set survive a picture made again', async () => {
    const { s, view } = await zoomed();
    s.nav.setTokens(DarkTokens);
    // another height: the ideogram fills it, made again
    const g = s.nav.frame({ size: { width: 600, height: 30 }, view: view.layout });
    expect(g.rebuilt).toBe(true);
    const p = g.lists.passes.window.paths;
    const quiet = ((DarkTokens.colors['seq-window'] & 0xffffff00) | Math.round(IdleAlpha * 255)) >>> 0;
    for (let i = 0; i < p.n; i++) expect(p.stroke[i]).toBe(quiet);
    expect(g.lists.passes.background.rects.fill[0]).toBe(DarkTokens.colors['msui-bg']);
  });

  it('a reversed section: the strip runs reversed, a drag maps back into the region’s units', async () => {
    const { s, f, strip } = await zoomed({ reversed: true });
    expect(s.nav.zoomHeight).toBe(16);
    expect(strip!.axis!.dir).toBe(-1);
    const w = f.windows[0];
    const box = strip!.placement.box;
    expect(w.x0).toBeCloseTo(box.x, 3);
    expect(w.x1).toBeCloseTo(box.x + box.width, 3);
    expect(s.nav.sectionUnitAt(R, box.x + box.width / 3)).toBeCloseTo(10, 0);
  });
});
