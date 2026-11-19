/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ruler rail and the export's navigators (_spec/README.md "Render", Navigators; Export): the rail on a ruler's
// bottom border (its grab zone the ruler, an alignment's lower share), a slim row without a ruler (a hidden ruler
// too), idle with no target, the window's parts, a wrap over a cyclic origin; an export leaves navigators out unless
// asked, and with `include` stacks them around the view as outlines, titled.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, stripNavigator } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type FrameInput, hitView, LightTokens, railUnitAt, SlimRailPx } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { P53 } from './fixtures.js';
import { settle } from './helpers.js';
import { ideogramSpec } from './navigator-fixtures.js';
import { readRects } from './small-helpers.js';

const Size = { width: 800, height: 300 };
const Section = '/views/v/protein';

function p53(options: { ruler?: boolean; hiddenRuler?: boolean; cyclic?: boolean; nav?: boolean } = {}): Node {
  const b = createSeqSpec({ title: 'p53 rail' });
  const p = b.data.sequence('p53', { alphabet: 'protein', letters: P53, ...(options.cyclic ? { cyclic: true } : {}) });
  b.views.view('v', (v) => {
    const s = v.section('protein', { of: p });
    if (options.ruler !== false)
      s.track('ruler', options.hiddenRuler ? { hidden: true } : {}).layer('ruler', { type: 'ruler' });
    s.track('letters').layer('letters', { type: 'letters' });
    if (options.nav) v.navigator('strip', stripNavigator(s));
  });
  return b.build();
}

const rails = { [Section]: true };
const active = { [Section]: true };
const idle = { [Section]: false };

describe('the ruler rail', () => {
  it('is drawn just below the ruler’s bottom border; its grab zone the ruler’s height; a drag there is the window’s', async () => {
    const s = await settle(p53());
    const input: FrameInput = { size: Size, viewport: { start: 100, end: 200 }, rails, railsActive: active };
    const f = await s.frame(input);
    const sec = f.layout.sections[0];
    const rail = sec.rail!;
    const ruler = sec.tracks[0].box;
    const border = ruler.y + ruler.height;
    // in a row of its own below the border (the ruler's well grown by it), the next track below that row
    expect(rail.y).toBe(border + SlimRailPx / 2);
    expect([rail.top, rail.bottom]).toEqual([border + SlimRailPx - ruler.height, border + SlimRailPx]);
    expect(sec.tracks[0].well.y + sec.tracks[0].well.height).toBe(border + SlimRailPx);
    expect(sec.tracks[1].box.y).toBeGreaterThanOrEqual(border + SlimRailPx);
    expect([rail.start, rail.end, rail.total]).toEqual([100, 200, 393]);
    const mid = rail.x + (150 / 393) * rail.width;
    const h = hitView(f.layout, mid, rail.y)!;
    expect(h).toMatchObject({ target: 'window', part: 'body', navigator: 'rail', section: 'protein' });
    expect(railUnitAt(rail, mid)).toBeCloseTo(150, 9);
    expect(hitView(f.layout, rail.x + (100 / 393) * rail.width - 2, rail.y)?.part).toBe('start');
    expect(hitView(f.layout, rail.x + 5, rail.y)?.part).toBe('outside');
    // the rail in the window pass: its line, its segment, its grips
    expect(readRects(f.lists, 'window').length).toBe(4);
  });

  it('idle: a quiet hairline, no grips and no target (a press falls through to the ruler), with no relayout', async () => {
    const s = await settle(p53());
    const busy = await s.frame({ size: Size, rails, railsActive: active });
    const key = busy.key;
    // idle and live again: the window pass alone is drawn again, the layout kept
    const n = s.renderer.frame({ size: Size, rails, railsActive: idle });
    expect(n.rebuilt).toBe(false);
    expect(n.key).toBe(key);
    expect(n.layout.sections[0].rail!.active).toBe(false);
    const ruler = n.layout.sections[0].tracks[0].box;
    expect(hitView(n.layout, 400, ruler.y + 2)?.target).not.toBe('window');
    // at the fit the window is the whole axis: a hairline at a quiet opacity over the line, never a heavy rule
    const rects = readRects(n.lists, 'window');
    expect(rects.length).toBe(2);
    expect(rects.every((r) => r.h <= 1)).toBe(true);
    expect(rects[1].fill & 0xff).toBeLessThan(0xff);
    const again = s.renderer.frame({ size: Size, rails, railsActive: active });
    expect(again.rebuilt).toBe(false);
    expect(again.layout.sections[0].rail!.active).toBe(true);
    expect(readRects(again.lists, 'window').length).toBe(4);
  });

  it('collapsed breaks are a notch each on the rail', async () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    const pieces = b.data.ranges('pieces', { target: p, rows: { start: [0, 200], end: [100, 300] } });
    b.views.view('v', (v) => {
      const s = v.section('protein', { of: p, pieces });
      s.track('ruler').layer('ruler', { type: 'ruler' });
      s.track('letters').layer('letters', { type: 'letters' });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: Size, viewport: { start: 0, end: 100 }, rails, railsActive: active });
    const rail = f.layout.sections[0].rail!;
    expect(rail.total).toBe(200);
    expect(rail.breaks).toEqual([100]);
    // the notch: 1 px wide, 4 px tall, at the break's x on the rail
    const x = rail.x + (100 / 200) * rail.width - 0.5;
    const notch = readRects(f.lists, 'window').filter((r) => r.w === 1 && r.h === 4);
    expect(notch.map((r) => r.x)).toEqual([x]);
  });

  it('by touch its grips hit wider', async () => {
    const s = await settle(p53());
    const f = await s.frame({ size: Size, viewport: { start: 100, end: 200 }, rails, railsActive: active });
    const rail = f.layout.sections[0].rail!;
    const x = rail.x + (100 / 393) * rail.width - 10;
    expect(hitView(f.layout, x, rail.y)?.part).toBe('outside');
    expect(hitView(f.layout, x, rail.y, true)?.part).toBe('start');
  });

  it('a section with no ruler gets a slim row under its tracks; a hidden ruler falls back to it', async () => {
    const plain = await settle(p53({ ruler: false }));
    // a frame is valid until the next: what is kept is read now
    const before = (await plain.frame({ size: Size })).layout.sections[0];
    const [y0, h0] = [before.tracks[0].box.y, before.height];
    const withRail = await plain.frame({ size: Size, rails, railsActive: active });
    const s = withRail.layout.sections[0];
    expect(s.tracks[0].box.y).toBe(y0);
    expect(s.height - h0).toBe(12);
    const rail = s.rail!;
    expect(rail.slim).toBe(true);
    expect(rail.bottom - rail.top).toBe(12);
    expect(rail.top).toBe(s.y + h0);
    const hidden = await settle(p53({ hiddenRuler: true }));
    expect((await hidden.frame({ size: Size, rails, railsActive: active })).layout.sections[0].rail!.slim).toBe(true);
  });

  it('without rails asked, the layout is as before', async () => {
    const s = await settle(p53());
    const f = await s.frame({ size: Size });
    expect(f.layout.sections[0].rail).toBeUndefined();
    expect(readRects(f.lists, 'window').length).toBe(0);
  });

  it('a window over a cyclic origin is drawn in two pieces, a grip at each true end; the axis’s ends are body', async () => {
    const s = await settle(p53({ cyclic: true }));
    const f = await s.frame({ size: Size, viewport: { start: 380, end: 410 }, rails, railsActive: active });
    const rail = f.layout.sections[0].rail!;
    expect(rail.end).toBeGreaterThan(rail.total);
    const segs = readRects(f.lists, 'window').filter((r) => r.h === 2);
    expect(segs.length).toBe(2);
    const k = rail.width / rail.total;
    const start = rail.x + 380 * k;
    const end = rail.x + 17 * k;
    const grips = readRects(f.lists, 'window').filter((r) => r.w === 2);
    const xs = grips.map((r) => r.x).sort((a, b) => a - b);
    expect(xs.length).toBe(2);
    expect(xs[0]).toBeCloseTo(end, 3);
    expect(xs[1]).toBeCloseTo(start - 2, 3);
    const y = rail.y;
    expect(hitView(f.layout, rail.x + rail.width - 1, y)?.part).toBe('body');
    expect(hitView(f.layout, rail.x + 1, y)?.part).toBe('body');
    expect(hitView(f.layout, start - 1, y)?.part).toBe('start');
    expect(hitView(f.layout, end + 1, y)?.part).toBe('end');
  });

  it('never covers the ruler’s ticks or labels, at rest and grown', async () => {
    const s = await settle(p53());
    const input: FrameInput = { size: Size, viewport: { start: 100, end: 200 }, rails, railsActive: active };
    const f = await s.frame(input);
    const rail = f.layout.sections[0].rail!;
    const ruler = f.layout.sections[0].tracks[0].box;
    const ticks = readRects(f.lists, 'data').filter((r) => r.y >= ruler.y && r.y + r.h <= ruler.y + ruler.height);
    expect(ticks.length).toBeGreaterThan(0);
    const lowest = Math.max(...ticks.map((r) => r.y + r.h));
    const h = { ...hitView(f.layout, 400, rail.y)! };
    for (const frame of [f, s.renderer.frame({ ...input, emphasis: { hover: h } })]) {
      const marks = readRects(frame.lists, 'window');
      expect(marks.length).toBe(4);
      for (const r of marks) expect(r.y).toBeGreaterThanOrEqual(lowest);
    }
  });

  it('grows while the pointer is over its grab zone, and while its drag runs', async () => {
    const s = await settle(p53());
    const input: FrameInput = { size: Size, viewport: { start: 100, end: 200 }, rails, railsActive: active };
    const f = await s.frame(input);
    const h = { ...hitView(f.layout, 400, f.layout.sections[0].rail!.y)! };
    const grown = s.renderer.frame({ ...input, emphasis: { hover: h } });
    expect(readRects(grown.lists, 'window').some((r) => r.h === 6)).toBe(true);
    const dragged = s.renderer.frame({ ...input, emphasis: { drag: h } });
    expect(readRects(dragged.lists, 'window').some((r) => r.h === 6)).toBe(true);
  });
});

describe('navigators in an export', () => {
  it('are left out by default; included, stacked around the view as outlines, each a titled group', async () => {
    const spec = ideogramSpec();
    const omit = await renderSvg(spec, { size: { width: 600, height: 120 } });
    expect(omit.svg).not.toContain('data-navigator');
    const include = await renderSvg(spec, { size: { width: 600, height: 120 }, navigators: 'include' });
    expect(include.svg).toContain('data-navigator="chr"');
    expect(include.svg).toContain('<title>Window: chr:301–360 of chr</title>');
    // the view's 120, the ideogram's 24 and its zoom strip's 16 (the region is 36 px of 600: under 48)
    expect(include.svg).toMatch(/^<svg xmlns="http:\/\/www.w3.org\/2000\/svg" width="600" height="160"/);
    expect(include.settled).toBe(true);
  });

  it('the view’s `export.navigators` says include, the call’s option overrides it', async () => {
    const b = createSeqSpec();
    const p = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    b.views.view('v', { export: { navigators: 'include' } }, (v) => {
      const s = v.section('protein', { of: p });
      s.track('letters').layer('letters', { type: 'letters' });
      v.navigator('strip', stripNavigator(s));
    });
    const spec = b.build();
    expect((await renderSvg(spec, { size: { width: 600, height: 80 } })).svg).toContain('data-navigator="strip"');
    expect((await renderSvg(spec, { size: { width: 600, height: 80 }, navigators: 'omit' })).svg).not.toContain(
      'data-navigator',
    );
  });

  it('included, a section no navigator binds draws its rail as a plain window: no grips, no idle state', async () => {
    const r = await renderSvg(p53(), { size: { width: 600, height: 120 }, navigators: 'include' });
    const rects = readRects(r.lists, 'window');
    expect(rects.length).toBe(2);
    expect(rects[1]).toMatchObject({ h: 2, fill: LightTokens.colors['seq-window'] });
    // omitted (the default), no rail at all
    expect(readRects((await renderSvg(p53(), { size: { width: 600, height: 120 } })).lists, 'window').length).toBe(0);
  });
});
