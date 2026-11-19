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
import { hitView, type FrameInput, railUnitAt, renderSvg } from '../index.js';
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
  it('is drawn on the ruler’s bottom border; the ruler row is its grab zone; a drag there is the window’s', async () => {
    const s = await settle(p53());
    const input: FrameInput = { size: Size, viewport: { start: 100, end: 200 }, rails, railsActive: active };
    const f = await s.frame(input);
    const sec = f.layout.sections[0];
    const rail = sec.rail!;
    const ruler = sec.tracks[0].box;
    expect(rail.y).toBe(ruler.y + ruler.height);
    expect([rail.top, rail.bottom]).toEqual([ruler.y, ruler.y + ruler.height]);
    expect([rail.start, rail.end, rail.total]).toEqual([100, 200, 393]);
    const mid = rail.x + (150 / 393) * rail.width;
    const h = hitView(f.layout, mid, ruler.y + 2)!;
    expect(h).toMatchObject({ target: 'window', part: 'body', navigator: 'rail', section: 'protein' });
    expect(railUnitAt(rail, mid)).toBeCloseTo(150, 9);
    expect(hitView(f.layout, rail.x + (100 / 393) * rail.width - 2, ruler.y + 2)?.part).toBe('start');
    expect(hitView(f.layout, rail.x + 5, ruler.y + 2)?.part).toBe('outside');
    // the rail in the window pass: its line, its segment, its grips
    expect(readRects(f.lists, 'window').length).toBe(4);
  });

  it('idle: a quiet segment, no grips and no target (a press falls through to the ruler), with no relayout', async () => {
    const s = await settle(p53());
    const busy = await s.frame({ size: Size, rails, railsActive: active });
    const n = s.renderer.frame({ size: Size, rails, railsActive: idle });
    expect(n.layout.sections[0]).toBe(busy.layout.sections[0]);
    const ruler = n.layout.sections[0].tracks[0].box;
    expect(hitView(n.layout, 400, ruler.y + 2)?.target).not.toBe('window');
    expect(readRects(n.lists, 'window').length).toBe(2);
  });

  it('a section with no ruler gets a slim row; a hidden ruler falls back to it', async () => {
    const plain = await settle(p53({ ruler: false }));
    // a frame is valid until the next: what is kept is read now
    const y0 = (await plain.frame({ size: Size })).layout.sections[0].tracks[0].box.y;
    const withRail = await plain.frame({ size: Size, rails, railsActive: active });
    expect(withRail.layout.sections[0].tracks[0].box.y - y0).toBe(12);
    const rail = withRail.layout.sections[0].rail!;
    expect(rail.slim).toBe(true);
    expect(rail.bottom - rail.top).toBe(12);
    const hidden = await settle(p53({ hiddenRuler: true }));
    expect((await hidden.frame({ size: Size, rails, railsActive: active })).layout.sections[0].rail!.slim).toBe(true);
  });

  it('without rails asked, the layout is as before', async () => {
    const s = await settle(p53());
    const f = await s.frame({ size: Size });
    expect(f.layout.sections[0].rail).toBeUndefined();
    expect(readRects(f.lists, 'window').length).toBe(0);
  });

  it('a window over a cyclic origin is drawn in two pieces', async () => {
    const s = await settle(p53({ cyclic: true }));
    const f = await s.frame({ size: Size, viewport: { start: 380, end: 410 }, rails, railsActive: active });
    const rail = f.layout.sections[0].rail!;
    expect(rail.end).toBeGreaterThan(rail.total);
    const segs = readRects(f.lists, 'window').filter((r) => r.h === 2);
    expect(segs.length).toBe(2);
  });

  it('grows while the pointer is over its grab zone, and while its drag runs', async () => {
    const s = await settle(p53());
    const input: FrameInput = { size: Size, viewport: { start: 100, end: 200 }, rails, railsActive: active };
    const f = await s.frame(input);
    const ruler = f.layout.sections[0].tracks[0].box;
    const h = { ...hitView(f.layout, 400, ruler.y + 2)! };
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
    expect(include.svg).toMatch(/^<svg xmlns="http:\/\/www.w3.org\/2000\/svg" width="600" height="144"/);
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

  it('included, a section no navigator binds draws its rail idle', async () => {
    const r = await renderSvg(p53(), { size: { width: 600, height: 120 }, navigators: 'include' });
    const pass = r.lists.passes.window;
    expect(pass.rects.n).toBe(2);
  });
});
