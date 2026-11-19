/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The section layout interface: groups on the axis (cells skip the gaps, spans bridge them), the tracks layout's
// placement as the inverse of its hits, the translate and polar warps, placed glyph runs in every backend, the window
// mark and its parts, the navigation's arithmetic, a section not laid out again while another moves, a layer with a
// second axis drawn in the emphasis passes, reused buffers against fresh ones.

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import {
  AxisTransform,
  Cut,
  DefaultRenderConfig,
  DrawLists,
  drawListsToSvg,
  drawWindowArgs,
  type Emphasis,
  type LayerType,
  LightTokens,
  type NavContext,
  PlacePoint,
  PolarWarp,
  recordDrawLists,
  StripPlacement,
  TracksNav,
  TranslateWarp,
  windowPart,
} from '../index.js';
import { BuiltInLayerTypes } from '../tracks.js';
import { P53, p53TracksSpec } from './fixtures.js';
import { settle } from './helpers.js';
import { rng } from './small-helpers.js';

const Size = { width: 640, height: 200 };

describe('groups on the axis', () => {
  // 25 positions, groups of 10 with a 6 px gap, 10 px a cell: 250 px of cells and two gaps
  const a = new AxisTransform('s', [0, 25], 0, 262, 1, undefined, { groups: { size: 10, gap: 6 } });

  it('cells skip the gaps: a group’s first cell starts past the gap before it', () => {
    expect(a.scale).toBeCloseTo(10);
    expect(a.cell(9)).toBe(true);
    expect([a.sx0, a.sx1]).toEqual([90, 100]);
    expect(a.cell(10)).toBe(true);
    expect([a.sx0, a.sx1]).toEqual([106, 116]);
    expect(a.cell(24)).toBe(true);
    expect(a.sx1).toBeCloseTo(262);
  });

  it('a span bridges the gaps inside it and stops at its own ends', () => {
    a.piecesIn(5, 20);
    expect(a.spanIn(0, 5, 20)).toBe(true);
    // from position 5's left edge to position 19's right edge: the gap at 10 inside, the one at 20 not
    expect([a.sx0, a.sx1]).toEqual([50, 206]);
    expect(a.scut).toBe(0);
  });

  it('x back to a unit: a cell’s unit, a gap the next group’s first', () => {
    for (let p = 0; p < 25; p++) {
      a.cell(p);
      expect(Math.floor(a.unitAt((a.sx0 + a.sx1) / 2))).toBe(p);
    }
    expect(a.inGap(102)).toBe(true);
    expect(a.unitAt(102)).toBe(10);
    expect(a.inGap(95)).toBe(false);
  });

  it('a line of its window: a range over the line’s end is cut there', () => {
    const line = new AxisTransform('s', [0, 40], 0, 206, 1, { start: 0, end: 20 }, { groups: { size: 10, gap: 6 } });
    line.piecesIn(15, 30);
    expect(line.spanIn(0, 15, 30)).toBe(true);
    expect(line.scut & Cut.right).toBe(Cut.right);
    expect(line.sx1).toBe(206);
  });
});

describe('the tracks layout on the interface', () => {
  it('placement is the inverse of hit-testing: every position in view placed once, `at` its cell’s center gives it back', async () => {
    const s = await settle(p53TracksSpec());
    const r = rng(7);
    for (let n = 0; n < 8; n++) {
      const start = Math.floor(r() * 300);
      const end = start + 20 + Math.floor(r() * 70);
      const f = await s.frame({ size: Size, viewport: { start, end } });
      const sec = f.layout.sections[0];
      const p = sec.placement;
      const pt = new PlacePoint();
      let last = -Infinity;
      for (let u = start; u < end; u++) {
        expect(p.place(u + 0.5, pt)).toBe(true);
        const x = pt.x;
        const y = sec.tracks[0].box.y + 1;
        expect(p.at([x, y], pt)).toBe(true);
        expect(pt.kind).toBe('data');
        expect(Math.floor(pt.offset)).toBe(u);
        const o = p.offsetAt(x, y);
        expect(o).toBeGreaterThan(last);
        last = o;
      }
      // the label column is the gutter
      expect(p.at([20, sec.tracks[0].box.y + 1], pt) && pt.kind).toBe('gutter');
    }
  });

  it('a section whose key holds is not laid out again while another moves', async () => {
    const b = createSeqSpec({ title: 'two' });
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    b.views.view('v', (v) => {
      v.section('a', { of: p53 }).track('l').layer('letters', { type: 'letters' });
      v.section('b', { of: p53 }).track('l').layer('letters', { type: 'letters' });
    });
    const s = await settle(b.build());
    const vp = (a: number, b: number) => ({
      '/views/v/a': { start: a, end: a + 40 },
      '/views/v/b': { start: b, end: b + 40 },
    });
    const f1 = await s.frame({ size: Size, viewports: vp(0, 100) });
    const [a1, b1] = f1.layout.sections;
    const lay = a1.lays[0][0].lay;
    const f2 = s.renderer.frame({ size: Size, viewports: vp(0, 120) });
    expect(f2.rebuilt).toBe(true);
    const [a2, b2] = f2.layout.sections;
    expect(a2).toBe(a1);
    expect(a2.lays[0][0].lay).toBe(lay);
    expect(b2).not.toBe(b1);
    // the reused section is drawn as a fresh layout draws it
    const fresh = await settle(b.build());
    const g = await fresh.frame({ size: Size, viewports: vp(0, 120) });
    expect(recordDrawLists(f2.lists)).toBe(recordDrawLists(g.lists));
  });

  it('reused buffers draw as fresh ones (a large window warms them, a small one is drawn)', async () => {
    const warm = await settle(p53TracksSpec());
    await warm.frame({ size: Size, viewport: { start: 0, end: 393 } });
    const a = await warm.frame({ size: Size, viewport: { start: 30, end: 60 } });
    const fresh = await settle(p53TracksSpec());
    const b = await fresh.frame({ size: Size, viewport: { start: 30, end: 60 } });
    expect(recordDrawLists(a.lists)).toBe(recordDrawLists(b.lists));
  });

  it('a layer with a second axis is drawn in every emphasis pass with the emphasis (the y path)', async () => {
    const seen: string[] = [];
    const ys: LayerType<null, null> = {
      name: 'ys',
      axes: 2,
      prepare: () => null,
      layout: () => null,
      draw(_lay, pass, out, emphasis) {
        // the rows of its y system the selection names, painted by the layer: a rect a row
        const part = emphasis?.selection?.parts.find((p) => p.system === 'rows');
        if (part) {
          seen.push(pass);
          out.fill(0, part.ranges[0] * 10, 10, 10, 0xff0000ff);
        }
      },
      hit: () => undefined,
    };
    const b = createSeqSpec({ title: 'y' });
    const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
    b.views.view('v', (v) => v.section('a', { of: p53 }).track('l').layer('l', { type: 'letters' }));
    const spec: Node = b.build();
    const s = await settle(spec);
    const { Renderer, bundledTextMetrics, renderConfig } = await import('../index.js');
    const r = new Renderer({
      model: s.model,
      view: s.view,
      data: s.store,
      tokens: LightTokens,
      text: bundledTextMetrics(),
      config: renderConfig({}),
      types: { ...BuiltInLayerTypes, letters: ys as LayerType<any, any> },
    });
    const emphasis: Emphasis = { selection: rangesLocation('rows', [2, 3]) };
    r.frame({ size: Size, emphasis });
    // the window pass is the layouts' own: no layer is drawn there
    expect(seen).toEqual(['data', 'highlight', 'selection', 'hover', 'cursor', 'labels']);
  });
});

describe('a placement of translated strips', () => {
  it('`at(place(o))` gives `o` back: on the page, the middle of the strip’s first band with height', () => {
    const p = new StripPlacement();
    p.box = { x: 0, y: 0, width: 400, height: 200 };
    const s = p.strips;
    s.reset(2);
    // two lines of 20 units, the second 60 px down; each line's band 0 collapsed, band 1 at 4 px, 10 high
    for (let l = 0; l < 2; l++) {
      const a = new AxisTransform('s', [0, 40], 50, 200, 1, { start: 20 * l, end: 20 * l + 20 });
      const i = s.add(a, 0, 0, 60 * l, 20 * l);
      s.bandY[2 * i + 1] = 4;
      s.bandH[2 * i + 1] = 10;
      s.height[i] = 14;
    }
    const pt = new PlacePoint();
    for (let o = 0.5; o < 40; o += 1) {
      expect(p.place(o, pt)).toBe(true);
      expect(pt.band).toBe(1);
      expect(p.at([pt.x, pt.y], pt)).toBe(true);
      expect(pt.kind).toBe('data');
      expect(pt.offset).toBeCloseTo(o);
    }
  });
});

describe('warps', () => {
  function drawSome(out: ReturnType<DrawLists['writer']>, dx = 0, dy = 0) {
    out.clip(10 + dx, 5 + dy, 100, 20);
    out.rect(12 + dx, 6 + dy, 30, 10, 0x112233ff, 0x445566ff, 1, 3, 3);
    out.beginGlyphs('10px mono', 15 + dy, 0x000000ff);
    out.letter(65, 20 + dx);
    out.letter(66, 30 + dx);
    out.endGlyphs();
    out.text('label', 40 + dx, 15 + dy, 0, '10px sans', 0x000000ff, 20, 30);
    out.path(
      (p) => {
        p.command(0);
        p.point(10 + dx, 10 + dy);
        p.command(1);
        p.point(50 + dx, 12 + dy);
      },
      0,
      0x000000ff,
      1,
    );
    out.unclip();
  }

  it('a translation moves what a strip drew as drawing it there would', () => {
    const a = new DrawLists();
    const w = a.writer('data');
    const t = new TranslateWarp();
    t.open(w);
    drawSome(w);
    t.close(w, Float64Array.of(7), Float64Array.of(120), 0);
    const b = new DrawLists();
    drawSome(b.writer('data'), 7, 120);
    expect(recordDrawLists(a, ['data'])).toBe(recordDrawLists(b, ['data']));
  });

  it('the polar warp bends rects into sectors, places letters, refuses images, and its inverse brings a point back', () => {
    const lists = new DrawLists();
    const out = lists.writer('data');
    const p = new PolarWarp();
    p.cx = 200;
    p.cy = 200;
    p.r0 = 150;
    p.x0 = 0;
    p.width = 360;
    p.tangent = true;
    const w = p.writer(out);
    w.rect(0, 0, 90, 20, 0x112233ff);
    w.beginGlyphs('10px mono', 10, 0x000000ff);
    w.letter(65, 90);
    w.endGlyphs();
    w.image(0, 0, 10, 10, { width: 1, height: 1, data: new Uint8Array(4) });
    p.close(out);
    expect(p.refused).toBe(1);
    const paths = out.pass.paths;
    expect(paths.n).toBe(1);
    // every point of the sector is on a radius between 130 and 150, and within the quarter it spans
    for (let i = 0; i < paths.pointCount; i++) {
      const x = paths.points[2 * i] - 200;
      const y = paths.points[2 * i + 1] - 200;
      const r = Math.hypot(x, y);
      expect(r).toBeGreaterThan(129.99);
      expect(r).toBeLessThan(150.01);
      expect(x).toBeGreaterThan(-0.01);
      expect(y).toBeLessThan(0.01);
    }
    // the arc is cut every few px: a quarter of a 150 px ring is many segments
    expect(paths.pointCount).toBeGreaterThan(40);
    const g = out.pass.glyphs;
    expect(g.placed[0]).toBe(1);
    expect(g.x[0]).toBeCloseTo(200 + 140);
    expect(g.ly[0]).toBeCloseTo(200);
    expect(g.angle[0]).toBeCloseTo(Math.PI / 2);
    for (const [x, y] of [
      [10, 5],
      [200, 30],
      [359, 0],
    ]) {
      p.map(x, y);
      expect(p.unmap(p.px, p.py)).toBe(true);
      expect(p.px).toBeCloseTo(x);
      expect(p.py).toBeCloseTo(y);
    }
  });

  it('a polar strip draws into a scratch pass named as the pass it warps (the emphasis reads the name)', () => {
    const lists = new DrawLists();
    const p = new PolarWarp();
    for (const name of ['data', 'highlight', 'selection', 'hover', 'cursor', 'labels'] as const)
      expect(p.writer(lists.writer(name)).pass.name).toBe(name);
  });

  it('a mixed run: letters added before or after it turns placed sit on its baseline, unturned, unscaled', () => {
    const lists = new DrawLists();
    const out = lists.writer('data');
    // stale numbers in the letters' slots, from an earlier placed run
    out.beginGlyphs('10px mono', 0, 0x000000ff);
    const g = out.args;
    for (let k = 0; k < 3; k++) {
      g[0] = k;
      g[1] = 99;
      g[2] = 1;
      g[3] = 3;
      g[4] = 0;
      out.placedLetterArgs(65);
    }
    out.endGlyphs();
    lists.reset();
    out.beginGlyphs('10px mono', 20, 0x000000ff);
    out.letter(65, 1);
    g[0] = 2;
    g[1] = 30;
    g[2] = 0;
    g[3] = 1;
    out.placedLetterArgs(66);
    out.letter(67, 3);
    out.endGlyphs();
    expect(recordDrawLists(lists, ['data'])).toContain(
      'glyphs "10px mono" placed #000000ff A 1 20 0 1, B 2 30 0 1, C 3 20 0 1',
    );
  });

  it('the polar warp clips a letter by its width, not its center alone', () => {
    const lists = new DrawLists();
    const out = lists.writer('data');
    const p = new PolarWarp();
    p.cx = 200;
    p.cy = 200;
    p.r0 = 150;
    p.width = 360;
    const w = p.writer(out);
    w.clip(10, 0, 100, 20);
    w.beginGlyphs('10px mono', 10, 0x000000ff);
    // 10 px mono: 3 px either side of its center; the first overlaps the clip, the second does not
    w.letter(65, 112);
    w.letter(66, 120);
    w.endGlyphs();
    w.unclip();
    p.close(out);
    expect(out.pass.glyphs.letters).toBe(1);
    expect(out.pass.glyphs.codes[0]).toBe(65);
  });

  it('a placed glyph run in the record and the SVG: each letter its place, angle and vertical scale', () => {
    const lists = new DrawLists();
    const out = lists.writer('data');
    out.beginGlyphs('10px mono', 0, 0x000000ff);
    const g = out.args;
    g[0] = 10;
    g[1] = 20;
    g[2] = Math.PI / 2;
    g[3] = 2;
    g[4] = 0;
    out.placedLetterArgs(77);
    out.endGlyphs();
    expect(recordDrawLists(lists, ['data'])).toContain('glyphs "10px mono" placed #000000ff M 10 20 90 2');
    expect(drawListsToSvg(lists)).toContain('transform="translate(10 20) rotate(90) scale(1 2)"');
  });
});

describe('the window mark', () => {
  it('its parts: the ends within half the hit width, the body between, nothing outside; a floor of 16 px', () => {
    // one geometry with the navigators' (`navigator.window`)
    const half = DefaultRenderConfig.navigator.window.hitPx / 2;
    expect(windowPart(100, 5, 100, 200, 0, 10)).toBe('start');
    expect(windowPart(100 - half, 5, 100, 200, 0, 10)).toBe('start');
    expect(windowPart(200 + half, 5, 100, 200, 0, 10)).toBe('end');
    expect(windowPart(150, 5, 100, 200, 0, 10)).toBe('body');
    expect(windowPart(150, 11, 100, 200, 0, 10)).toBeUndefined();
    expect(windowPart(80, 5, 100, 200, 0, 10)).toBeUndefined();
    // a narrow window grows to the floor about its middle; by touch the ends take more
    expect(windowPart(92, 5, 100, 101, 0, 10)).toBe('start');
    expect(windowPart(100 - 11, 5, 100, 200, 0, 10, true)).toBe('start');
  });

  it('active: a body and two grips; idle: a quiet outline alone', () => {
    const lists = new DrawLists();
    const out = lists.writer('window');
    const draw = (active: boolean, lit: 'start' | undefined) => {
      out.args[0] = 100;
      out.args[1] = 200;
      out.args[2] = 0;
      out.args[3] = 30;
      drawWindowArgs(out, active, lit, LightTokens);
    };
    draw(true, 'start');
    expect(out.pass.rects.n).toBe(3);
    draw(false, undefined);
    expect(out.pass.rects.n).toBe(4);
    expect(recordDrawLists(lists, ['window'])).toMatch(/^size 0 0\npass window\nrect 100 0 100 30 #/);
  });
});

describe('the tracks layout’s navigation', () => {
  const c: NavContext = { total: 400, width: 400, height: 100, cyclic: false, maxScale: 32 };

  it('pan by px, zoom about a unit, fit, reveal, clamp', () => {
    const v = { start: 100, end: 200 };
    expect(TracksNav.pan(c, v, 40, 0)).toMatchObject({ start: 110, end: 210 });
    expect(TracksNav.pan(c, v, -4_000, 0)).toMatchObject({ start: 0, end: 100 });
    expect(TracksNav.zoom(c, v, 2, 150)).toMatchObject({ start: 125, end: 175 });
    expect(TracksNav.zoom(c, v, 2, 100)).toMatchObject({ start: 100, end: 150 });
    // never past the most px a position
    expect(TracksNav.zoom(c, v, 1_000, 150).end - TracksNav.zoom(c, v, 1_000, 150).start).toBeCloseTo(400 / 32);
    expect(TracksNav.fit(c, v)).toMatchObject({ start: 0, end: 400 });
    expect(TracksNav.reveal(c, v, 120, 130)).toBe(v);
    expect(TracksNav.reveal(c, v, 300, 310)).toMatchObject({ start: 210, end: 310 });
    expect(TracksNav.reveal(c, v, 50, 60)).toMatchObject({ start: 50, end: 150 });
    expect(TracksNav.clamp(c, { start: -50, end: 50 })).toMatchObject({ start: 0, end: 100 });
  });

  it('steps: along the order, the ends, none past them; a cyclic axis wraps', () => {
    expect(TracksNav.step(c, 10, 'right')).toBe(11);
    expect(TracksNav.step(c, 0, 'left')).toBe(-1);
    expect(TracksNav.step(c, 10, 'home')).toBe(0);
    // one line: up and down leave it (the cursor crosses to the section above or below)
    expect(TracksNav.step(c, 10, 'up')).toBe(-1);
    expect(TracksNav.step(c, 10, 'down')).toBe(-1);
    expect(TracksNav.step(c, 10, 'end')).toBe(399);
    expect(TracksNav.step({ ...c, cyclic: true }, 399, 'right')).toBe(0);
    expect(TracksNav.clamp({ ...c, cyclic: true }, { start: 390, end: 420 })).toMatchObject({ start: 390, end: 420 });
  });
});
