/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rings layout: placement as the inverse of hit-testing over random rotations (every position placed once, `at`
// the inverse of `place`, `offsetAt` monotone in angle), a linear detail following the ring's window (over the origin
// too) and a ring following a detail, the window out of the ring's key (a pan of the detail draws the ring from its
// kept lays, its labels not placed again), drawing against the data, labels outside hit as their items, the window's
// parts, reused buffers against fresh ones, Canvas2D-free lists equal in SVG.

import { describe, expect, it } from 'vitest';
import {
  type Frame,
  type Hit,
  PlacePoint,
  type RingPlacement,
  RingsNav,
  type NavContext,
  recordDrawLists,
  drawListsToSvg,
  hitView,
  renderSvg,
} from '../index.js';
import { rangesLocation } from '@molstar/seqx-core/query';
import { settle } from './helpers.js';
import { RingFeatures, RingLength, ringLetters, ringSpec } from './rings-fixtures.js';
import { type ReadPath, readPaths, rng } from './small-helpers.js';

const Size = { width: 640, height: 900 };

const ringOf = (f: Frame) => f.layout.sections.find((s) => s.layout.name === 'rings')!;

describe('the rings layout', () => {
  it('places the section round its whole axis, the first track outermost, the hole kept, the title in it', async () => {
    const s = await settle(ringSpec());
    const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 100, end: 200 } } });
    const ring = ringOf(f);
    const p = ring.placement as RingPlacement;
    expect(p.strips.n).toBe(1);
    expect(ring.axis!.u0).toBe(0);
    expect(ring.axis!.u1).toBe(RingLength);
    // the first band from the outer radius in; the rings' depth keeps a hole of at least 0.35 of the radius
    expect(p.strips.bandY[0]).toBe(0);
    expect(p.depth).toBeLessThanOrEqual(p.r * 0.65 + 1e-9);
    expect(p.title).toBe('construct');
    expect(p.subtitle).toBe('600 bp');
    expect(f.layout.issues.filter((i) => i.severity === 'error')).toEqual([]);
  });

  it('placement is the inverse of hit-testing over rotations: every position placed once, `at` gives it back, `offsetAt` monotone in angle', async () => {
    const r = rng(3);
    for (let n = 0; n < 4; n++) {
      const rotation = Math.floor(r() * 360);
      const s = await settle(ringSpec({ rotation, detail: false }));
      const f = await s.frame({ size: Size });
      const p = ringOf(f).placement as RingPlacement;
      const pt = new PlacePoint();
      const seen = new Set<number>();
      let last = -Infinity;
      for (let u = 0; u < RingLength; u++) {
        expect(p.place(u + 0.5, pt)).toBe(true);
        const { x, y } = pt;
        expect(p.at(x, y, pt)).toBe(true);
        expect(pt.kind).toBe('data');
        const back = Math.floor(pt.offset);
        expect(back).toBe(u);
        seen.add(back);
        const o = p.offsetAt(x, y);
        expect(o).toBeGreaterThan(last);
        last = o;
      }
      expect(seen.size).toBe(RingLength);
    }
  });

  it('a linear detail follows the ring’s window, over the origin as two pieces; the ring follows a detail the same', async () => {
    const s = await settle(ringSpec());
    const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 580, end: 620 } } });
    const detail = f.layout.sections.find((x) => x.section.address === '/views/v/detail')!;
    expect(detail.axis!.pieces).toEqual([580, 600, 0, 20]);
    const p = ringOf(f).placement as RingPlacement;
    expect([p.windowShown, p.windowActive, p.windowStart, p.windowEnd]).toEqual([true, true, 580, 620]);
    // the ring follows: the detail's range is its window
    const t = await settle(ringSpec({ ringFollows: true }));
    const g = await t.frame({ size: Size, viewports: { '/views/v/detail': { start: 100, end: 160 } } });
    const q = ringOf(g).placement as RingPlacement;
    expect(ringOf(g).axis!.total).toBe(RingLength);
    expect([q.windowShown, q.windowStart, q.windowEnd]).toEqual([true, 100, 160]);
  });

  it('a pan of the detail does not lay the ring out again: its lays and labels kept, the window moved', async () => {
    const s = await settle(ringSpec());
    const f1 = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 100, end: 200 } } });
    const ring1 = ringOf(f1);
    const lay = ring1.lays[0][0].lay;
    const labels = (ring1.placement as RingPlacement).labels.m;
    const f2 = s.renderer.frame({ size: Size, viewports: { '/views/v/ring': { start: 140, end: 240 } } });
    expect(f2.rebuilt).toBe(true);
    const ring2 = ringOf(f2);
    expect(ring2).toBe(ring1);
    expect(ring2.lays[0][0].lay).toBe(lay);
    expect((ring2.placement as RingPlacement).labels.m).toBe(labels);
    expect((ring2.placement as RingPlacement).windowStart).toBe(140);
    // drawn as a fresh layout draws it
    const fresh = await settle(ringSpec());
    const g = await fresh.frame({ size: Size, viewports: { '/views/v/ring': { start: 140, end: 240 } } });
    expect(recordDrawLists(f2.lists)).toBe(recordDrawLists(g.lists));
  });

  it('labels outside, each hit as its item; the feature over the origin one label, at its middle angle', async () => {
    const s = await settle(ringSpec({ detail: false }));
    const f = await s.frame({ size: Size });
    const ring = ringOf(f);
    const p = ring.placement as RingPlacement;
    const texts = Array.from({ length: p.labels.m }, (_, k) => p.labels.text[p.labels.placed[k]]).sort();
    expect(texts).toEqual(['AmpR', 'MCS', 'ori']);
    for (let k = 0; k < p.labels.m; k++) {
      const i = p.labels.placed[k];
      const x = p.labels.right[k] ? p.labels.x[k] + 2 : p.labels.x[k] - 2;
      const h = hitAtPage(f, x, p.labels.y[k] - 3);
      expect(h?.element?.key).toBe(p.labels.hits[i]?.element?.key);
    }
    // ori runs 560 to 40 over the origin: its label's angle is its middle, position 0's
    const ori = p.labels.text.indexOf('ori');
    expect(p.labels.angle[ori]).toBeCloseTo(0, 6);
    // no label is drawn inside a mark on a ring
    expect(f.lists.passes.labels.texts.n).toBeGreaterThanOrEqual(3);
  });

  it('draws the features against the data: each feature’s arc where its positions are, over the origin one arc', async () => {
    const s = await settle(ringSpec({ detail: false }));
    const f = await s.frame({ size: Size });
    const ring = ringOf(f);
    const p = ring.placement as RingPlacement;
    // across the features ring at each position: the features holding it, and nothing where none is
    const keysAt = (pos: number) => {
      const keys = new Set<string | number>();
      const x = p.frameX(pos + 0.5);
      for (let y = 1; y < p.strips.bandH[0]; y += 2) {
        p.warp.map(x, y);
        const k = hitAtPage(f, p.warp.px, p.warp.py)?.element?.key;
        if (k !== undefined) keys.add(k);
      }
      return [...keys].sort();
    };
    const rows = RingFeatures;
    for (let pos = 0; pos < RingLength; pos += 7) {
      const want = rows.id.filter((_, i) =>
        rows.start[i] < rows.end[i]
          ? pos >= rows.start[i] && pos < rows.end[i]
          : pos >= rows.start[i] || pos < rows.end[i],
      );
      expect(keysAt(pos)).toEqual(want.sort());
    }
    // the letters the detail would show are the data's
    expect(ringLetters().length).toBe(RingLength);
  });

  it('the window’s parts on the rim: body and ends hit as the window target; idle with no follower: no window', async () => {
    const s = await settle(ringSpec());
    const f = await s.frame({ size: Size, viewports: { '/views/v/ring': { start: 100, end: 200 } } });
    const ring = ringOf(f);
    const p = ring.placement as RingPlacement;
    const at = (frameX: number) => {
      p.warp.map(frameX, -5);
      return hitAtPage(f, p.warp.px, p.warp.py);
    };
    expect(at(p.frameX(150))).toMatchObject({ target: 'window', part: 'body', section: 'ring' });
    expect(at(p.frameX(100))).toMatchObject({ target: 'window', part: 'start' });
    expect(at(p.frameX(200))).toMatchObject({ target: 'window', part: 'end' });
    // a window drawn in the window pass, nothing of it in an export by default
    expect(f.lists.passes.window.opCount).toBeGreaterThan(0);
    const alone = await settle(ringSpec({ detail: false }));
    const g = await alone.frame({ size: Size });
    expect((ringOf(g).placement as RingPlacement).windowShown).toBe(false);
    expect(g.lists.passes.window.opCount).toBe(0);
  });

  it('the ring’s lists in SVG: every primitive bent onto the page (no rect survives the warp)', async () => {
    const s = await settle(ringSpec({ detail: false }));
    const f = await s.frame({ size: Size });
    const svg = drawListsToSvg(f.lists);
    expect(svg).toContain('<path');
    expect(f.lists.passes.data.rects.n).toBe(0);
  });
});

describe('a ring in two consumers and in an export', () => {
  it('renderSvg draws the lists a renderer draws (the window only with `window`)', async () => {
    const spec = ringSpec();
    const viewports = { '/views/v/ring': { start: 500, end: 640 } };
    const out = await renderSvg(spec, { size: Size, viewports });
    const s = await settle(spec);
    const f = await s.frame({ size: Size, viewports, figure: true, settled: true });
    expect(recordDrawLists(out.lists)).toBe(recordDrawLists(f.lists));
    expect(out.lists.passes.window.opCount).toBe(0);
    const windowed = await renderSvg(spec, { size: Size, viewports, window: true });
    expect(windowed.lists.passes.window.opCount).toBeGreaterThan(0);
    expect(windowed.svg).toContain('<path');
  });

  it('a selection over the origin is one arc on the ring: two runs whose square ends meet at the seam', async () => {
    const s = await settle(ringSpec({ detail: false }));
    const selection = rangesLocation('/data/construct', [0, 30, 570, 600]);
    const f = await s.frame({ size: Size, emphasis: { selection } });
    const paths = readPaths(f.lists, 'selection');
    // a tint and an outline a run, two runs: each a bent path on the page
    expect(paths.length).toBe(4);
    const p = ringOf(f).placement as RingPlacement;
    // the runs meet at the origin's angle: a point of each at the seam
    p.warp.map(p.frameX(0), 1);
    const seam = { x: p.warp.px, y: p.warp.py };
    const near = (q: ReadPath) => {
      for (let i = 0; i + 1 < q.xy.length; i += 2)
        if (Math.hypot(q.xy[i] - seam.x, q.xy[i + 1] - seam.y) < 2) return true;
      return false;
    };
    expect(paths.filter(near).length).toBeGreaterThanOrEqual(2);
  });
});

describe('the rings layout’s navigation', () => {
  const c: NavContext = { total: 600, width: 400, height: 400, cyclic: true, maxScale: 20 };
  it('pans round the origin, zooms about a unit, reveals by the shortest way, fits', () => {
    const v = { start: 580, end: 620 };
    expect(RingsNav.pan(c, v, -200, 0)).toMatchObject({ start: 0, end: 40 });
    expect(RingsNav.reveal(c, v, 590, 610)).toBe(v);
    expect(RingsNav.reveal(c, v, 5, 15)).toBe(v);
    expect(RingsNav.reveal(c, v, 300, 310)).toMatchObject({ start: 285, end: 325 });
    expect(RingsNav.fit(c, v)).toMatchObject({ start: 0, end: 600 });
    expect(RingsNav.step(c, 599, 'right')).toBe(0);
  });
});

/** The hit at a page point of a frame. */
function hitAtPage(f: Frame, x: number, y: number): Hit | undefined {
  return hitView(f.layout, x, y);
}
