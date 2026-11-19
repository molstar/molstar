/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The navigators' stories, headless on a fake clock, composed as the stories compose them (a viewer of the view and its
// SVG with the navigators included): chr17's ideogram over the gene's region (a band's name under the pointer, the
// region's window), p53's tracks with a strip (it replaces the rail; a drag on its window pans), the protein's strip
// with its gene's window (a drag on the gene's moves the protein), the Globins' overview and miniature (a drag moves the
// rows), the plasmid's strip over the detail (a window over the origin in two pieces); and the feature's rail.

import { GlobinsViews, P53Views, PlasmidViews } from '@molstar/seqx-examples';
import { regionSeqText } from '@molstar/seqx-examples/p53/region';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from '../../../examples/_test/draw-checks.js';
import { openStrip } from '../compose.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

/** A navigator of the viewer drawn, its picture's wants answered on the fake clock. */
async function navigator(v: SeqViewer, id: string) {
  const nav = v.navigatorModel(id)!;
  nav.draw();
  await vi.advanceTimersByTimeAsync(600);
  v.draw();
  return { nav, f: nav.draw()! };
}

/** The SVG as the story draws it: the navigators included. */
async function svg(v: SeqViewer) {
  const p = v.renderSvg({ navigators: 'include' });
  await vi.advanceTimersByTimeAsync(600);
  return p;
}

describe('the navigators’ examples, headless', () => {
  it('p53’s gene region on chr17: the region’s window on the ideogram, a band named under the pointer', async () => {
    viewer = await openExample(P53Views.find((v) => v.id === 'gene-region-ideogram')!, { gene: regionSeqText() });
    const { nav, f } = await navigator(viewer, 'chr17');
    const w = f.windows[0];
    expect(w.shown).toBe(true);
    // the region is a sliver of the chromosome: its window on the zoom strip under the ideogram
    expect(nav.height).toBe(24 + 16);
    expect(w.ly).toBeGreaterThanOrEqual(24);
    nav.move(nav.width * 0.6, 12);
    expect(nav.readout.value).toMatch(/^chr17:[\d,]+ · q\d/);
    // over the window: what it shows
    nav.move((w.dx0 + w.dx1) / 2, (w.ly + w.hy) / 2);
    expect(nav.readout.value).toBe('chr17:7,666,403–7,672,402 of chr17');
    const r = await svg(viewer);
    expect(r.settled).toBe(true);
    expect(r.svg).toContain('data-navigator="chr17"');
    expect(r.svg).toContain('<title>Window: chr17:');
  });

  it('p53’s tracks with a strip: the strip replaces the rail; a drag on its window pans', async () => {
    viewer = await openExample(P53Views.find((v) => v.id === 'tracks-strip')!, { gene: '' });
    const v = viewer;
    expect(v.frame!.layout.sections[0].rail).toBeUndefined();
    v.view.setRange('/views/main/protein', 100, 150);
    v.draw();
    const { nav, f } = await navigator(v, 'strip');
    const w = f.windows[0];
    const mid = (w.dx0 + w.dx1) / 2;
    expect(nav.down(mid, 20)).toBe(true);
    nav.move(mid + (40 / 393) * nav.width, 20);
    nav.up();
    expect(v.view.get('/views/main/protein')!.start).toBeCloseTo(140, 0);
    const r = await svg(v);
    expect(r.svg).toContain('data-navigator="strip"');
  });

  it('the Strip story opens zoomed in: its window is a part of the strip, so it can be dragged', async () => {
    viewer = await openExample(P53Views.find((v) => v.id === 'tracks-strip')!, { gene: '' });
    const v = viewer;
    const whole = await navigator(v, 'strip');
    const full = whole.f.windows[0];
    expect(full.dx1 - full.dx0).toBeGreaterThan(whole.nav.width * 0.95);
    openStrip(v);
    v.draw();
    const { nav, f } = await navigator(v, 'strip');
    const w = f.windows[0];
    expect(w.dx1 - w.dx0).toBeLessThan(nav.width * 0.5);
    const before = v.view.get('/views/main/protein')!.start;
    const mid = (w.dx0 + w.dx1) / 2;
    expect(nav.down(mid, 20)).toBe(true);
    nav.move(mid + 30, 20);
    nav.up();
    expect(v.view.get('/views/main/protein')!.start).toBeGreaterThan(before);
  });

  it('the protein’s strip with its gene’s window: a drag on the gene’s window moves the protein, which leads', async () => {
    viewer = await openExample(P53Views.find((v) => v.id === 'protein-gene-strip')!, { gene: '' });
    const v = viewer;
    v.view.setRange('/views/main/protein', 100, 150);
    v.draw();
    const { nav, f } = await navigator(v, 'strip');
    const [p, g] = f.windows;
    expect(p.shown && g.shown).toBe(true);
    const mid = (g.dx0 + g.dx1) / 2;
    expect(nav.hit(mid, 20)?.section).toBeDefined();
    nav.down(mid, 20);
    nav.move(mid + (30 / 393) * nav.width, 20);
    nav.up();
    expect(v.view.get('/views/main/protein')!.start).toBeGreaterThan(110);
  });

  it('the Globins’ 10,000 rows: a drag on the overview’s window moves the rows; the miniature covers both axes', {
    tags: ['medium'],
  }, async () => {
    viewer = await openExample(GlobinsViews.find((v) => v.id === 'scale-navigators')!, { fixtures: true });
    const v = viewer;
    const { nav, f } = await navigator(v, 'overview');
    const w = f.windows[0];
    expect(w.rows).toBe(10_000);
    const x = (w.dx0 + w.dx1) / 2;
    const y = (w.dy0 + w.dy1) / 2;
    expect(nav.down(x, y)).toBe(true);
    nav.move(x, y + 100);
    nav.up();
    expect(v.view.band('/views/main/rows/rows')!.bandTop).toBeGreaterThan(1000);
    // the rows scrolled: the overview's window is the rows in view, at their place in the column
    v.view.setBand('/views/main/rows/rows', { bandTop: 4000 });
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    const scrolled = nav.draw()!.windows[0];
    const shown = scrolled.r1 - scrolled.r0;
    expect(scrolled.r0).toBeCloseTo(4000, 0);
    expect(shown).toBeGreaterThan(0);
    const column = scrolled.hy - scrolled.ly;
    expect((scrolled.y0 - scrolled.ly) / column).toBeCloseTo(0.4, 2);
    expect((scrolled.y1 - scrolled.y0) / column).toBeCloseTo(shown / 10_000, 4);
    // drawn: the band across the column at its rows (floored), its grips the rows axis's, no positions hairline
    expect(scrolled.axes).toBe(2);
    expect(scrolled.dy0).toBeLessThanOrEqual(scrolled.y0 + 1e-6);
    expect(scrolled.dy1).toBeGreaterThanOrEqual(scrolled.y1 - 1e-6);
    nav.move(nav.width / 2, (scrolled.dy0 + scrolled.dy1) / 2);
    // over the window the readout names the rows shown
    expect(nav.readout.value).toMatch(/^rows 4,001–[\d,]+ of 10,000$/);
    const marks = nav.draw()!.lists.passes.hover.rects;
    // the hover is the row's hairline across the column alone, never a positions hairline down it
    expect(marks.n).toBe(1);
    expect(marks.h[0]).toBeLessThan(2);
    const mini = await navigator(v, 'miniature');
    expect(mini.f.windows[0].axes).toBe(3);
    const r = await svg(v);
    expect(r.svg).toContain('data-navigator="overview"');
    expect(r.svg).toContain('data-navigator="miniature"');
  });

  it('the plasmid’s strip over the detail: a window over the origin in two pieces', async () => {
    viewer = await openExample(PlasmidViews.find((v) => v.id === 'detail-strip')!, { gene: '' });
    const v = viewer;
    v.view.setRange('/views/main/ring', 2600, 2760);
    await vi.advanceTimersByTimeAsync(50);
    v.draw();
    const { f } = await navigator(v, 'strip');
    expect(f.windows[0].wrap).toBe(true);
    const r = await svg(v);
    expect(r.svg).toContain('data-navigator="strip"');
  });
});

describe('the zoom strip of chr17’s ideogram over the gene’s region, headless', () => {
  const Region = '/views/main/region';

  /** The gene region's viewer at a range of the region, its ideogram drawn; the region's length and the strip's px. */
  async function zoomed(start: number, end: number) {
    viewer = await openExample(P53Views.find((v) => v.id === 'gene-region-ideogram')!, { gene: regionSeqText() });
    const v = viewer;
    v.view.setRange(Region, start, end);
    v.draw();
    const { nav, f } = await navigator(v, 'chr17');
    const total = v.frame!.layout.sections.find((s) => s.section.address === Region)!.axis!.total;
    return { v, nav, f, total, perPx: total / nav.width };
  }

  it('a drag on the zoom strip’s body pans the view by the zoom strip’s scale', async () => {
    const { v, nav, f, perPx } = await zoomed(1000, 2000);
    const w = f.windows[0];
    const x = (w.dx0 + w.dx1) / 2;
    const y = (w.ly + w.hy) / 2;
    expect(nav.hit(x, y)?.part).toBe('body');
    expect(nav.down(x, y)).toBe(true);
    nav.move(x + 40, y);
    nav.up();
    const r = v.view.get(Region)!;
    expect(r.start).toBeCloseTo(1000 + 40 * perPx, -1);
    expect(r.end - r.start).toBeCloseTo(1000, 0);
  });

  it('a grip zooms: the end dragged, the start held', async () => {
    const { v, nav, f, perPx } = await zoomed(1000, 2000);
    const w = f.windows[0];
    const y = (w.ly + w.hy) / 2;
    const x = w.dx1 + 2;
    expect(nav.hit(x, y)?.part).toBe('end');
    expect(nav.down(x, y)).toBe(true);
    nav.move(x + 40, y);
    nav.up();
    const r = v.view.get(Region)!;
    expect(r.start).toBeCloseTo(1000, 0);
    // the end where the pointer is, at the zoom strip's scale (a px of it is `perPx` positions, not chr17's 100 kb)
    expect(Math.abs(r.end - (x + 40) * perPx)).toBeLessThan(perPx);
  });

  it('a view zoom or pan moves the zoom strip’s window', async () => {
    const { v, nav, f, perPx } = await zoomed(1000, 2000);
    const before = f.windows[0];
    expect(before.x0).toBeCloseTo(1000 / perPx, 0);
    v.view.setRange(Region, 3000, 3500);
    v.draw();
    const after = nav.draw()!.windows[0];
    expect(after.x0).toBeCloseTo(3000 / perPx, 0);
    expect(after.x1 - after.x0).toBeCloseTo(500 / perPx, 0);
  });

  it('a press on the mark is no drag, its hover names the region; the context is no target', async () => {
    const { v, nav, f } = await zoomed(1000, 2000);
    const win = f.lists.passes.window.rects;
    // the mark: the outline over the ideogram's height
    const i = Array.from({ length: win.n }, (_, k) => k).find((k) => win.y[k] === 0 && win.h[k] === 24)!;
    const mx = win.x[i] + win.w[i] / 2;
    nav.move(mx, 12);
    expect(nav.readout.value).toMatch(/^chr17:7,666,403–[\d,]+$/);
    expect(nav.cursor).toBe('default');
    const before = v.view.get(Region)!;
    expect(nav.down(mx, 12)).toBe(false);
    nav.move(mx + 40, 12);
    nav.up();
    expect(nav.click(mx + 40, 12)).toBe(false);
    expect(v.view.get(Region)).toEqual(before);
    // elsewhere on the ideogram: a band's name, no centring
    nav.move(nav.width * 0.6, 12);
    expect(nav.readout.value).toMatch(/^chr17:[\d,]+ · q\d/);
    expect(nav.down(nav.width * 0.6, 12)).toBe(false);
    expect(nav.doubleClick(nav.width * 0.6, 12)).toBe(false);
    expect(v.view.get(Region)).toEqual(before);
  });
});

describe('the navigators feature, headless', () => {
  it('the rail: every linear section has one, live; an export leaves it out', async () => {
    viewer = await openExample(P53Views.find((v) => v.id === 'tracks')!, { gene: '' });
    const v = viewer;
    v.view.setRange('/views/main/protein', 100, 150);
    v.draw();
    const rail = v.frame!.layout.sections[0].rail!;
    expect(rail.active).toBe(true);
    const x = rail.x + (125 / rail.total) * rail.width;
    expect(v.hit(x, rail.bottom - 2)).toMatchObject({ target: 'window', part: 'body', navigator: 'rail' });
    const p = v.renderSvg();
    await vi.advanceTimersByTimeAsync(400);
    expect((await p).lists.passes.window.rects.n).toBe(0);
  });
});
