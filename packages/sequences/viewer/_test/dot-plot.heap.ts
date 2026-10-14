/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot's heap guard (the memory rule; seqspec:viewer#performance): its draw (grid tiles at level 0 over two
// views: the fit, about 3 device px a cell, and zoomed in; the repeat's segments), its draw zoomed out past a device
// pixel a cell (2,400 residues at 0.45 px: level 2 or past it, the pyramid's bins, over two pans), and the pointer
// swept over its cells with an overlay frame each (the partner's hit, both positions lit, the rows and the hovered
// cell) allocate nothing once warm. The positive control: the draw's frames allocating an object each (`leaking`),
// through the same warm-up and measurement, are caught. Its own file: cases in one file share V8's feedback.

import { describe, expect, it, vi } from 'vitest';
import {
  controlCaught,
  describeScaling,
  FRAME_SETTLE_MS,
  FRAME_WARM,
  FRAME_WARM_MS,
  heapScaling,
  isFlat,
  leaking,
  settleFrame,
} from './heap.js';

/** 1,000 steps: they catch about 22 bytes a step, the trimmed 500 about 43 (heap.ts). */
const N = 1_000;
/** The guards' measurement: warmed by count and by time, then `N` frames and `4N` measured. */
const measure = (frame: (i: number) => void) => {
  const t0 = performance.now();
  for (let i = 0; i < FRAME_WARM || performance.now() - t0 < FRAME_WARM_MS; i++) frame(i);
  // Settled after the timed warm-up, as every timed warm-up is (lib/seqstar/AGENTS.md, Tests; LC-20).
  let settled = 0;
  settleFrame(() => frame(settled++), FRAME_SETTLE_MS);
  return heapScaling(N, frame);
};
const guard = (frame: () => void) => {
  const r = measure(frame);
  return isFlat(r) ? '' : describeScaling(r, N);
};

/** A 600-residue sequence with its repeat every 150: the grid's 360,000 cells and the matches' segments. */
async function plot() {
  const { dotPlotDoc, REPEAT } = await import('./dot-plot-fixtures.js');
  const { lit } = await import('./overlay-fixtures.js');
  const { recording } = await import('./composer-fixtures.js');
  const { Dirty } = await import('../src/draw/index.js');
  const seq = Array.from({ length: 20 }, (_, k) => (k % 5 === 4 ? 'GSGSGSGSGSGSGSGSGSGSGSGSGSGSGS' : REPEAT)).join('');
  const t = await lit(dotPlotDoc(seq), { width: 900, height: 600, dpr: 2 });
  t.ctx.behaviors.get('hoverLabel')?.remove();
  const be = recording();
  let i = 0;
  /** Two views, a frame each: the fit, and zoomed in past it. */
  const draw = () => {
    if (i++ & 1) t.vp?.setView(1.4, 0, 'api');
    else t.vp?.setView(6, 600, 'api');
    t.frame(be, Dirty.all);
  };
  const v = t.ctx.layout.axis(t.section);
  const l = t.ctx.layout.get(t.section)?.current;
  if (!v || !l) throw new Error('no layout');
  t.vp?.setView(1.4, 0, 'api');
  t.frame(be, Dirty.all);
  const plotEntry = l.entryCount - 1;
  const top = l.scrolledTop + l.entryTop[plotEntry];
  const xs = Float64Array.of(v.left + 50, v.left + 300, v.left + 520, v.left + 700);
  const ys = Float64Array.of(top + 30, top + 200, top + 90, top + 400);
  let j = 0;
  /** The pointer on a cell and an overlay frame. */
  const hover = () => {
    t.hover(xs[j & 3], ys[(j >> 2) & 3]);
    j++;
    t.overlay(be);
  };
  const { Op } = await import('../src/testing/index.js');
  /** The images the last frame drew. */
  const images = () => {
    let n = 0;
    for (let k = 0; k < be.opCount; k++) if (be.ops[4 * k] === Op.image) n++;
    return n;
  };
  return { t, be, draw, hover, images, dispose: () => t.ctx.dispose() };
}

/** 2,400 residues, its grid one Float32Array, zoomed out to 0.45 px a cell (level 2 or past it), two pans. */
async function zoomedOut() {
  const { dotPlotDoc } = await import('./dot-plot-fixtures.js');
  const { lit } = await import('./overlay-fixtures.js');
  const { recording } = await import('./composer-fixtures.js');
  const { Dirty } = await import('../src/draw/index.js');
  const t = await lit(dotPlotDoc('MKTAYIAKQR'.repeat(240), true), { width: 900, height: 600, dpr: 1 });
  const be = recording();
  let i = 0;
  const draw = () => {
    t.vp?.setView(0.45, i++ & 1 ? 200 : 0, 'api');
    t.frame(be, Dirty.all);
  };
  return { be, draw, dispose: () => t.ctx.dispose() };
}

describe('the dot plot allocates nothing per frame', () => {
  it('its draw over two views: grid tiles at two levels and the segments', async () => {
    vi.resetModules();
    const g = await plot();
    g.draw();
    g.draw();
    // Not vacuous: the tiles and the segments are drawn.
    expect(g.be.pathItems.n).toBeGreaterThan(10);
    expect(g.images()).toBeGreaterThan(0);
    expect(guard(g.draw)).toBe('');
    g.dispose();
  }, 120_000);

  it("its draw zoomed out: past level 1, the pyramid's bins, over two pans", async () => {
    vi.resetModules();
    const g = await zoomedOut();
    g.draw();
    g.draw();
    const { Op } = await import('../src/testing/index.js');
    let images = 0;
    for (let k = 0; k < g.be.opCount; k++) if (g.be.ops[4 * k] === Op.image) images++;
    expect(images).toBeGreaterThan(0);
    expect(guard(g.draw)).toBe('');
    g.dispose();
  }, 180_000);

  it("the pointer over its cells: the partner's hit, both positions lit, the rows and the cell outlined", async () => {
    vi.resetModules();
    const g = await plot();
    g.hover();
    expect(g.t.ctx.hover.live.current.pairIndex).toBeGreaterThanOrEqual(0);
    expect(guard(g.hover)).toBe('');
    g.dispose();
  }, 120_000);

  it("the draw's frames allocating an object each are caught (positive control)", async () => {
    vi.resetModules();
    const g = await plot();
    expect(controlCaught(() => measure(leaking(g.draw)))).toBe(true);
    g.dispose();
  }, 120_000);
});
