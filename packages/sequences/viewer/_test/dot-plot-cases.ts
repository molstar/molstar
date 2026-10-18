/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The dot plot's cases, for its heap guard and its quick twin (dot-plot.heap.ts, dot-plot.quick.ts): `plot`, a
// 600-residue sequence with its repeat every 150 (the grid's 360,000 cells and the matches' segments), its draw over
// two views (the fit and zoomed in) and the pointer over its cells with an overlay frame each; `zoomedOut`, 2,400
// residues drawn at 0.45 px a cell (level 2 or past it) over two pans. `images` counts the image ops of the last frame.
// Each call loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';
import type { AllocOptions } from '../../core/_test/model/alloc.js';

/** The draws' quick counts: a whole frame of the plot is about 0.5 ms (M4), so windows of 250 steps. */
export const DRAW: AllocOptions = { warmSteps: 1_000, steps: 250, windows: 3 };

/** The image ops a recording backend's last frame drew. */
async function imageCount() {
  const { Op } = await import('../src/testing/index.js');
  return (be: { opCount: number; ops: ArrayLike<number> }) => {
    let n = 0;
    for (let k = 0; k < be.opCount; k++) if (be.ops[4 * k] === Op.image) n++;
    return n;
  };
}

/** A 600-residue sequence with its repeat every 150: the grid's 360,000 cells and the matches' segments. */
export async function plot() {
  vi.resetModules();
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
  const images = await imageCount();
  return { t, be, draw, hover, images: () => images(be), dispose: () => t.ctx.dispose() };
}

/** 2,400 residues, its grid one Float32Array, zoomed out to 0.45 px a cell (level 2 or past it), two pans. */
export async function zoomedOut() {
  vi.resetModules();
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
  const images = await imageCount();
  return { be, draw, images: () => images(be), dispose: () => t.ctx.dispose() };
}
