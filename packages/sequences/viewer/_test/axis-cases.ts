/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The cases of the axis and viewport guards (VC3; axis.heap.ts) and their quick twin (axis.quick.ts), a step a move:
// the per-move axis lookups (`indexOf`, `shownOf`, `breakAt`, the view's `x`) over a forward ranged unit, a filtered, a
// reversed and a cyclic piece; both run converters; viewport pans and anchored zooms with a live subscriber;
// programmatic (`api`) moves, each committed at once, with no `committed` subscriber, and with one (its records are its
// own, and show). The heap guard and the quick twin share them.

import { AxisGeometry, SectionAxisView, type ShownIndex } from '../src/axis.js';
import { axisOf, open, rangedDoc, resolved, setPiece, vp } from './axis-fixtures.js';

/** The heap guard's steps in the shorter run; the longer runs four times as many. */
export const N = 50_000;

let sink = 0;
/** What the steps summed, for the guards' non-vacuity. */
export const sunk = () => sink;
const out = new Int32Array(16);

/** The four axes: forward, filtered, reversed and cyclic. */
async function axes() {
  const [v, f, w, c] = await Promise.all([
    resolved(rangedDoc(), 'v'),
    resolved(rangedDoc(), 'f'),
    resolved(rangedDoc(), 'w'),
    resolved(rangedDoc(), 'w'),
  ]);
  setPiece(w.view, 0, 0, 20, -1, 20);
  setPiece(c.view, 0, 15, 5, 1, 20);
  return [axisOf(v), axisOf(f), axisOf(w), axisOf(c)];
}

/** `indexOf`, `shownOf`, `breakAt` and `x` over the four axes. */
export async function lookups(): Promise<(i: number) => void> {
  const views = (await axes()).map((ax) => {
    const v = new SectionAxisView();
    const l = { geometry: AxisGeometry.linear, dataLeft: 10, dataRight: 300, scrolledTop: 0, height: 100 };
    const s = { cell: 9, scrollX: 3, scrollY: 0, dpr: 2, lod: 0, contentHeight: 0 };
    v.sync(ax, l, { ...s, visLo: new Int32Array(1), visHi: new Int32Array(1) } as never);
    return v;
  });
  return (i) => {
    for (let a = 0; a < views.length; a++) {
      const v = views[a];
      const k = (i % 12) as ShownIndex;
      sink = (sink + v.indexOf(0, k) + v.shownOf(0, i % 21) + (v.breakAt(0, k) ? 1 : 0) + v.x(0, k)) | 0;
    }
  };
}

/** `runToShown` and `shownToRuns` over the four axes. */
export async function runs(): Promise<(i: number) => void> {
  const all = await axes();
  return (i) => {
    for (let a = 0; a < all.length; a++) {
      const ax = all[a];
      sink = (sink + ax.runToShown(0, i % 20, (i % 20) + 4, out) + ax.shownToRuns(0, i % 8, (i % 8) + 3, out)) | 0;
    }
  };
}

/** A viewport on the managers' document, zoomed to cell 28 (the most), centered: every move changes the view. */
async function viewport() {
  const { ctx } = await open();
  ctx.views.setCurrent('g');
  await ctx.graph.update();
  const v = vp(ctx, 'g.one');
  v.zoomTo(0, 8, 12);
  return v;
}

/** Pans and anchored zooms at pointer rate, with a live subscriber; `moves` counts the versions bumped. */
export async function panning() {
  const v = await viewport();
  v.live.version.subscribe(() => {
    sink = (sink + v.live.current.scrollX + v.live.current.visLo[0]) & 0xffff;
  });
  const before = v.live.version.value;
  const step = (k: number) => {
    // Constants per branch: a double chosen by `?:` is boxed by the caller before the call.
    if (!(k & 64)) v.panBy(k & 1 ? 3 : -3, 0, 'pointer');
    else if (k & 1) v.zoomAt(1.25, 40, 'pointer');
    else v.zoomAt(0.8, 40, 'pointer');
  };
  return { step, v, moves: () => v.live.version.value - before };
}

/** Pans and views set with origin api: committed each, built only when read (nobody reads `committed`). */
export async function api() {
  const v = await viewport();
  const before = v.committed.value;
  const step = (k: number) => {
    if (k & 1) v.panBy(k & 2 ? 3 : -3, 0, 'api');
    else v.setView(k & 2 ? 20 : 24, 40, 'api');
  };
  return { step, v, before };
}

/** An api pan with a `committed` subscriber: each move's record is delivered, and shows. */
export async function subscribed() {
  const v = await viewport();
  let n = 0;
  v.committed.subscribe(() => n++);
  return { step: (k: number) => v.panBy(k & 1 ? 3 : -3, 0, 'api'), delivered: () => n };
}
