/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A grid's y in the viewer (_spec/README.md "ViewModel", Grids): with a window of its own (`follow_y: none`) the view
// keeps it apart from x's, a vertical drag and a zoom move it, a fit shows all of it, and the frame draws it.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { afterEach, describe, expect, it } from 'vitest';
import { DefaultSeqViewerSpec } from '../default-spec.js';
import { createSeqViewer, type SeqViewer } from '../viewer.js';
import { type Frames, installFrames } from './frames.js';

let frames: Frames | undefined;
const viewers: SeqViewer[] = [];

afterEach(() => {
  for (const v of viewers.splice(0)) v.retire();
  frames?.restore();
  frames = undefined;
});

const Plot = '/views/v/plot';

/** A sequence of 60 over x and another of 40 over y, the grid's y its own window. */
function spec() {
  const b = createSeqSpec({ title: 'grid y' });
  const x = b.data.sequence('x', { alphabet: 'rna', letters: 'ACGU'.repeat(15) });
  const y = b.data.sequence('y', { alphabet: 'rna', letters: 'GUCA'.repeat(10) });
  const g = b.data.grid('g', { x, y, values: Float32Array.from({ length: 60 * 40 }, (_, i) => i) });
  b.views.view('v', (v) => {
    v.section('plot', { layout: 'grid', of: x, of_y: y, follow_y: 'none' })
      .track('plot')
      .layer('grid', { type: 'grid', data: g });
  });
  return b.build();
}

describe('a grid’s own y window', () => {
  it('is kept apart from x: a vertical drag pans it, a zoom narrows it with x, a fit shows all of it', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(viewer);
    await viewer.load(spec());
    viewer.resize(800);
    await viewer.settle();
    frames.flush();
    const view = viewer.view;
    expect(view.movesY(Plot)).toBe(true);
    const ys = () => {
      viewer.draw();
      const y = viewer.frame!.layout.sections[0].placement.strips.ys[0]!;
      return [y.u0, y.u1];
    };
    expect(ys()).toEqual([0, 40]);
    // a zoom about y 10: both axes narrower, y about its anchor
    expect(view.zoom(Plot, 2, 30, 10)).toBe(true);
    expect(view.get(Plot)!.y).toEqual({ start: 5, end: 25 });
    expect(ys()).toEqual([5, 25]);
    const x = view.get(Plot)!;
    expect(x.end - x.start).toBe(30);
    // a vertical pan by a quarter of the plot's height: a quarter of y's window
    const height = viewer.placementOf(Plot)!.box.height;
    expect(view.panPx(Plot, 0, height / 4)).toBe(true);
    expect(view.get(Plot)!.y).toEqual({ start: 10, end: 30 });
    expect(view.get(Plot)!.start).toBe(x.start);
    // clamped at y's end
    view.panPx(Plot, 0, 10 * height);
    expect(view.get(Plot)!.y).toEqual({ start: 20, end: 40 });
    view.fit(Plot);
    expect(view.get(Plot)!.y).toBeUndefined();
    expect(ys()).toEqual([0, 40]);
  });

  it('is kept by the history: a snapshot before a y move is recorded, and back returns to it', async () => {
    frames = installFrames();
    const viewer = createSeqViewer({ ...DefaultSeqViewerSpec(), ui: undefined });
    viewers.push(viewer);
    await viewer.load(spec());
    viewer.resize(800);
    await viewer.settle();
    frames.flush();
    viewer.draw();
    const view = viewer.view;
    view.zoom(Plot, 2, 30, 10);
    const before = view.snapshot();
    const x = view.get(Plot)!;
    // a y-only move is a change of the view: recorded (a navigator's drag records its press on release)
    view.panPx(Plot, 0, viewer.placementOf(Plot)!.box.height / 4);
    expect(view.get(Plot)!.y).toEqual({ start: 10, end: 30 });
    expect(view.record(before)).toBe(true);
    expect(view.back()).toBe(true);
    expect(view.get(Plot)!.y).toEqual({ start: 5, end: 25 });
    expect(view.get(Plot)!.start).toBe(x.start);
    expect(view.forward()).toBe(true);
    expect(view.get(Plot)!.y).toEqual({ start: 10, end: 30 });
    // a jump that moves only y is a move, and back undoes it
    expect(view.jump(() => view.panPx(Plot, 0, -viewer.placementOf(Plot)!.box.height / 4))).toBe('moved');
    expect(view.get(Plot)!.y).toEqual({ start: 5, end: 25 });
    expect(view.back()).toBe(true);
    expect(view.get(Plot)!.y).toEqual({ start: 10, end: 30 });
  });
});
