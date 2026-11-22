/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tRNA and the Complex's PAE and contacts stories, headless on a fake clock, composed as the stories compose them
// (a viewer of the view and its SVG export), and the interaction matrix's cells they fill, as bindings: a pair (click,
// Cmd/Ctrl+click, double-click), a grid's cell (click, Shift+click, Cmd/Ctrl+click), a rectangle (drag, Cmd/Ctrl+drag
// adds, Alt+drag subtracts), the crosshair hover, a grid's zoom moving its leaders, the arcs the canvas and the export
// draw (equal lists).

import { type GridLay, type GridPlacement, type PairsLay, recordDrawLists } from '@molstar/seqx-core/render';
import { ComplexViews, TrnaViews } from '@molstar/seqx-examples';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { openExample } from './draw-checks.js';

let viewer: SeqViewer | undefined;

afterEach(() => {
  viewer?.retire();
  viewer = undefined;
  vi.useRealTimers();
});

const trna = (id: string) => TrnaViews.find((v) => v.id === id)!;
const complex = (id: string) => ComplexViews.find((v) => v.id === id)!;

/** The point of a cubic at t. */
function curveAt(c: Float64Array, o: number, t: number): [number, number] {
  const u = 1 - t;
  const f = (k: number) =>
    u * u * u * c[o + k] + 3 * u * u * t * c[o + 2 + k] + 3 * u * t * t * c[o + 4 + k] + t * t * t * c[o + 6 + k];
  return [f(0), f(1)];
}

/** The arcs' lay of a section (its first `section` layer). */
const arcsOf = (v: SeqViewer, section = 0) => v.frame!.layout.sections[section].sectionLays[0].lay as PairsLay;

/** A grid section's layer lay (its one track's `k`th layer). */
const gridOf = (v: SeqViewer, section: number, k = 0) => v.frame!.layout.sections[section].lays[0][k].lay as GridLay;

/** The page centre of a grid's cell. */
function cellCentre(g: GridLay, i: number, j: number): [number, number] {
  g.axis!.cell(i);
  const x = (g.axis!.sx0 + g.axis!.sx1) / 2;
  g.axisY!.cell(j);
  return [x, (g.axisY!.sx0 + g.axisY!.sx1) / 2];
}

/**
 * The grid layout's size rule at aspect 1: a plot's box is fixed, the data box's width (as tall as wide) capped by the
 * room (the larger of a grid track's default height and nine tenths of the view below the plot's top), whatever the
 * windows; cells square only where x and y show spans in that ratio.
 */
function expectAspect(v: SeqViewer, section: number) {
  const st = v.frame!.layout.sections[section].placement.strips;
  const x = st.axes[0];
  const room = Math.max(240, (v.height - st.top[0]) * 0.9);
  expect(st.height[0]).toBeCloseTo(Math.max(40, Math.min(x.width, room)), 6);
}

/** The paths of the data pass, as recorded (the arcs). */
const paths = (lists: Parameters<typeof recordDrawLists>[0]) =>
  recordDrawLists(lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('path'));

describe('tRNA: arcs, headless', () => {
  it('every pair has its arc; the export draws the canvas’s arcs; a base hovered lights its partner’s arc', async () => {
    viewer = await openExample(trna('arcs'), { gene: '' });
    const v = viewer;
    // the rails off: an export has none, nor their row under the ruler
    v.setRailsShown(false);
    v.draw();
    const lay = arcsOf(v);
    expect(lay.n).toBe(27);
    const svg = v.renderSvg({ idPrefix: 't-' });
    await vi.advanceTimersByTimeAsync(400);
    const r = await svg;
    expect(r.issues).toEqual([]);
    expect(paths(r.lists)).toEqual(paths(v.frame!.lists));
    // hover the letters at position 9 (1-based 10): the arc 10–25 lights
    const letters = v.frame!.layout.sections[0].tracks.find((t) => t.track.address.endsWith('/letters'))!.box;
    const a = v.frame!.layout.sections[0].axis!;
    a.cell(9);
    v.input.move((a.sx0 + a.sx1) / 2, letters.y + letters.height / 2);
    v.draw();
    expect(v.hover.current.position).toBe(9);
    expect(v.frame!.lists.passes.hover.paths.n).toBe(1);
  });

  it('a pair: click selects both ends, Cmd/Ctrl+click toggles it out, double-click focuses both ends', async () => {
    viewer = await openExample(trna('arcs'), { gene: '' });
    const v = viewer;
    v.input.mac = false;
    const lay = arcsOf(v);
    // the acceptor stem's outermost pair, 1–72
    const i = Array.from({ length: lay.n }, (_, k) => k).find((k) => lay.prep!.data!.keys[lay.pair[k]] === '1-72')!;
    const [x, y] = curveAt(lay.curve, 8 * i, 0.5);
    expect(v.input.click(x, y)).toBe(true);
    const sel = v.selection.value.value;
    expect(sel.elements.keys.map((k) => k.key)).toEqual(['1-72']);
    expect([...sel.shown.parts[0].ranges]).toEqual([0, 1, 71, 72]);
    v.input.click(x, y, { ctrlKey: true });
    expect(v.selection.value.value.elements.keys).toEqual([]);
    // the anticodon stem's pair 30–40, then focus: the window its ends' hull
    const j = Array.from({ length: lay.n }, (_, k) => k).find((k) => lay.prep!.data!.keys[lay.pair[k]] === '30-40')!;
    const [x2, y2] = curveAt(lay.curve, 8 * j, 0.5);
    v.input.doubleClick(x2, y2);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    const w = v.view.get('/views/main/trna')!;
    expect(w.start).toBeLessThanOrEqual(29);
    expect(w.end).toBeGreaterThanOrEqual(40);
    expect(w.end - w.start).toBeLessThan(76);
  });

  it('split and chords: the crossing pairs below, the ring’s chords through its hole', async () => {
    viewer = await openExample(trna('split'), { gene: '' });
    expect(arcsOf(viewer).prep!.other!.reduce((s, x) => s + x, 0)).toBeGreaterThan(0);
    viewer.retire();
    viewer = await openExample(trna('ring'), { gene: '' });
    const lay = arcsOf(viewer);
    expect(lay.n).toBe(27);
    const [x, y] = curveAt(lay.curve, 0, 0.5);
    expect(viewer.hit(x, y)?.element?.owner).toBe('/data/pairs');
  });
});

describe('tRNA: the dot plot, headless', () => {
  it('a cell hovered: the crosshair, both bases lit on the letters and the arcs ending at either (projected)', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    const v = viewer;
    const g = gridOf(v, 1);
    expect(v.frame!.layout.issues).toEqual([]);
    const [x, y] = cellCentre(g, 30, 39);
    v.input.move(x, y);
    v.draw();
    const h = v.hover.current;
    expect([h.position, h.yPosition, h.system, h.ySystem]).toEqual([30, 39, '/data/trna', '/data/trna']);
    const all = v.hover.projectedAll()!;
    expect([...all.parts[0].ranges]).toEqual([30, 31, 39, 40]);
    // the crosshair: a column, a row and the cell on the plot, the two bases on the letters (projected outlines)
    const hover = recordDrawLists(v.frame!.lists, ['hover'])
      .split('\n')
      .filter((l) => l.startsWith('rect'));
    // the crosshair (a column, a row, the cell outlined), the two bases on the letters
    expect(hover.length).toBe(3 + 2);
    // the arcs ending at either base lit (the projected hover), and none other
    const arcs = arcsOf(v);
    const d = arcs.prep!.data!;
    const ends = (k: number, p: number) => (d.aStart[k] <= p && p < d.aEnd[k]) || (d.bStart[k] <= p && p < d.bEnd[k]);
    let lit = 0;
    for (let i = 0; i < arcs.n; i++) if (ends(arcs.pair[i], 30) || ends(arcs.pair[i], 39)) lit++;
    expect(lit).toBeGreaterThan(0);
    expect(v.frame!.lists.passes.hover.paths.n).toBe(lit);
  });

  it('the plot’s x is the letters’ x (a hovered cell’s column under its base); its cells square where the room allows, the height capped otherwise', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    const v = viewer;
    const top = v.frame!.layout.sections[0].axis!;
    const st = v.frame!.layout.sections[1].placement.strips;
    const x = st.axes[0];
    expect([x.x0, x.width, x.u0, x.u1]).toEqual([top.x0, top.width, top.u0, top.u1]);
    expectAspect(v, 1);
    // a hovered cell's column is its base's cell on the letters
    const g = gridOf(v, 1);
    const [cx, cy] = cellCentre(g, 40, 47);
    v.input.move(cx, cy);
    v.draw();
    expect(v.hover.current.position).toBe(40);
    top.cell(40);
    expect(cx).toBeCloseTo((top.sx0 + top.sx1) / 2, 6);
  });

  it('a cell: click a one-cell rect, Shift+click the rect from it, Cmd/Ctrl+click toggles a cell out', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    const v = viewer;
    v.input.mac = false;
    const g = gridOf(v, 1);
    const [x, y] = cellCentre(g, 10, 20);
    v.input.click(x, y);
    const rect = (k = 0) => {
      const r = v.selection.value.value.location.rects[k];
      return [...r.a.ranges, ...r.b.ranges];
    };
    expect(rect()).toEqual([10, 11, 20, 21]);
    const [x2, y2] = cellCentre(g, 14, 25);
    v.input.click(x2, y2, { shiftKey: true });
    expect(rect()).toEqual([10, 15, 20, 26]);
    // a cell inside toggled out: the rect less it, exact (at most four rects of single ranges, here two of several)
    const [x3, y3] = cellCentre(g, 12, 22);
    v.input.click(x3, y3, { ctrlKey: true });
    const loc = v.selection.value.value.location;
    let cells = 0;
    for (const r of loc.rects) {
      let w = 0;
      let h = 0;
      for (let k = 0; k < r.a.ranges.length; k += 2) w += r.a.ranges[k + 1] - r.a.ranges[k];
      for (let k = 0; k < r.b.ranges.length; k += 2) h += r.b.ranges[k + 1] - r.b.ranges[k];
      cells += w * h;
    }
    expect(cells).toBe(5 * 6 - 1);
    // the rect's sides light the letters (the projected selection)
    v.draw();
    expect(
      recordDrawLists(v.frame!.lists, ['selection'])
        .split('\n')
        .filter((l) => l.startsWith('rect')).length,
    ).toBeGreaterThan(2);
  });

  it('a rectangle: a drag selects it, Cmd/Ctrl+drag adds one, Alt+drag takes one out', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    const v = viewer;
    v.input.mac = false;
    const g = gridOf(v, 1);
    const drag = (from: [number, number], to: [number, number], mods = {}) => {
      expect(v.input.down(from[0], from[1], mods)).toBe(true);
      v.input.move(to[0], to[1], mods);
      v.input.up();
    };
    drag(cellCentre(g, 5, 5), cellCentre(g, 14, 14));
    let loc = v.selection.value.value.location;
    expect(loc.rects.map((r) => [...r.a.ranges, ...r.b.ranges])).toEqual([[5, 15, 5, 15]]);
    drag(cellCentre(g, 40, 40), cellCentre(g, 44, 44), { ctrlKey: true });
    loc = v.selection.value.value.location;
    expect(loc.rects).toHaveLength(2);
    drag(cellCentre(g, 8, 8), cellCentre(g, 10, 10), { altKey: true });
    loc = v.selection.value.value.location;
    let cells = 0;
    for (const r of loc.rects) {
      let w = 0;
      let h = 0;
      for (let k = 0; k < r.a.ranges.length; k += 2) w += r.a.ranges[k + 1] - r.a.ranges[k];
      for (let k = 0; k < r.b.ranges.length; k += 2) h += r.b.ranges[k + 1] - r.b.ranges[k];
      cells += w * h;
    }
    expect(cells).toBe(100 + 25 - 9);
    expect(loc.rects.length).toBeLessThanOrEqual(1 + 4);
  });

  it('a zoom on the plot moves its leader: x and y both, about the pointer', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    const v = viewer;
    v.input.mac = false;
    const g = gridOf(v, 1);
    const [x, y] = cellCentre(g, 38, 38);
    v.input.wheel(x, y, { deltaX: 0, deltaY: -300, ctrlKey: true, deltaMode: 0 });
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    const w = v.view.get('/views/main/trna')!;
    expect(w.end - w.start).toBeLessThan(76);
    const strips = v.frame!.layout.sections[1].placement.strips;
    expect([strips.ys[0]!.u0, strips.ys[0]!.u1]).toEqual([strips.axes[0].u0, strips.axes[0].u1]);
  });
});

describe('a double-click on a grid whose y follows x', () => {
  it('tRNA’s dot plot: the step keeps the point under the pointer on x, y following x’s window', async () => {
    viewer = await openExample(trna('dot-plot'), { gene: '' });
    const v = viewer;
    const range = () => {
      const r = v.view.get('/views/main/trna')!;
      return [r.start, r.end];
    };
    const [s0, e0] = range();
    const g = gridOf(v, 1);
    // off the diagonal: x 10, y 60
    const [x, y] = cellCentre(g, 10, 60);
    const at = v.placementOf(v.frame!.layout.sections[1].section.address)!.offsetAt(x, y);
    const fx = (at - s0) / (e0 - s0);
    expect(v.input.doubleClick(x, y)).toBe(true);
    await vi.advanceTimersByTimeAsync(800);
    v.draw();
    const [s1, e1] = range();
    expect(e1 - s1).toBeCloseTo((e0 - s0) / v.config.navigate.doubleClickZoomStep, 6);
    expect((at - s1) / (e1 - s1)).toBeCloseTo(fx, 6);
    const st = v.frame!.layout.sections[1].placement.strips;
    expect([st.ys[0]!.u0, st.ys[0]!.u1]).toEqual([st.axes[0].u0, st.axes[0].u1]);
  });

  it('over the complex: a double-click in the gap between blocks zooms about that point as on a cell, never fits', async () => {
    viewer = await openExample(complex('complex-grid'), { fixtures: true });
    const v = viewer;
    const range = () => {
      const r = v.view.get('/views/main/chains')!;
      return [r.start, r.end];
    };
    const [s0, e0] = range();
    for (const step of [1, 2]) {
      v.draw();
      const p = v.placementOf('/views/main/pae') as GridPlacement;
      const xr = p.xRow!;
      // the gap before the second block (x), at the first y block's middle
      const x = (xr.x[0] + xr.w[0] + xr.x[1]) / 2;
      const y = p.yRow!.x[0] + p.yRow!.w[0] / 2;
      expect(xr.x[1] - (xr.x[0] + xr.w[0])).toBeGreaterThan(2);
      const at = p.offsetAt(x, y);
      const [a, b] = range();
      const f = (at - a) / (b - a);
      expect(v.input.doubleClick(x, y)).toBe(true);
      await vi.advanceTimersByTimeAsync(800);
      const [c, d] = range();
      expect(d - c, `step ${step}`).toBeCloseTo((e0 - s0) / v.config.navigate.doubleClickZoomStep ** step, 6);
      expect((at - c) / (d - c)).toBeCloseTo(f, 6);
    }
  });
});

describe('Complex: PAE and contacts, headless', () => {
  it('the grid’s x follows CDK2’s section, its y cyclin A’s; a cell hovered lights a residue of each', async () => {
    viewer = await openExample(complex('pae'), { fixtures: true });
    const v = viewer;
    expect(v.frame!.layout.issues).toEqual([]);
    v.view.setRange('/views/main/cdk2', 40, 80);
    v.view.setRange('/views/main/cyclin', 100, 140);
    await vi.advanceTimersByTimeAsync(500);
    v.draw();
    const strips = v.frame!.layout.sections[2].placement.strips;
    expect([strips.axes[0].u0, strips.axes[0].u1]).toEqual([40, 80]);
    expect([strips.ys[0]!.u0, strips.ys[0]!.u1]).toEqual([100, 140]);
    const g = gridOf(v, 2);
    const [x, y] = cellCentre(g, 50, 120);
    v.input.move(x, y);
    v.draw();
    const all = v.hover.projectedAll()!;
    expect(all.parts.map((p) => [p.system, ...p.ranges])).toEqual([
      ['/data/cdk2-protein', 50, 51],
      ['/data/cyclin-protein', 120, 121],
    ]);
    const svg = v.renderSvg({ idPrefix: 'pae-' });
    await vi.advanceTimersByTimeAsync(400);
    expect((await svg).issues).toEqual([]);
  });

  it('the PAE’s x is CDK2’s x whatever the two windows; its cells square where the room allows, capped otherwise', async () => {
    viewer = await openExample(complex('pae'), { fixtures: true });
    const v = viewer;
    for (const [c, y] of [
      [80, 40],
      [40, 120],
      [298, 432],
    ]) {
      v.view.setRange('/views/main/cdk2', 0, c);
      v.view.setRange('/views/main/cyclin', 0, y);
      await vi.advanceTimersByTimeAsync(500);
      v.draw();
      const cdk2 = v.frame!.layout.sections[0].axis!;
      const x = v.frame!.layout.sections[2].placement.strips.axes[0];
      expect([x.x0, x.width, x.u0, x.u1]).toEqual([cdk2.x0, cdk2.width, cdk2.u0, cdk2.u1]);
      expectAspect(v, 2);
    }
  });

  it('a zoom on the plot moves both leaders about the pointer; a vertical wheel pans y’s leader (cyclin A)', async () => {
    viewer = await openExample(complex('pae'), { fixtures: true });
    const v = viewer;
    v.view.setRange('/views/main/cdk2', 40, 120);
    v.view.setRange('/views/main/cyclin', 100, 180);
    await vi.advanceTimersByTimeAsync(500);
    v.draw();
    const range = (s: string) => {
      const r = v.view.get(`/views/main/${s}`)!;
      return [r.start, r.end];
    };
    const g = gridOf(v, 2);
    const [x, y] = cellCentre(g, 60, 150);
    // a pinch (Ctrl+wheel) in: both windows narrower, each about the pointer's position on it
    v.input.wheel(x, y, { deltaX: 0, deltaY: -40, ctrlKey: true, deltaMode: 0 });
    await vi.advanceTimersByTimeAsync(100);
    const [c0, c1] = range('cdk2');
    const [y0, y1] = range('cyclin');
    expect(c1 - c0).toBeLessThan(80);
    expect(y1 - y0).toBeLessThan(80);
    expect(c0 < 60.5 && 60.5 < c1).toBe(true);
    expect(y0 < 150.5 && 150.5 < y1).toBe(true);
    // the anchor stays where it was: the pointer's fraction of each window
    expect((150.5 - y0) / (y1 - y0)).toBeCloseTo((150.5 - 100) / 80, 1);
    v.draw();
    // a vertical wheel pans cyclin A, x stays
    v.input.wheel(x, y, { deltaX: 0, deltaY: 60, deltaMode: 0 });
    await vi.advanceTimersByTimeAsync(100);
    expect(range('cdk2')).toEqual([c0, c1]);
    expect(range('cyclin')[0]).toBeGreaterThan(y0);
  });

  it('a plot’s box is fixed: a zoom of one axis (y following cyclin A’s section) leaves its height and every track’s place', async () => {
    viewer = await openExample(complex('pae'), { fixtures: true });
    const v = viewer;
    const places = () => {
      v.draw();
      return v.frame!.layout.sections.map((s) => [
        s.placement.box.y,
        s.placement.box.height,
        ...s.tracks.map((t) => t.box.y),
      ]);
    };
    const before = places();
    const height = v.height;
    for (const [s, e] of [
      [0, 40],
      [380, 432],
      [100, 300],
    ]) {
      v.view.setRange('/views/main/cyclin', s, e);
      await vi.advanceTimersByTimeAsync(500);
      expect(places()).toEqual(before);
      expect(v.height).toBe(height);
    }
  });

  it('the wheel over the plot is the plot’s: it pans y to its end and stops there, consumed (no page scroll); outside it, the page’s', async () => {
    viewer = await openExample(complex('pae'), { fixtures: true });
    const v = viewer;
    v.input.mac = false;
    v.view.setRange('/views/main/cdk2', 40, 120);
    v.view.setRange('/views/main/cyclin', 100, 180);
    await vi.advanceTimersByTimeAsync(500);
    v.draw();
    const box = { ...v.placementOf('/views/main/pae')!.box };
    const [x, y] = [v.dataX + 100, box.y + box.height / 2];
    const cyclin = () => {
      const r = v.view.get('/views/main/cyclin')!;
      return [r.start, r.end];
    };
    // down past y's end: every wheel taken; y at its end, its width kept, the plot where it was
    for (let i = 0; i < 40; i++) {
      expect(v.input.wheel(x, y, { deltaX: 0, deltaY: 120, deltaMode: 0 }), `wheel ${i}`).toBe(true);
      await vi.advanceTimersByTimeAsync(20);
      v.draw();
      expect(v.placementOf('/views/main/pae')!.box).toEqual(box);
    }
    const [s, e] = cyclin();
    expect(e).toBeCloseTo(432, 6);
    expect(e - s).toBeCloseTo(80, 6);
    // at the end, still taken (the page does not scroll from inside the plot), nothing moves
    expect(v.input.wheel(x, y, { deltaX: 0, deltaY: 120, deltaMode: 0 })).toBe(true);
    await vi.advanceTimersByTimeAsync(20);
    expect(cyclin()).toEqual([s, e]);
    // and back up to y's start: the same, the width kept
    for (let i = 0; i < 60; i++) v.input.wheel(x, y, { deltaX: 0, deltaY: -120, deltaMode: 0 });
    await vi.advanceTimersByTimeAsync(20);
    const [s2, e2] = cyclin();
    expect(s2).toBeCloseTo(0, 6);
    expect(e2 - s2).toBeCloseTo(80, 6);
    // below the view's sections (outside the plot and every section): the page's wheel, not taken
    expect(v.input.wheel(x, v.height + 20, { deltaX: 0, deltaY: 120, deltaMode: 0 })).toBe(false);
  });

  it('a double-click on a cell zooms a step about it on both axes, again further (never re-centred); back undoes each, `0` fits', async () => {
    viewer = await openExample(complex('pae'), { fixtures: true });
    const v = viewer;
    v.input.mac = false;
    const range = (s: string) => {
      const r = v.view.get(`/views/main/${s}`)!;
      return [r.start, r.end];
    };
    const fit = [range('cdk2'), range('cyclin')];
    const step = v.config.navigate.doubleClickZoomStep;
    // the fraction of each window a position is at
    const frac = (s: string, u: number) => {
      const [a, b] = range(s);
      return (u - a) / (b - a);
    };
    let g = gridOf(v, 2);
    let [x, y] = cellCentre(g, 100, 300);
    const fx = frac('cdk2', 100.5);
    const fy = frac('cyclin', 300.5);
    expect(v.input.doubleClick(x, y)).toBe(true);
    await vi.advanceTimersByTimeAsync(800);
    v.draw();
    const once = [range('cdk2'), range('cyclin')];
    // a step about the cell on both axes: each window by the step, the cell where it was (their ratio kept)
    expect(once[0][1] - once[0][0]).toBeCloseTo((fit[0][1] - fit[0][0]) / step, 6);
    expect(once[1][1] - once[1][0]).toBeCloseTo((fit[1][1] - fit[1][0]) / step, 6);
    expect(frac('cdk2', 100.5)).toBeCloseTo(fx, 6);
    expect(frac('cyclin', 300.5)).toBeCloseTo(fy, 6);
    // again on another cell: further, about that one
    g = gridOf(v, 2);
    [x, y] = cellCentre(g, 110, 320);
    const gx = frac('cdk2', 110.5);
    const gy = frac('cyclin', 320.5);
    expect(v.input.doubleClick(x, y)).toBe(true);
    await vi.advanceTimersByTimeAsync(800);
    v.draw();
    const twice = [range('cdk2'), range('cyclin')];
    expect(twice[0][1] - twice[0][0]).toBeCloseTo((fit[0][1] - fit[0][0]) / step ** 2, 6);
    expect(twice[1][1] - twice[1][0]).toBeCloseTo((fit[1][1] - fit[1][0]) / step ** 2, 6);
    expect(frac('cdk2', 110.5)).toBeCloseTo(gx, 6);
    expect(frac('cyclin', 320.5)).toBeCloseTo(gy, 6);
    // Alt+Left back: once to the first zoom, again to the fit
    expect(v.input.key({ key: 'ArrowLeft', altKey: true })).toBe(true);
    expect([range('cdk2'), range('cyclin')]).toEqual(once);
    expect(v.input.key({ key: 'ArrowLeft', altKey: true })).toBe(true);
    expect([range('cdk2'), range('cyclin')]).toEqual(fit);
    // forward again, and `0` fits the plot
    v.view.forward();
    v.view.forward();
    expect([range('cdk2'), range('cyclin')]).toEqual(twice);
    await vi.advanceTimersByTimeAsync(800);
    v.draw();
    const [px, py] = [v.dataX + 50, v.placementOf('/views/main/pae')!.box.y + 50];
    v.input.move(px, py);
    expect(v.input.key({ key: '0' })).toBe(true);
    expect([range('cdk2'), range('cyclin')]).toEqual(fit);
  });

  it('the contacts as arcs between the chains’ blocks: a hover lights both ends', async () => {
    viewer = await openExample(complex('contacts'), { fixtures: true });
    const v = viewer;
    // 1FIN's contacts of both copies (a pair of residues once a copy), folded onto the merged blocks
    const lay = arcsOf(v);
    expect(lay.n).toBeGreaterThan(100);
    // an arc no ATP arc (its own layer, over the chains' arcs) lies on at its middle
    const systems = (h: ReturnType<typeof v.hit>) => h?.extent?.parts.map((p) => p.system);
    let at: [number, number] | undefined;
    for (let k = 0; k < lay.n && !at; k++) {
      const [x, y] = curveAt(lay.curve, k, 0.5);
      if (JSON.stringify(systems(v.hit(x, y))) === JSON.stringify(['/data/cdk2', '/data/cyclin'])) at = [x, y];
    }
    expect(at).toBeDefined();
    v.input.move(at![0], at![1]);
    v.draw();
    expect(systems(v.hover.current)).toEqual(['/data/cdk2', '/data/cyclin']);
  });
});
