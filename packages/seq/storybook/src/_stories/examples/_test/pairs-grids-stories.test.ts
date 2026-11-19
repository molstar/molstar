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

import { type GridLay, type PairsLay, recordDrawLists } from '@molstar/seqx-core/render';
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

/** The paths of the data pass, as recorded (the arcs). */
const paths = (lists: Parameters<typeof recordDrawLists>[0]) =>
  recordDrawLists(lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('path'));

describe('tRNA: arcs, headless', () => {
  it('every pair has its arc; the export draws the canvas’s arcs; a base hovered lights its partner’s arc', async () => {
    viewer = await openExample(trna('arcs'), { gene: '' });
    const v = viewer;
    const lay = arcsOf(v);
    expect(lay.n).toBe(27);
    const svg = v.renderSvg({ idPrefix: 't-' });
    await vi.advanceTimersByTimeAsync(400);
    const r = await svg;
    expect(r.issues).toEqual([]);
    expect(paths(r.lists)).toEqual(paths(v.frame!.lists));
    // hover the letters at position 9 (1-based 10): the arc 10–25 lights
    const letters = v.frame!.layout.sections[0].tracks[1].box;
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
    // the inner stem's pair 31–39, then focus: the window its ends' hull
    const j = Array.from({ length: lay.n }, (_, k) => k).find((k) => lay.prep!.data!.keys[lay.pair[k]] === '31-39')!;
    const [x2, y2] = curveAt(lay.curve, 8 * j, 0.5);
    v.input.doubleClick(x2, y2);
    await vi.advanceTimersByTimeAsync(1000);
    v.draw();
    const w = v.view.get('/views/main/trna')!;
    expect(w.start).toBeLessThanOrEqual(30);
    expect(w.end).toBeGreaterThanOrEqual(39);
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
  it('a cell hovered: the crosshair, both bases lit on the letters (projected)', async () => {
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
    // the crosshair (a column, a row, the cell outlined), the two bases on the letters; the pairs at those bases not lit
    expect(hover.length).toBe(3 + 2);
    expect(v.frame!.lists.passes.hover.paths.n).toBe(0);
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

describe('Complex: PAE and contacts, headless', () => {
  it('the grid’s x follows CDK2’s section, its y cyclin A’s; a cell hovered lights a residue of each', async () => {
    viewer = await openExample(complex('pae'), { gene: '' });
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
      ['/data/cdk2', 50, 51],
      ['/data/cyclin', 120, 121],
    ]);
    const svg = v.renderSvg({ idPrefix: 'pae-' });
    await vi.advanceTimersByTimeAsync(400);
    expect((await svg).issues).toEqual([]);
  });

  it('the contacts as arcs between the chains’ blocks: a hover lights both ends', async () => {
    viewer = await openExample(complex('contacts'), { gene: '' });
    const v = viewer;
    const lay = arcsOf(v);
    expect(lay.n).toBe(24);
    const [x, y] = curveAt(lay.curve, 0, 0.5);
    v.input.move(x, y);
    v.draw();
    expect(v.hover.current.extent?.parts.map((p) => p.system)).toEqual(['/data/cdk2', '/data/cyclin']);
  });
});
