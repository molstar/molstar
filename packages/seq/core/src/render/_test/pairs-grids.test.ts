/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Pairs and grids drawn (_spec/README.md "Render", The grid layout, The grid layer, The pairs layer): the grid's y
// following x; cells drawn against the data at the fit and across an image tile's edge, hit as the inverse of
// placement, kept by the triangle; the build budget while moving, the rest when settled; the crosshair and a selected
// rect on y; every pair in view with its arc, every arc hit on its curve, arcs in their band, nested arcs nesting, an end
// past the window taken `reach` past it so an arc moves continuously; `split` against a brute force; `maxArcs`; chords
// on a ring; arcs between blocks and a caret for an end not shown; the pairs as dots, mirrored and by kind.

import { location, rangesLocation } from '@molstar/seqx-core/query';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import {
  crosses,
  DrawLists,
  drawView,
  type Emphasis,
  type GridLay,
  GridPlacement,
  hitView,
  type PairsLay,
  PairsType,
  type SectionLayout,
  splitCrossing,
  type ViewLayout,
  viewLegends,
} from '../index.js';
import { ImageTile } from '../images.js';
import { settle } from './helpers.js';
import { arcsSpec, cellValues, gridSpec, letters, pairRows, randomPairs } from './pairs-fixtures.js';

const Size = { width: 720, height: 900 };

function gridLay(layout: ViewLayout, section = 1, layer = 0): GridLay {
  const s = layout.sections[section];
  return s.lays[0 * s.placement.strips.bands + 0][layer].lay as GridLay;
}

function arcsLay(s: SectionLayout): PairsLay {
  return s.sectionLays[0].lay as PairsLay;
}

/** The point of a cubic at t. */
function at(c: Float64Array, o: number, t: number): [number, number] {
  const u = 1 - t;
  const f = (k: number) =>
    u * u * u * c[o + k] + 3 * u * u * t * c[o + 2 + k] + 3 * u * t * t * c[o + 4 + k] + t * t * t * c[o + 6 + k];
  return [f(0), f(1)];
}

describe('the grid section', () => {
  it('places one strip, the plot, with y following x’s window: square cells at aspect 1', async () => {
    const s = await settle(gridSpec({ n: 40 }));
    const f = await s.renderer.frame({
      size: Size,
      viewports: { '/views/v/s': { start: 10, end: 30 } },
      settled: true,
    });
    const g = f.layout.sections[1];
    expect(g.placement).toBeInstanceOf(GridPlacement);
    const strips = g.placement.strips;
    const x = strips.axes[0];
    const y = strips.ys[0]!;
    expect([x.u0, x.u1, y.u0, y.u1]).toEqual([10, 30, 10, 30]);
    expect(y.x0).toBe(strips.top[0]);
    expect(y.width).toBeCloseTo(x.width);
    expect(g.placement.box.height).toBeCloseTo(x.width);
    expect(f.layout.issues).toEqual([]);
  });

  it('a y of another system keeps a window of its own (`viewport.y`), all of it by default', async () => {
    const s = await settle(gridSpec({ n: 30, ny: 12 }));
    let f = await s.renderer.frame({ size: Size, settled: true });
    let y = f.layout.sections[1].placement.strips.ys[0]!;
    expect([y.system, y.u0, y.u1]).toEqual(['/data/y', 0, 12]);
    f = await s.renderer.frame({
      size: Size,
      settled: true,
      viewports: { '/views/v/plot': { start: 0, end: 30, y: { start: 4, end: 8 } } },
    });
    y = f.layout.sections[1].placement.strips.ys[0]!;
    expect([y.u0, y.u1]).toEqual([4, 8]);
  });
});

describe('a grid across two sections', () => {
  it('x follows one section, y another: a pan of y’s leader alone lays the grid out again with y moved', async () => {
    const b = createSeqSpec();
    const x = b.data.sequence('x', { alphabet: 'rna', letters: letters(40) });
    const y = b.data.sequence('y', { alphabet: 'rna', letters: letters(30, 5) });
    const g = b.data.grid('g', { x, y, values: cellValues(40, 30) });
    b.views.view('v', (v) => {
      const sx = v.section('x', { of: x });
      sx.track('t').layer('l', { type: 'letters' });
      const sy = v.section('y', { of: y });
      sy.track('t').layer('l', { type: 'letters' });
      v.section('plot', { layout: 'grid', of: x, of_y: y, follow: sx.address, follow_y: sy.address })
        .track('t')
        .layer('g', { type: 'grid', data: g });
    });
    const st = await settle(b.build());
    const vx = { start: 5, end: 25 };
    let f = await st.renderer.frame({
      size: Size,
      settled: true,
      viewports: { '/views/v/x': vx, '/views/v/y': { start: 0, end: 30 } },
    });
    let plot = f.layout.sections[2].placement.strips;
    expect([plot.axes[0].u0, plot.axes[0].u1]).toEqual([5, 25]);
    expect([plot.ys[0]!.u0, plot.ys[0]!.u1]).toEqual([0, 30]);
    f = await st.renderer.frame({
      size: Size,
      settled: true,
      viewports: { '/views/v/x': vx, '/views/v/y': { start: 10, end: 20 } },
    });
    plot = f.layout.sections[2].placement.strips;
    expect([plot.ys[0]!.u0, plot.ys[0]!.u1]).toEqual([10, 20]);
    expect(f.layout.issues).toEqual([]);
  });
});

describe('a plot’s box is fixed', () => {
  it('a zoom of one axis (y following a section) leaves the plot’s height and every track’s place; cells rectangles', async () => {
    const b = createSeqSpec();
    const x = b.data.sequence('x', { alphabet: 'rna', letters: letters(40) });
    const y = b.data.sequence('y', { alphabet: 'rna', letters: letters(30, 5) });
    const g = b.data.grid('g', { x, y, values: cellValues(40, 30) });
    b.views.view('v', (v) => {
      const sx = v.section('x', { of: x });
      sx.track('t').layer('l', { type: 'letters' });
      const sy = v.section('y', { of: y });
      sy.track('t').layer('l', { type: 'letters' });
      v.section('plot', { layout: 'grid', of: x, of_y: y, follow: sx.address, follow_y: sy.address })
        .track('t')
        .layer('g', { type: 'grid', data: g });
      v.section('below', { of: x }).track('t').layer('l', { type: 'letters' });
    });
    const st = await settle(b.build());
    const places = (l: ViewLayout) =>
      l.sections.map((s) => [s.placement.box.y, s.placement.box.height, ...s.tracks.map((t) => t.box.y)]);
    const frame = (vy: { start: number; end: number }) =>
      st.renderer.frame({
        size: Size,
        settled: true,
        viewports: { '/views/v/x': { start: 0, end: 40 }, '/views/v/y': vy },
      });
    const fit = await frame({ start: 0, end: 30 });
    const plot = fit.layout.sections[2].placement;
    // the data box's width over the aspect (1): square as a box, whatever the windows
    expect(plot.box.height).toBeCloseTo(plot.strips.axes[0].width, 6);
    const before = places(fit.layout);
    for (const vy of [
      { start: 10, end: 20 },
      { start: 0, end: 4 },
      { start: 26, end: 30 },
    ]) {
      const f = await frame(vy);
      expect(places(f.layout)).toEqual(before);
      const s = f.layout.sections[2].placement.strips;
      expect([s.ys[0]!.u0, s.ys[0]!.u1]).toEqual([vy.start, vy.end]);
      // the cells rectangles: x's 40 across the width, y's few down the same height
      expect(s.ys[0]!.scale).toBeGreaterThan(s.axes[0].scale);
    }
  });
});

describe('the grid layer', () => {
  it('at the fit and zoomed in: every cell in view drawn with its value’s color, hit back as its two positions', async () => {
    const n = 24;
    const s = await settle(gridSpec({ n }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const lay = gridLay(f.layout);
    // 24 cells across ~600 px: rects
    expect(lay.images.length).toBe(0);
    expect(lay.n).toBe(n * n);
    const a = lay.axis!;
    const ay = lay.axisY!;
    for (let i = 0; i < n; i += 5)
      for (let j = 0; j < n; j += 7) {
        a.cell(i);
        const cx = (a.sx0 + a.sx1) / 2;
        ay.cell(j);
        const cy = (ay.sx0 + ay.sx1) / 2;
        const h = hitView(f.layout, cx, cy);
        expect(h).toMatchObject({ target: 'cell', position: i, yPosition: j, system: '/data/s', ySystem: '/data/s' });
        // the rect there has the color of i * 1000 + j: values rise along the grid, so colors are distinct and ordered
        const k = i * n + j;
        expect(lay.rx[4 * k]).toBeCloseTo(Math.min(a.sx0, a.sx1));
        expect(lay.rx[4 * k + 1]).toBeCloseTo(Math.min(ay.sx0, ay.sx1));
      }
    expect(lay.rc[0]).not.toBe(lay.rc[n * n - 1]);
  });

  it('keeps to a triangle: upper is y before x; a missing value is transparent', async () => {
    const n = 10;
    const values = cellValues(n, n);
    values[3 * n + 1] = Number.NaN;
    const s = await settle(gridSpec({ n, values, triangle: 'upper' }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const lay = gridLay(f.layout);
    // the upper triangle less one missing cell
    expect(lay.n).toBe((n * (n - 1)) / 2 - 1);
    const a = lay.axis!;
    const ay = lay.axisY!;
    a.cell(2);
    const x = (a.sx0 + a.sx1) / 2;
    // a cell not drawn is still the plot's cell under the point (the crosshair hovers the whole plot)
    ay.cell(5);
    expect(hitView(f.layout, x, (ay.sx0 + ay.sx1) / 2)).toMatchObject({ position: 2, yPosition: 5 });
    ay.cell(1);
    expect(hitView(f.layout, x, (ay.sx0 + ay.sx1) / 2)).toMatchObject({ position: 2, yPosition: 1 });
  });

  it('below rectMinPx: image tiles of the shared reduction, across a tile edge with no gap, each sample its cell’s color', async () => {
    const n = 500;
    const s = await settle(gridSpec({ n }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const lay = gridLay(f.layout);
    expect(lay.n).toBe(0);
    // ~584 px for 500 cells at dpr 1: level 0, two tiles a side
    expect(lay.level).toBe(0);
    expect(lay.images.length).toBe(4);
    const a = lay.axis!;
    // tiles meet: the first tile's right edge is the second's left (the 256th cell's boundary)
    const xs = new Set<number>();
    for (let k = 0; k < lay.images.length; k++) xs.add(Math.round(lay.ib[4 * k] * 1000) / 1000);
    expect([...xs].sort((p, q) => p - q)).toEqual([0, ImageTile].map((p) => Math.round(a.x(p))));
    // on the device grid at any ratio, so neighbours share a whole pixel edge (no antialiased seam between tiles)
    for (const dpr of [1, 1.5, 2]) {
      const g = gridLay((await s.renderer.frame({ size: Size, dpr, settled: true })).layout);
      for (let k = 0; k < 4 * g.images.length; k++)
        expect(Number.isInteger(Math.round(g.ib[k] * dpr * 1e6) / 1e6)).toBe(true);
      const right = g.ib[0] + g.ib[2];
      const next = g.images.findIndex((_, k) => g.ib[4 * k] > g.ib[0] && g.ib[4 * k + 1] === g.ib[1]);
      expect(g.ib[4 * next]).toBeCloseTo(right, 9);
    }
    // across the edge: the samples at x 255 and 256 (y 0) are their cells' colors
    const first = lay.images[0];
    const second = lay.images.find((_, k) => lay.ib[4 * k] > lay.ib[0] && lay.ib[4 * k + 1] === lay.ib[1])!;
    const colorAt = (img: typeof first, c: number, r: number) => {
      const q = 4 * (r * img.width + c);
      return ((img.data[q] << 24) | (img.data[q + 1] << 16) | (img.data[q + 2] << 8) | img.data[q + 3]) >>> 0;
    };
    lay.datum.value = 255 * 1000;
    const c255 = lay.color!.at(lay.datum);
    lay.datum.value = 256 * 1000;
    const c256 = lay.color!.at(lay.datum);
    expect(colorAt(first, 255, 0)).toBe(c255);
    expect(colorAt(second, 0, 0)).toBe(c256);
  });

  it('zoomed out past a device pixel a cell: a coarser level, its samples by `aggregate`', async () => {
    const n = 2000;
    const s = await settle(gridSpec({ n, aggregate: 'min' }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const lay = gridLay(f.layout);
    // ~600 px for 2000 cells: a sample is 4 cells a side (level 2)
    expect(lay.level).toBe(2);
    const img = lay.images[0];
    // the sample (0, 1) holds cells x 0–3, y 4–7: their least is x 0, y 4
    lay.datum.value = 4;
    const want = lay.color!.at(lay.datum);
    const q = 4 * (1 * img.width + 0);
    expect(((img.data[q] << 24) | (img.data[q + 1] << 16) | (img.data[q + 2] << 8) | img.data[q + 3]) >>> 0).toBe(want);
  });

  it('a triangle reduces its own cells alone: coarse samples on the diagonal drawn, none mixing in the other half', async () => {
    const n = 1024;
    const values = (own: number, other: number, upper: boolean) => {
      const v = new Float32Array(n * n);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) v[i * n + j] = (upper ? j < i : j > i) ? own : other;
      return v;
    };
    const b = createSeqSpec();
    const sq = b.data.sequence('s', { alphabet: 'rna', letters: letters(n) });
    const ga = b.data.grid('a', { x: sq, values: values(1, 100, true) });
    const gb = b.data.grid('b', { x: sq, values: values(2, 200, false) });
    b.views.view('v', (v) => {
      v.section('plot', { layout: 'grid', of: sq }).track('plot', {}, (t) => {
        t.layer('upper', { type: 'grid', data: ga, triangle: 'upper', aggregate: 'mean' });
        t.layer('lower', { type: 'grid', data: gb, triangle: 'lower', aggregate: 'mean' });
      });
    });
    const s = await settle(b.build());
    const f = await s.renderer.frame({ size: Size, settled: true });
    for (const k of [0, 1]) {
      const lay = gridLay(f.layout, 0, k);
      expect(lay.level).toBeGreaterThan(0);
      let first = -1;
      let diagonal = 0;
      let mixed = 0; // counted, asserted once: an `expect` a pixel costs more than the pixel
      for (let t = 0; t < lay.images.length; t++) {
        const img = lay.images[t];
        // a tile on the diagonal: as far right of the plot's left as down from its top
        const onDiagonal = Math.abs(lay.ib[4 * t] - lay.box.x - (lay.ib[4 * t + 1] - lay.box.y)) < 1;
        for (let r = 0; r < img.height; r++)
          for (let c = 0; c < img.width; c++) {
            const q = 4 * (r * img.width + c);
            if (img.data[q + 3] === 0) continue;
            const rgba =
              ((img.data[q] << 24) | (img.data[q + 1] << 16) | (img.data[q + 2] << 8) | img.data[q + 3]) >>> 0;
            // every sample its own half's one value: one color
            if (first < 0) first = rgba;
            if (rgba !== first) mixed++;
            if (onDiagonal && r === c) diagonal++;
          }
      }
      expect(mixed).toBe(0);
      // the diagonal's samples hold cells of both halves: each layer draws them
      expect(diagonal).toBeGreaterThan(0);
    }
  });

  it('builds at most `grid.buildBudget` samples a frame while moving (a coarser tile meanwhile, never empty); the rest when settled', async () => {
    const n = 1200;
    const s = await settle(gridSpec({ n }));
    const images = s.renderer.images;
    const before = images.size;
    const moving = await s.renderer.frame({ size: Size, settled: false });
    const lay = gridLay(moving.layout);
    // a tile (65,536 samples) within the budget, the others drawn from a coarser level
    expect(images.size - before).toBeLessThanOrEqual(1 + 4);
    expect(lay.images.length).toBeGreaterThan(0);
    const settled = await s.renderer.frame({ size: Size, settled: true });
    const all = gridLay(settled.layout);
    expect(all.images.every((img) => img.width <= ImageTile && img.height <= ImageTile)).toBe(true);
    // settled: the level's every tile in view
    const side = Math.ceil(n / 2 ** all.level / ImageTile);
    expect(all.images.length).toBe(side * side);
  });

  it('a grid on other systems than its section’s is not drawn, and says so', async () => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'rna', letters: letters(8) });
    const t = b.data.sequence('t', { alphabet: 'rna', letters: letters(8, 3) });
    const g = b.data.grid('g', { x: s, values: cellValues(8, 8) });
    b.views.view('v', (v) => {
      v.section('plot', { layout: 'grid', of: s, of_y: t, follow_y: 'none' })
        .track('t')
        .layer('g', { type: 'grid', data: g });
      v.section('tracks', { of: s }).track('t').layer('g', { type: 'grid', data: g });
    });
    const st = await settle(b.build());
    const f = await st.renderer.frame({ size: Size, settled: true });
    expect(f.layout.issues.map((i) => i.code).sort()).toEqual(['seq.render.grid.systems', 'seq.render.layer.axes']);
  });
});

describe('the emphasis along y', () => {
  it('a hover on a cell: a crosshair (a column, a row) and the cell outlined; a selected rect drawn on the plot', async () => {
    const n = 20;
    const s = await settle(gridSpec({ n }));
    const f0 = await s.renderer.frame({ size: Size, settled: true });
    const lay = gridLay(f0.layout);
    const a = lay.axis!;
    const ay = lay.axisY!;
    const hover = { layer: lay.key, section: 'plot', system: '/data/s', position: 4, ySystem: '/data/s', yPosition: 9 };
    const lists = new DrawLists();
    drawView(f0.layout, lists, { hover } satisfies Emphasis, ['hover']);
    const r = lists.passes.hover.rects;
    const box = lay.box;
    a.cell(4);
    ay.cell(9);
    const rects: number[][] = [];
    for (let k = 0; k < r.n; k++) rects.push([r.x[k], r.y[k], r.w[k], r.h[k]]);
    // the column down the plot at x 4, the row across it at y 9, the cell
    expect(rects.some((q) => Math.abs(q[0] - a.sx0) < 1e-3 && Math.abs(q[3] - box.height) < 1e-3)).toBe(true);
    expect(rects.some((q) => Math.abs(q[1] - ay.sx0) < 1e-3 && Math.abs(q[2] - box.width) < 1e-3)).toBe(true);
    expect(rects.some((q) => Math.abs(q[0] - a.sx0) < 1e-3 && Math.abs(q[1] - ay.sx0) < 1e-3)).toBe(true);
    // a selected rect x 2–5 × y 10–13
    const sel = location({
      rects: [{ a: { system: '/data/s', ranges: [2, 6] }, b: { system: '/data/s', ranges: [10, 14] } }],
    });
    drawView(f0.layout, lists, { selection: sel }, ['selection']);
    const sr = lists.passes.selection.rects;
    let found = false;
    for (let k = 0; k < sr.n; k++)
      if (
        Math.abs(sr.x[k] - a.x(2)) < 1e-3 &&
        Math.abs(sr.y[k] - ay.x(10)) < 1e-3 &&
        Math.abs(sr.w[k] - (a.x(6) - a.x(2))) < 1e-3
      )
        found = true;
    expect(found).toBe(true);
    // a range selected on the tracks: its rows on the plot, beside the strip's columns
    drawView(f0.layout, lists, { selection: rangesLocation('/data/s', [3, 5]) }, ['selection']);
    const rows = lists.passes.selection.rects;
    let row = false;
    for (let k = 0; k < rows.n; k++)
      if (Math.abs(rows.y[k] - ay.x(3)) < 1e-3 && Math.abs(rows.w[k] - box.width) < 1e-3) row = true;
    expect(row).toBe(true);
  });

  it('a selected cell is drawn on the plot as its rect alone: its sides, projected onto the sections of its systems, are not painted back onto it', async () => {
    const s = await settle(gridSpec({ n: 60 }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const lay = gridLay(f.layout);
    const b = lay.box;
    const selection = location({
      rects: [{ a: { system: '/data/s', ranges: [10, 11] }, b: { system: '/data/s', ranges: [40, 41] } }],
    });
    // what the viewer passes: the rect's sides projected (the tracks section lights 10 and 40)
    const projected = { selection: rangesLocation('/data/s', [10, 11, 40, 41]) };
    const lists = new DrawLists();
    drawView(f.layout, lists, { selection, projected }, ['selection']);
    const r = lists.passes.selection.rects;
    let onPlot = 0;
    let onTracks = 0;
    for (let k = 0; k < r.n; k++) {
      const inside = r.y[k] >= b.y - 1e-3 && r.y[k] + r.h[k] <= b.y + b.height + 1e-3;
      if (inside) onPlot++;
      else onTracks++;
    }
    expect(onPlot).toBe(1);
    expect(onTracks).toBeGreaterThan(0);
  });
});

/** The tRNA's shape: four stems, nested, and pairs crossing them. */
const Stems: [number, number][] = [
  [0, 71],
  [1, 70],
  [2, 69],
  [3, 68],
  [4, 67],
  [5, 66],
  [6, 65],
  [9, 24],
  [10, 23],
  [11, 22],
  [12, 21],
  [26, 42],
  [27, 41],
  [28, 40],
  [29, 39],
  [30, 38],
  [48, 64],
  [49, 63],
  [50, 62],
  [51, 61],
  [52, 60],
];
const Crossing: [number, number][] = [
  [7, 13],
  [14, 47],
  [17, 54],
  [18, 55],
  [25, 43],
  [53, 57],
];

describe('the pairs layer: arcs', () => {
  it('every pair in view has its arc in its band, hit at its apex as the pair; hovering a base lights its partner', async () => {
    const s = await settle(arcsSpec({ pairs: [...Stems, ...Crossing] }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const sec = f.layout.sections[0];
    const lay = arcsLay(sec);
    expect(lay.n).toBe(27);
    const strips = sec.placement.strips;
    const top = strips.bandY[1];
    const bottom = top + strips.bandH[1];
    for (let i = 0; i < lay.n; i++) {
      const [x, y] = at(lay.curve, 8 * i, 0.5);
      // the whole curve in the band, its feet on the band's bottom edge (inset by the stroke), its apex above them
      for (let t = 0; t <= 1; t += 0.125) {
        const q = at(lay.curve, 8 * i, t)[1];
        expect(q).toBeGreaterThanOrEqual(top);
        expect(q).toBeLessThanOrEqual(bottom - 1 + 1e-6);
      }
      expect(lay.curve[8 * i + 1]).toBeCloseTo(bottom - 1);
      expect(y).toBeLessThan(bottom - 1);
      const h = hitView(f.layout, x, y);
      expect(h?.element?.key).toBe(`${lay.prep!.data!.aStart[lay.pair[i]]}-${lay.prep!.data!.bStart[lay.pair[i]]}`);
      expect(h?.extent?.parts[0].ranges.length).toBe(4);
    }
    // a hover at a base lights the arcs with an end there (its partner's)
    const lists = new DrawLists();
    drawView(f.layout, lists, { hover: { layer: 'x', section: 's', system: '/data/s', position: 9 } }, ['hover']);
    expect(lists.passes.hover.paths.n).toBe(1);
  });

  it('through the themes and the data filter: kinds colored by the channel’s default, a uniform color, a filter leaving pairs out', async () => {
    const pairs: [number, number][] = [
      [0, 20],
      [2, 18],
      [5, 9],
      [11, 15],
    ];
    const types = ['stem', 'stem', 'tertiary', 'wobble'];
    const build = (query = '', color?: string) => {
      const b = createSeqSpec();
      const sq = b.data.sequence('s', { alphabet: 'rna', letters: letters(24) });
      const p = b.data.pairs('p', { a: sq, rows: pairRows(pairs, types) });
      b.views.view('v', (v) => {
        const sec = v.section('s', { of: sq });
        sec.track('letters').layer('letters', { type: 'letters' });
        sec.track('pairs').layer('arcs', { type: 'pairs', data: `${p}${query}`, ...(color ? { color } : {}) });
      });
      return b.build();
    };
    const lay = async (spec: ReturnType<typeof build>) => {
      const st = await settle(spec);
      const f = await st.renderer.frame({ size: Size, settled: true });
      return { f, lay: arcsLay(f.layout.sections[0]) };
    };
    // by kind: one color a kind, the same kind the same
    const all = (await lay(build())).lay;
    expect(all.n).toBe(4);
    const colorOf = (l: PairsLay, k: number) => l.color[[...l.pair.slice(0, l.n)].indexOf(k)];
    expect(colorOf(all, 0)).toBe(colorOf(all, 1));
    expect(new Set([colorOf(all, 0), colorOf(all, 2), colorOf(all, 3)]).size).toBe(3);
    // a uniform color
    const red = (await lay(build('', '#ff0000'))).lay;
    for (let i = 0; i < red.n; i++) expect(red.color[i]).toBe(0xff0000ff);
    // the filter: the tertiary pair left out, counted
    const { f, lay: some } = await lay(build('/(not :of (item :type tertiary))'));
    expect(some.n).toBe(3);
    expect([...some.pair.slice(0, some.n)].sort()).toEqual([0, 1, 3]);
    expect(f.layout.issues).toEqual([]);
    expect(PairsType.filtered!(some)).toBe(1);
    expect(f.layout.sections[0].tracks[1].filter?.count).toBe(1);
  });

  it('a cell’s hover lights the arcs ending at either position; a hover projected here lights them too', async () => {
    const pairs: [number, number][] = [
      [0, 20],
      [2, 18],
      [5, 9],
      [11, 15],
    ];
    const st = await settle(arcsSpec({ n: 24, pairs }));
    const f = await st.renderer.frame({ size: Size, settled: true });
    const lists = new DrawLists();
    const hover = { layer: 'x', section: 'plot', system: '/data/s', position: 2, ySystem: '/data/s', yPosition: 15 };
    drawView(f.layout, lists, { hover }, ['hover']);
    expect(lists.passes.hover.paths.n).toBe(2);
    drawView(f.layout, lists, { projected: { hover: rangesLocation('/data/s', [9, 10]), hoverSection: '/x' } }, [
      'hover',
    ]);
    expect(lists.passes.hover.paths.n).toBe(1);
  });

  it('nested arcs nest: an outer pair’s apex above an inner one’s (random nested sets)', async () => {
    for (let seed = 1; seed <= 6; seed++) {
      const pairs = randomPairs(60, 12, seed);
      const kept = splitCrossing({
        n: pairs.length,
        aStart: Int32Array.from(pairs.map((p) => p[0])),
        aEnd: Int32Array.from(pairs.map((p) => p[0] + 1)),
        bStart: Int32Array.from(pairs.map((p) => p[1])),
        bEnd: Int32Array.from(pairs.map((p) => p[1] + 1)),
      } as never);
      const nested = pairs.filter((_, k) => kept[k] === 0);
      const s = await settle(arcsSpec({ n: 60, pairs: nested }));
      const f = await s.renderer.frame({ size: Size, settled: true });
      const lay = arcsLay(f.layout.sections[0]);
      expect(lay.n).toBe(nested.length);
      const apex = new Map<number, number>();
      for (let i = 0; i < lay.n; i++) apex.set(lay.pair[i], at(lay.curve, 8 * i, 0.5)[1]);
      for (let p = 0; p < nested.length; p++)
        for (let q = 0; q < nested.length; q++)
          if (nested[p][0] < nested[q][0] && nested[q][1] < nested[p][1])
            expect(apex.get(p)!).toBeLessThan(apex.get(q)!);
    }
  });

  it('an end past the window is taken `reach` past its edge: an arc moves continuously as the window pans', async () => {
    const s = await settle(arcsSpec({ n: 200, pairs: [[20, 180]] }));
    let last: number | undefined;
    for (let start = 0; start <= 60; start += 2) {
      const f = await s.renderer.frame({
        size: Size,
        settled: true,
        viewports: { '/views/v/s': { start, end: start + 50 } },
      });
      const lay = arcsLay(f.layout.sections[0]);
      // both ends past the window on either side: still drawn, each a quarter window past its edge at most
      expect(lay.n).toBe(1);
      const x = lay.curve[0];
      if (last !== undefined) expect(Math.abs(x - last)).toBeLessThan(40);
      last = x;
      // the far end: a quarter window past the right edge
      const a = f.layout.sections[0].axis!;
      expect(lay.curve[6]).toBeCloseTo(a.xOfUnit(a.u1 + 0.25 * (a.u1 - a.u0)));
    }
  });

  it('`split`: the largest nested set on its side, the rest on the other; exact against a brute force', () => {
    for (let seed = 1; seed <= 25; seed++) {
      const pairs = randomPairs(30, 9, 100 + seed);
      const d = {
        n: pairs.length,
        aStart: Int32Array.from(pairs.map((p) => p[0])),
        aEnd: Int32Array.from(pairs.map((p) => p[0] + 1)),
        bStart: Int32Array.from(pairs.map((p) => p[1])),
        bEnd: Int32Array.from(pairs.map((p) => p[1] + 1)),
      } as never;
      const other = splitCrossing(d);
      const kept = pairs.filter((_, k) => !other[k]);
      for (const p of kept) for (const q of kept) expect(crosses(p[0], p[1], q[0], q[1])).toBe(false);
      let best = 0;
      for (let mask = 0; mask < 1 << pairs.length; mask++) {
        const set = pairs.filter((_, k) => mask & (1 << k));
        if (set.length <= best) continue;
        if (set.every((p) => set.every((q) => !crosses(p[0], p[1], q[0], q[1])))) best = set.length;
      }
      expect(kept.length).toBe(best);
    }
  });

  it('`split` keeps self-pairs and equal or end-sharing pairs on their side (they cross nothing): exact against a brute force', () => {
    let seed = 7;
    const next = () => {
      seed = (Math.imul(seed, 1103515245) + 12345) >>> 0;
      return seed >>> 16;
    };
    for (let round = 0; round < 40; round++) {
      // few positions, so self-pairs, duplicates and shared ends are common
      const pairs: [number, number][] = [];
      for (let k = 0; k < 9; k++) {
        const a = next() % 8;
        pairs.push([a, k % 3 === 0 ? a : next() % 8]);
      }
      const d = {
        n: pairs.length,
        aStart: Int32Array.from(pairs.map((p) => p[0])),
        aEnd: Int32Array.from(pairs.map((p) => p[0] + 1)),
        bStart: Int32Array.from(pairs.map((p) => p[1])),
        bEnd: Int32Array.from(pairs.map((p) => p[1] + 1)),
      } as never;
      const other = splitCrossing(d);
      const span = (p: [number, number]) => [Math.min(p[0], p[1]), Math.max(p[0], p[1])] as const;
      const free = (set: [number, number][]) =>
        set.every((p) => set.every((q) => !crosses(span(p)[0], span(p)[1], span(q)[0], span(q)[1])));
      const kept = pairs.filter((_, k) => !other[k]);
      expect(free(kept)).toBe(true);
      for (let k = 0; k < pairs.length; k++) if (pairs[k][0] === pairs[k][1]) expect(other[k]).toBe(0);
      let best = 0;
      for (let mask = 0; mask < 1 << pairs.length; mask++) {
        const set = pairs.filter((_, k) => mask & (1 << k));
        if (set.length > best && free(set)) best = set.length;
      }
      expect(kept.length).toBe(best);
    }
  });

  it('`split`: the feet on the band’s middle, the nested set above it, the crossing pairs below, all in the band', async () => {
    const s = await settle(arcsSpec({ pairs: [...Stems, ...Crossing], crossing: 'split' }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const sec = f.layout.sections[0];
    const lay = arcsLay(sec);
    const top = sec.placement.strips.bandY[1];
    const h = sec.placement.strips.bandH[1];
    let below = 0;
    for (let i = 0; i < lay.n; i++) {
      expect(lay.curve[8 * i + 1]).toBeCloseTo(top + h / 2);
      const apex = at(lay.curve, 8 * i, 0.5)[1];
      expect(apex).toBeGreaterThanOrEqual(top);
      expect(apex).toBeLessThanOrEqual(top + h);
      if (apex > top + h / 2) below++;
    }
    // the tertiary pairs that cross the stems: at least one each crossing group, never a stem
    expect(below).toBeGreaterThan(0);
    expect(below).toBeLessThanOrEqual(Crossing.length);
  });

  it('past `pairs.maxArcs`: one arc a bin, the rest a density strip, said once', async () => {
    const pairs = randomPairs(400, 150, 5);
    const s = await settle(arcsSpec({ n: 400, pairs }), [], { pairs: { maxArcs: 20 } });
    const f = await s.renderer.frame({ size: Size, settled: true });
    const lay = arcsLay(f.layout.sections[0]);
    expect(lay.n).toBeLessThanOrEqual(20);
    expect(lay.densityMax).toBeGreaterThan(0);
    expect(f.layout.issues.map((i) => i.code)).toEqual(['seq.render.pairs.limit']);
  });

  it('on a ring: chords through the hole, ends on the band’s outer edge', async () => {
    const s = await settle(arcsSpec({ pairs: Stems, layout: 'rings', side: 'below' }));
    const f = await s.renderer.frame({ size: { width: 600, height: 600 }, settled: true });
    const sec = f.layout.sections[0];
    const lay = arcsLay(sec);
    expect(lay.n).toBe(Stems.length);
    const w = sec.placement.warpOf!(0);
    const r = w.r0 - sec.placement.strips.bandY[1];
    for (let i = 0; i < lay.n; i++) {
      const o = 8 * i;
      expect(Math.hypot(lay.curve[o] - w.cx, lay.curve[o + 1] - w.cy)).toBeCloseTo(r, 4);
      // the apex nearer the centre than the ends
      const [x, y] = at(lay.curve, o, 0.5);
      expect(Math.hypot(x - w.cx, y - w.cy)).toBeLessThan(r);
      const h = hitView(f.layout, x, y);
      expect(h?.element?.owner).toBe('/data/p');
    }
  });

  it('between blocks: an arc from one chain to the other; a chain not shown leaves a caret at the end shown', async () => {
    const b = createSeqSpec();
    const x = b.data.sequence('x', { alphabet: 'protein', letters: 'MKLVAEGHIKLMNPQ' });
    const y = b.data.sequence('y', { alphabet: 'protein', letters: 'MEEPQSDPSV' });
    const z = b.data.sequence('z', { alphabet: 'protein', letters: 'MAAA' });
    const p = b.data.pairs('c', {
      a: x,
      b: y,
      rows: pairRows([
        [3, 5],
        [10, 2],
      ]),
    });
    const q = b.data.pairs('d', { a: x, b: z, rows: pairRows([[4, 1]]) });
    b.views.view('v', (v) => {
      const sec = v.section('s', {});
      sec.block('x', { of: x });
      sec.block('y', { of: y });
      sec.track('arcs', {}, (t) => {
        t.layer('c', { type: 'pairs', data: p });
        t.layer('d', { type: 'pairs', data: q });
      });
    });
    const st = await settle(b.build());
    const f = await st.renderer.frame({ size: Size, settled: true });
    const sec = f.layout.sections[0];
    const c = sec.sectionLays[0].lay as PairsLay;
    const d = sec.sectionLays[1].lay as PairsLay;
    expect(c.n).toBe(2);
    // the first end on block x, the second on block y
    expect(c.curve[0]).toBeLessThan(sec.placement.strips.axes[1].x0);
    expect(c.curve[6]).toBeGreaterThan(sec.placement.strips.axes[1].x0);
    expect(d.n).toBe(0);
    expect(d.carets.length).toBe(4);
  });
});

describe('an RNA’s letters legend', () => {
  it('lists the nucleotides drawn: U, never T; no gap where there is none', async () => {
    const st = await settle(arcsSpec({ n: 24, pairs: [[0, 20]] }));
    const f = await st.frame({ size: Size, settled: true });
    const legend = viewLegends(f.layout).legends.find((l) => l.layers.some((k) => k.endsWith('/letters/letters')))!;
    const labels = legend.entries.kind === 'swatches' ? legend.entries.items.map((i) => i.label) : [];
    expect(labels.sort()).toEqual(['A', 'C', 'G', 'U']);
  });
});

describe('the pairs layer: dots on a grid', () => {
  it('a dot a pair at (a, b), mirrored at (b, a) too; hit as the pair; kept by the triangle', async () => {
    const pairs: [number, number][] = [
      [2, 15],
      [4, 12],
      [6, 9],
    ];
    for (const mirror of [false, true]) {
      const s = await settle(gridSpec({ n: 20, pairs, mirror }));
      const f = await s.renderer.frame({ size: Size, settled: true });
      const lay = gridLay(f.layout, 1, 1) as unknown as PairsLay;
      expect(lay.nd).toBe(mirror ? 6 : 3);
      const g = gridLay(f.layout);
      const a = g.axis!;
      const ay = g.axisY!;
      a.cell(4);
      const x = (a.sx0 + a.sx1) / 2;
      ay.cell(12);
      expect(hitView(f.layout, x, (ay.sx0 + ay.sx1) / 2)?.element?.key).toBe('4-12');
      ay.cell(4);
      a.cell(12);
      const back = hitView(f.layout, (a.sx0 + a.sx1) / 2, (ay.sx0 + ay.sx1) / 2);
      if (mirror) expect(back?.element?.key).toBe('4-12');
      else expect(back?.element).toBeUndefined();
    }
    const s = await settle(gridSpec({ n: 20, pairs, mirror: true, triangle: 'lower' }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    expect((gridLay(f.layout, 1, 1) as unknown as PairsLay).nd).toBe(3);
  });

  it('a self-pair mirrored on its own grid is one dot, the cell whole', async () => {
    const s = await settle(gridSpec({ n: 20, pairs: [[5, 5]], mirror: true }));
    const f = await s.renderer.frame({ size: Size, settled: true });
    const lay = gridLay(f.layout, 1, 1) as unknown as PairsLay;
    const a = gridLay(f.layout).axis!;
    expect(lay.nd).toBe(1);
    a.cell(5);
    expect(lay.dots[2]).toBeCloseTo(Math.max(Math.abs(a.sx1 - a.sx0), 2 / a.dpr));
  });

  it('a stretch runs with its partner when both ends are on one strand, against it otherwise', async () => {
    for (const strand of ['+', '-']) {
      const rows = { a_start: [2], a_end: [6], b_start: [10], b_end: [14], a_strand: ['+'], b_strand: [strand] };
      const s = await settle(gridSpec({ n: 20, rows }));
      const f = await s.renderer.frame({ size: Size, settled: true });
      const lay = gridLay(f.layout, 1, 1) as unknown as PairsLay;
      const ay = gridLay(f.layout).axisY!;
      expect(lay.ns).toBe(1);
      // from x's first: y's first (parallel) or y's last (antiparallel)
      ay.cell(strand === '+' ? 10 : 13);
      expect(lay.segs[1]).toBeCloseTo((ay.sx0 + ay.sx1) / 2);
    }
  });
});
