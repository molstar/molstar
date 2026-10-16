/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// XR1a (seqspec:viewer#extension-points): every renderer a registry entry, one module per kernel type under
// `renderers/` on the renderer API, registered by `DefaultRegistry`; the placement helpers in each geometry (block
// coordinate system runs placed through the axis's converters; a filter and a reversed piece; a wrapped line; rings).
// A renderer module is imported first here: the API module evaluates before the root entry with no cycle.

import { readdirSync, readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { SectionAxisView } from '../src/axis.js';
import { AxisGeometry, band, CHARS, DrawView, mark, placeRun, RectBatch, type ShownIndex } from '../src/draw/api.js';
import { createSeqViewer, KernelRenderers } from '../src/index.js';
import { DefaultSeqViewerSpec } from '../src/spec.js';
import { RangesRenderer } from '../src/renderers/ranges.js';

import { layoutRuns, lineTable, prepareRun } from '../src/runs.js';
import { axisOf, rangedDoc, resolved, setPiece } from './axis-fixtures.js';

// Vite resolves a static import's `.js` to the `.ts`; a computed one it cannot.
const LOADERS: Record<string, () => Promise<unknown>> = {
  letters: () => import('../src/renderers/letters.js'),
  ranges: () => import('../src/renderers/ranges.js'),
  points: () => import('../src/renderers/points.js'),
  emphasis: () => import('../src/renderers/emphasis.js'),
  variants: () => import('../src/renderers/variants.js'),
  secondary: () => import('../src/renderers/secondary.js'),
  values: () => import('../src/renderers/values.js'),
  color: () => import('../src/renderers/color.js'),
  ruler: () => import('../src/renderers/ruler.js'),
  pairs: () => import('../src/renderers/pairs.js'),
  rows: () => import('../src/renderers/rows.js'),
  logo: () => import('../src/renderers/logo.js'),
  matrix: () => import('../src/renderers/matrix.js'),
  map: () => import('../src/renderers/map.js'),
};
const MODULES = ['letters', 'ranges', 'points', 'emphasis', 'variants', 'secondary', 'values', 'color', 'ruler'];
const MORE = ['pairs', 'rows', 'logo', 'matrix', 'map'];

describe('renderer modules', () => {
  it('one module per kernel type and seqv:map, each its own /renderers/* subpath, all in DefaultRegistry', async () => {
    const ctx = createSeqViewer(DefaultSeqViewerSpec());
    for (const name of [...MODULES, ...MORE]) {
      const m = (await LOADERS[name]()) as Record<string, unknown>;
      const entries = Object.values(m);
      expect(entries).toHaveLength(1);
      const e = entries[0] as { name: string };
      expect(e.name).toBe(name === 'map' ? 'seqv:map' : name);
      expect(ctx.registry.get('renderers', e.name)).toBe(e);
      expect(KernelRenderers).toContain(e);
    }
    expect(ctx.registry.get('renderers', 'ranges')).toBe(RangesRenderer);
    ctx.dispose();
  });

  it('one character table (`CHARS`, a string per 7-bit code) in the renderer API, used by every renderer that draws letters', () => {
    expect(CHARS).toHaveLength(128);
    expect(CHARS[65]).toBe('A');
    const table = /String\.fromCharCode\(i\)/;
    const src = (f: string) => readFileSync(new URL(`../src/${f}`, import.meta.url), 'utf8');
    const own = readdirSync(new URL('../src/renderers/', import.meta.url))
      .map((f) => `renderers/${f}`)
      .concat('draw/points.ts')
      .filter((f) => table.test(src(f)));
    expect(own).toEqual([]);
    expect(table.test(src('draw/backend.ts'))).toBe(true);
    for (const f of [
      'renderers/letters.ts',
      'renderers/emphasis.ts',
      'renderers/logo.ts',
      'renderers/rows.ts',
      'draw/points.ts',
    ])
      expect(src(f), f).toMatch(/\bCHARS\[/);
  });
});

/** A draw view over a section's axis: data area from x 10, `cell` CSS px, `dpr`, shown [lo, hi) of block 0 in view. */
async function view(o: { cell?: number; dpr?: number; lo?: number; hi?: number; reversed?: boolean; v?: string } = {}) {
  const r = await resolved(rangedDoc(), o.v ?? 'f');
  if (o.reversed) setPiece(r.view, 0, 5, 15, -1, 20);
  const ax = axisOf(r);
  const s = new SectionAxisView();
  const n = ax.shown[0];
  s.sync(ax, { geometry: AxisGeometry.linear, dataLeft: 10, dataRight: 1000, scrolledTop: 0, height: 100 }, {
    visLo: Int32Array.of(o.lo ?? 0),
    visHi: Int32Array.of(o.hi ?? n),
    cell: o.cell ?? 10,
    lod: 0,
    dpr: o.dpr ?? 2,
    scrollY: 0,
    scrollX: 0,
    contentHeight: 0,
  } as never);
  return { ax, v: new DrawView().from(s) };
}
const out = new Float64Array(16);
const pieces = (n: number) => Array.from({ length: n }, (_, p) => [out[2 * p], out[2 * p + 1]]);
const k = (n: number) => n as ShownIndex;

describe('placement helpers', () => {
  it('linear: a run over a filter (block 7..10 and 12..14 shown) as its shown runs, device px', async () => {
    const { v } = await view();
    // View f shows unit [5, 15)'s 7, 8, 9, 12, 13: shown 0..4.
    const n = placeRun(v, 0, 7, 14, out);
    const px = (s: number) => Math.round(v.x(0, k(s)) * 2);
    // Shown runs that touch merge (the converters' rule): one piece across the break.
    expect(pieces(n)).toEqual([[px(0), px(5)]]);
    expect(placeRun(v, 0, 10, 12, out)).toBe(0);
    expect(pieces(placeRun(v, 0, 8, 13, out)).flat()[0]).toBe(px(1));
  });

  it("linear: clipped to the shown indices in view; another block's lane (stacked) places nothing", async () => {
    const { v } = await view({ lo: 1, hi: 3 });
    const px = (s: number) => Math.round(v.x(0, k(s)) * 2);
    expect(pieces(placeRun(v, 0, 5, 15, out))).toEqual([[px(1), px(3)]]);
    v.onlyBlock = 1;
    expect(placeRun(v, 0, 5, 15, out)).toBe(0);
  });

  it("a reversed piece: the run's x order flips, the edges stay ascending", async () => {
    const { v, ax } = await view({ reversed: true, v: 'v' });
    // Unit [5, 15) reversed: shown 0 is index 14.
    expect(ax.indexOf(0, k(0))).toBe(14);
    const n = placeRun(v, 0, 12, 14, out);
    expect(n).toBe(1);
    expect([out[0], out[1]]).toEqual([Math.round(v.x(0, k(1)) * 2), Math.round(v.x(0, k(3)) * 2)]);
    expect(out[0]).toBeLessThan(out[1]);
  });

  it("band: a rect per piece in device px; mark: a cell's center, false where not shown", async () => {
    const { v } = await view({ v: 'v' });
    const rects = new RectBatch();
    rects.begin(0x336699);
    band(v, 0, 6, 9, 20, 30, rects);
    expect(rects.n).toBe(1);
    expect([rects.x[0], rects.y[0], rects.w[0], rects.h[0]]).toEqual([
      Math.round(v.x(0, k(1)) * 2),
      40,
      Math.round(v.x(0, k(4)) * 2) - Math.round(v.x(0, k(1)) * 2),
      20,
    ]);
    expect(mark(v, 0, 6, 25, out)).toBe(true);
    expect([out[0], out[1]]).toEqual([(v.x(0, k(1)) + 5) * 2, 50]);
    expect(mark(v, 0, 2, 25, out)).toBe(false);
  });

  it("wrapped: pieces on the view's line only, the last cell's right edge across a group gap", async () => {
    const { ax, v } = await view({ v: 'v' });
    const run = prepareRun(0, ax, 10, 4);
    const lines = lineTable();
    // 10 positions in groups of 4, 60 px a line (gap 6): lines [0, 4), [4, 8), [8, 10).
    layoutRuns([run], 60, 6, lines);
    expect(Array.from(lines.start.subarray(0, lines.count))).toEqual([0, 4, 8]);
    v.geometry = AxisGeometry.wrapped;
    v.runs = [run];
    v.gap = 6;
    v.onlyBlock = 0;
    v.lineStart = 4;
    v.lineStop = 8;
    // Unit [5, 15): indices 7..12 are shown 2..7; the line holds shown 4..7.
    expect(pieces(placeRun(v, 0, 7, 13, out))).toEqual([[20, 100]]);
    v.lineStart = 0;
    v.lineStop = 4;
    expect(pieces(placeRun(v, 0, 7, 13, out))).toEqual([[60, 100]]);
  });
});
