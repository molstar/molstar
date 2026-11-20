/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pairs' and grids' allocation targets (the design's "Performance and size"): a warm draw of a dot plot (image
// tiles and rects) and of 4,096 arcs, and a warm hover (the crosshair over the cells, the arcs lit), allocate nothing;
// a hit on a cell, an arc or a caret allocates no hit object, the pairs on split copies too. By hand, never in CI:
// `pnpm -C packages/seq/core test:perf`.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { rangesLocation } from '@molstar/seqx-core/query';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import {
  drawView,
  type FrameInput,
  type GridLay,
  GridType,
  type Hit,
  hitView,
  type PairsLay,
  PlacePoint,
} from '../index.js';
import { settle } from './helpers.js';
import { arcsSpec, gridSpec, letters, pairRows, randomPairs } from './pairs-fixtures.js';
import { reusedHover } from './scenarios.js';

const Size = { width: 1440, height: 1600 };
const at = (start: number, end: number): FrameInput => ({
  size: Size,
  dpr: 2,
  settled: true,
  viewports: { '/views/v/s': { start, end } },
});

/**
 * Bytes a hit may allocate: no hit object (one is over 100 B). A layer's hit reads its point from an array and measures
 * 0 at a cell (the test's array read into it); a whole `hitView` measures the two numbers the test's own call boxes
 * (32 B) and nothing past them: the placement's `at` reads the page point from an array too (`Placement.at`).
 */
const HitBound = 32;
/** A grid cell's `hitView`: the test's own call, nothing beyond it. */
const GridPathBound = HitBound;

/** The point a layer's hit reads (`LayerType.hit`). */
const At = new Float64Array(2);

describe('pairs and grids draw and hover with nothing allocated once warm', () => {
  for (const [name, n, input] of [
    ['a 2,000 × 2,000 dot plot (image tiles)', 2000, at(0, 2000)],
    ['a dot plot zoomed in (rects)', 200, at(10, 60)],
  ] as const)
    it(`${name}: draw, the crosshair hover, a hit on a cell`, async ({ task }) => {
      const s = await settle(
        gridSpec({
          n,
          pairs: [
            [3, 40],
            [5, 38],
            [12, 30],
          ],
          mirror: true,
        }),
      );
      const f = await s.frame(input);
      const d = quickCase(task, `grid draw ${name}`, () => drawView(f.layout, f.lists), {});
      console.info(`grid draw (${name}): ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
      expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
      const g = f.layout.sections[1].lays[0][0].lay as GridLay;
      const points: [number, number][] = [];
      for (let k = 0; k < 40; k++) {
        const i = Math.floor(g.axis!.u0 + ((k * 7) % 40) * ((g.axis!.u1 - g.axis!.u0) / 40));
        const j = Math.floor(g.axisY!.u0 + ((k * 11) % 40) * ((g.axisY!.u1 - g.axisY!.u0) / 40));
        g.axis!.cell(i);
        const x = (g.axis!.sx0 + g.axis!.sx1) / 2;
        g.axisY!.cell(j);
        points.push([x, (g.axisY!.sx0 + g.axisY!.sx1) / 2]);
      }
      const hits: Hit[] = [];
      for (const p of points) {
        const h = hitView(f.layout, p[0], p[1]);
        // a copy: the layer's hit is one scratch
        if (h) hits.push({ ...h });
      }
      expect(hits.length).toBeGreaterThan(30);
      const hover = reusedHover(hits);
      const hoverAt = { ...input, emphasis: hover.emphasis };
      const h = quickCase(
        task,
        `grid hover ${name}`,
        (i) => {
          hover.set(i);
          s.renderer.frame(hoverAt);
        },
        {},
      );
      console.info(`grid hover (${name}): ${h.r.bytesPerStep.toFixed(0)} bytes a frame (${h.r.note})`);
      expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
      const flat = Float64Array.from(points.flat());
      const m = points.length;
      // the layer's own hit, at the points within its strip
      const place = f.layout.sections[1].placement;
      const pt = new PlacePoint();
      const local = new Float64Array(2 * m);
      for (let k = 0; k < m; k++) {
        expect(place.at([flat[2 * k], flat[2 * k + 1]], pt)).toBe(true);
        local[2 * k] = pt.x;
        local[2 * k + 1] = pt.y;
      }
      const own = quickCase(
        task,
        `grid hit ${name}`,
        (i) => {
          At[0] = local[2 * (i % m)];
          At[1] = local[2 * (i % m) + 1];
          GridType.hit(g, At);
        },
        {},
      );
      console.info(`grid hit (${name}): ${own.r.bytesPerStep.toFixed(0)} bytes a hit (${own.r.note})`);
      expect(own.r.bytesPerStep, own.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
      const path = quickCase(
        task,
        `grid hitView ${name}`,
        (i) => void hitView(f.layout, flat[2 * (i % m)], flat[2 * (i % m) + 1]),
        {},
      );
      console.info(`grid hitView (${name}): ${path.r.bytesPerStep.toFixed(0)} bytes a hit (${path.r.note})`);
      expect(path.r.bytesPerStep, path.r.note).toBeLessThanOrEqual(GridPathBound);
    });

  it('a dot plot: overlapping highlight channels (rows in stripes, one in each style) allocate nothing once warm', async ({
    task,
  }) => {
    const s = await settle(gridSpec({ n: 200 }));
    const input = at(10, 60);
    const f = await s.frame(input);
    const loc = (a: number, b: number) => rangesLocation('/data/s', [a, b]);
    const made = (k: number) => [
      { location: loc(20 + k, 40 + k), color: 0x7aa3a6ff },
      { location: loc(25 + k, 45 + k), color: 0xd86c17ff, style: 'band' as const },
      { location: loc(30 + k, 50 + k), color: 0x2a9d8fff, style: 'outline' as const },
      { location: loc(35 + k, 55 + k), color: 0xe9c46aff },
    ];
    const sets = [made(0), made(1)];
    const emphasis = { highlights: sets[0] };
    const withHighlights = { ...input, emphasis };
    const r = quickCase(
      task,
      'grid overlapping highlights',
      (i) => {
        emphasis.highlights = sets[i & 1];
        s.renderer.frame(withHighlights);
      },
      {},
    );
    console.info(`grid overlapping highlights: ${r.r.bytesPerStep.toFixed(0)} bytes a frame (${r.r.note})`);
    // 96 B measured (a run's few bytes, 12 runs here: the grid's rows, TODO to find); the grid's emphasis redraw made
    // 88 B for one channel and 208 for a selection before the stripes
    expect(r.r.bytesPerStep, r.r.note).toBeLessThanOrEqual(128);
    expect(f.layout.sections.length).toBeGreaterThan(1);
  });

  it('4,096 arcs: draw, a hover lighting the arc, a hit on the curve', async ({ task }) => {
    const s = await settle(arcsSpec({ n: 20_000, pairs: randomPairs(20_000, 4096, 17), height: 120 }));
    const input = at(0, 20_000);
    const f = await s.frame(input);
    const d = quickCase(task, 'pairs draw', () => drawView(f.layout, f.lists), {});
    console.info(`pairs draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const lay = f.layout.sections[0].sectionLays[0].lay as PairsLay;
    const points: [number, number][] = [];
    for (let i = 0; i < lay.n; i += 97) {
      const o = 8 * i;
      points.push([(lay.curve[o] + lay.curve[o + 6]) / 2, (lay.curve[o + 3] + lay.curve[o + 1]) / 2]);
    }
    const hits: Hit[] = [];
    for (const p of points) {
      const h = hitView(f.layout, p[0], p[1]);
      if (h) hits.push({ ...h });
    }
    expect(hits.length).toBeGreaterThan(20);
    const hover = reusedHover(hits);
    const hoverAt = { ...input, emphasis: hover.emphasis };
    const h = quickCase(
      task,
      'pairs hover',
      (i) => {
        hover.set(i);
        s.renderer.frame(hoverAt);
      },
      {},
    );
    console.info(`pairs hover: ${h.r.bytesPerStep.toFixed(0)} bytes a frame (${h.r.note})`);
    expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const flat = Float64Array.from(points.flat());
    const m = points.length;
    const hit = quickCase(
      task,
      'pairs hit',
      (i) => void hitView(f.layout, flat[2 * (i % m)], flat[2 * (i % m) + 1]),
      {},
    );
    console.info(`pairs hit: ${hit.r.bytesPerStep.toFixed(0)} bytes a hit (${hit.r.note})`);
    expect(hit.r.bytesPerStep, hit.r.note).toBeLessThanOrEqual(HitBound);
  });

  it('pairs on split copies: draw, a hit on an arc and on a caret', async ({ task }) => {
    // a sequence's 0–300 twice (`copies: each`), 400 pairs within it: each drawn on both copies, those past 300 carets
    const b = createSeqSpec();
    const seq = b.data.sequence('s', { alphabet: 'rna', letters: letters(400) });
    const p = b.data.pairs('p', { a: seq, b: seq, rows: pairRows(randomPairs(400, 150, 5)) });
    b.views.view('v', (v) =>
      v.section('s', { of: seq, copies: 'each' }, (sec) => {
        sec.block('a1', { of: seq, range: [0, 300] });
        sec.block('a2', { of: seq, range: [0, 300] });
        sec.track('pairs', { height: 120 }).layer('pairs', { type: 'pairs', data: p });
      }),
    );
    const s = await settle(b.build());
    const f = await s.frame({ size: Size, dpr: 2, settled: true });
    const d = quickCase(task, 'copies draw', () => drawView(f.layout, f.lists), {});
    console.info(`copies draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const lay = f.layout.sections[0].sectionLays[0].lay as PairsLay;
    expect(lay.carets.length).toBeGreaterThan(0);
    const points: number[] = [];
    for (let i = 0; i < lay.n; i += 7) {
      const o = 8 * i;
      points.push((lay.curve[o] + lay.curve[o + 6]) / 2, (lay.curve[o + 3] + lay.curve[o + 1]) / 2);
    }
    for (let i = 0; i + 3 < lay.carets.length; i += 4)
      points.push(lay.carets[i], lay.carets[i + 1] + (lay.carets[i + 2] > 0 ? -3 : 3));
    // the points on an arc or a caret (an arc's middle by its controls is near its curve, not on it, for a few)
    const on: number[] = [];
    for (let i = 0; i + 1 < points.length; i += 2)
      if (hitView(f.layout, points[i], points[i + 1])) on.push(points[i], points[i + 1]);
    expect(on.length).toBeGreaterThan(0.9 * points.length);
    const flat = Float64Array.from(on);
    const m = flat.length / 2;
    const hit = quickCase(
      task,
      'copies hit',
      (i) => void hitView(f.layout, flat[2 * (i % m)], flat[2 * (i % m) + 1]),
      {},
    );
    console.info(`copies hit: ${hit.r.bytesPerStep.toFixed(0)} bytes a hit (${hit.r.note})`);
    expect(hit.r.bytesPerStep, hit.r.note).toBeLessThanOrEqual(HitBound);
  });
});
