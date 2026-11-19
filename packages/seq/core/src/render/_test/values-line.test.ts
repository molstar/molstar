/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Values as a line and as an area (_spec/README.md "Render", Layer types; EC-57): the line through the cells' centers,
// held flat over the half cells at its ends and at gaps, broken by missing values, its color changing halfway between
// two positions, the area under it; the same over a tile edge, at the device-column step, and on a reversed, collapsed
// axis.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { DrawLists, LightTokens, type ValuesLay, ValuesType } from '../index.js';
import { BigRuntime, conservation, GENE_EXONS, GENE_LENGTH, geneLetters, P53, seeded } from './fixtures.js';
import { settle } from './helpers.js';
import { layerLay, readPaths, readRects } from './small-helpers.js';

const MISSING = LightTokens.colors['seq-missing'];

function p53Values(encoding: 'bars' | 'line' | 'area') {
  const b = createSeqSpec({ title: 'line' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const c = b.data.values('conservation', { target: p53, values: conservation() });
  b.views.view('main', (v) => {
    v.section('protein', { of: p53 }).track('c').layer('c', { type: 'values', data: c, encoding });
  });
  return b.build();
}

/** A value's y in the layer's box, through its domain. */
function yOf(lay: ValuesLay, v: number): number {
  const [lo, hi] = lay.domain!;
  const y1 = lay.box.y + lay.box.height;
  return y1 - ((Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo)) * lay.box.height;
}

describe('values as a line', () => {
  it('is a path a run of values: half a cell held flat at each end, the centers between; missing breaks it', async () => {
    const s = await settle(p53Values('line'));
    const f = await s.frame({ size: { width: 800, height: 60 } });
    const lay = layerLay<ValuesLay>(f, 0);
    const paths = readPaths(f.lists).filter((p) => !p.closed);
    // two missing positions (100 and 101) split the line in two; their marks stay at the baseline
    expect(paths.length).toBe(2);
    const v = conservation();
    const axis = lay.axis;
    const runs = [
      [0, 100],
      [102, 393],
    ];
    paths.forEach((p, k) => {
      const [a, b] = runs[k];
      // one center a position, a point at each end of the run
      expect(p.xy.length / 2).toBe(b - a + 2);
      expect(p.stroke).toBe(lay.fill);
      expect(p.fill).toBe(0);
      expect(p.strokeWidth).toBe(1.5);
      axis.cell(a);
      const x0 = axis.snap(axis.sx0);
      axis.cell(b - 1);
      const x1 = axis.snap(axis.sx1);
      expect(p.xy[0]).toBeCloseTo(x0, 4);
      expect(p.xy[p.xy.length - 2]).toBeCloseTo(x1, 4);
      // held flat over the first half cell and the last
      expect(p.xy[1]).toBeCloseTo(yOf(lay, v[a]), 3);
      expect(p.xy[3]).toBeCloseTo(p.xy[1], 5);
      expect(p.xy[p.xy.length - 1]).toBeCloseTo(p.xy[p.xy.length - 3], 5);
      for (let i = 0; i < b - a; i++) {
        expect(p.xy[2 * (i + 1) + 1]).toBeCloseTo(yOf(lay, v[a + i]), 3);
        if (i) expect(p.xy[2 * (i + 1)]).toBeGreaterThan(p.xy[2 * i]);
      }
    });
    // a neutral mark at the baseline for each missing position: never a zero
    const marks = readRects(f.lists).filter((r) => r.fill === MISSING);
    expect(marks.length).toBe(2);
    for (const m of marks) expect(m.y + m.h).toBeCloseTo(lay.box.y + lay.box.height, 4);
  });

  it('a color that changes does so halfway between two positions, the line going on to the boundary', async () => {
    const s = await settle(p53Values('line'));
    const f = await s.frame({ size: { width: 800, height: 60 } });
    const lay = layerLay<ValuesLay>(f, 0);
    // positions 200 on drawn in another color
    const other = 0xc2410cff;
    const first = lay.ord.subarray(0, lay.n).findIndex((j) => lay.x0[j] >= lay.axis.x(200) - 1e-6);
    for (let k = first; k < lay.n; k++) lay.col[lay.ord[k]] = other;
    const lists = new DrawLists();
    lists.reset();
    ValuesType.draw(lay, 'data', lists.writer('data'));
    const paths = readPaths(lists).filter((p) => !p.closed);
    // the run before the gap, then the run after it split at position 200
    expect(paths.map((p) => p.stroke)).toEqual([lay.fill, lay.fill, other]);
    const before = paths[1];
    const after = paths[2];
    const v = conservation();
    // the first ends at the boundary of cell 199 and 200, the second starts there, at the same height: halfway between
    // the centers of the two cells
    const bx = lay.axis.snap(lay.axis.x(200));
    expect(before.xy[before.xy.length - 2]).toBeCloseTo(bx, 4);
    expect(after.xy[0]).toBeCloseTo(bx, 4);
    const mid = (yOf(lay, v[199]) + yOf(lay, v[200])) / 2;
    expect(before.xy[before.xy.length - 1]).toBeCloseTo(mid, 3);
    expect(after.xy[1]).toBeCloseTo(mid, 3);
  });

  it('is drawn with no bars: the layer’s marks are its paths and the missing marks', async () => {
    const bars = await (await settle(p53Values('bars'))).frame({ size: { width: 800, height: 60 } });
    const line = await (await settle(p53Values('line'))).frame({ size: { width: 800, height: 60 } });
    expect(bars.lists.passes.data.paths.n).toBe(0);
    expect(line.lists.passes.data.paths.n).toBe(2);
    expect(line.lists.passes.data.rects.n).toBeLessThan(bars.lists.passes.data.rects.n);
  });
});

describe('values as an area', () => {
  it('is a polygon under each run of the line, down to the baseline, at 0.35 opacity, behind the line', async () => {
    const s = await settle(p53Values('area'));
    const f = await s.frame({ size: { width: 800, height: 60 } });
    const lay = layerLay<ValuesLay>(f, 0);
    const paths = readPaths(f.lists);
    // for each of the two runs, the area then the line
    expect(paths.map((p) => p.closed)).toEqual([true, false, true, false]);
    const base = lay.box.y + lay.box.height;
    paths.forEach((p, i) => {
      if (i % 2) return;
      const line = paths[i + 1];
      expect(p.fill & 0xff).toBe(89);
      expect(p.stroke).toBe(0);
      expect(p.xy.length).toBe(line.xy.length + 4);
      // the line's points, then down to the baseline at its last x and back at its first
      expect(p.xy.slice(0, line.xy.length)).toEqual(line.xy);
      expect(p.xy[p.xy.length - 1]).toBeCloseTo(base, 5);
      expect(p.xy[p.xy.length - 3]).toBeCloseTo(base, 5);
      expect(p.xy[p.xy.length - 4]).toBeCloseTo(line.xy[line.xy.length - 2], 5);
      expect(p.xy[p.xy.length - 2]).toBeCloseTo(line.xy[0], 5);
    });
  });
});

describe('values as a line: across a tile edge and at the device-column step', () => {
  function big(encoding: 'line' | 'area') {
    const b = createSeqSpec({ title: 'big line' });
    const src = b.data.source('synthetic', { format: 'synthetic', uri: 'synthetic:chr' });
    const chr = b.data.sequence('chr', { alphabet: 'dna', source: src });
    const v = b.data.values('signal', { target: chr, source: src });
    b.views.view('main', (view) => {
      view.section('chr', { of: chr }).track('signal').layer('signal', { type: 'values', data: v, encoding });
    });
    return b.build();
  }

  it('over the edge between two tiles the line is one run', async () => {
    const s = await settle(big('line'), [BigRuntime]);
    // the tile edge at 1024
    const f = await s.frame({ size: { width: 1000, height: 60 }, viewport: { start: 1000, end: 1050 } });
    const lay = layerLay<ValuesLay>(f, 0);
    expect(lay.n).toBe(50);
    const paths = readPaths(f.lists);
    expect(paths.length).toBe(1);
    expect(paths[0].xy.length / 2).toBe(52);
  });

  it('at the device-column step it follows each column’s maximum, the column’s extent drawn under it', async () => {
    const s = await settle(big('line'), [BigRuntime]);
    const f = await s.frame({ size: { width: 800, height: 60 }, viewport: { start: 0, end: 8000 } });
    const lay = layerLay<ValuesLay>(f, 0);
    expect(['bars', 'density']).toContain(lay.step);
    const paths = readPaths(f.lists);
    expect(paths.length).toBeGreaterThan(0);
    // a point a column in view, at its maximum
    const points = paths.reduce((a, p) => a + p.xy.length / 2 - 2, 0);
    let valid = 0;
    for (let j = 0; j < lay.n; j++) if (!Number.isNaN(lay.v[j])) valid++;
    expect(points).toBe(valid);
    for (const p of paths)
      for (let i = 0; i < p.xy.length; i += 2)
        expect(p.xy[i + 1]).toBeLessThanOrEqual(lay.box.y + lay.box.height + 1e-6);
  });
});

describe('values as a line: a reversed, collapsed axis', () => {
  it('a run a piece, in screen order, each point at its position’s value', async () => {
    const b = createSeqSpec({ title: 'reversed line' });
    const g = b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters() });
    const r = seeded(3);
    const values = Float32Array.from({ length: GENE_LENGTH }, () => r());
    const signal = b.data.values('signal', { target: g, values });
    const exons = b.data.ranges('exons', { target: g, rows: GENE_EXONS });
    b.views.view('main', (v) => {
      v.section('gene', { of: g, direction: 'reversed', pieces: exons })
        .track('s')
        .layer('s', { type: 'values', data: signal, encoding: 'area' });
    });
    const s = await settle(b.build());
    const f = await s.frame({ size: { width: 2000, height: 60 } });
    const lay = layerLay<ValuesLay>(f, 0);
    const lines = readPaths(f.lists).filter((p) => !p.closed);
    // three exons: three runs, the last exon (highest positions) first on the screen
    expect(lines.length).toBe(3);
    const xs = lines.map((p) => p.xy[0]);
    expect([...xs].sort((a, b) => a - b)).toEqual(xs);
    lines.forEach((p, k) => {
      const exon = 2 - k;
      for (let i = 1; i + 1 < p.xy.length / 2; i++) {
        const pos = lay.axis.positionAt(p.xy[2 * i]);
        expect(pos).toBeGreaterThanOrEqual(GENE_EXONS.start[exon]);
        expect(pos).toBeLessThan(GENE_EXONS.end[exon]);
        expect(p.xy[2 * i + 1]).toBeCloseTo(yOf(lay, values[pos]), 3);
      }
    });
    // and each run has its area
    expect(readPaths(f.lists).filter((p) => p.closed).length).toBe(3);
  });
});
