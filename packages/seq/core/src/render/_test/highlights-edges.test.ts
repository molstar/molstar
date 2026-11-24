/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Overlapping highlight channels, the edges (_spec/README.md "Render", Emphasis): nested partial overlaps (counts and
// slots change mid-range), adjacency, a reversed axis, a cyclic origin, a break (one bar a channel's hidden range), a
// channel on another system, an outline whole over another channel's boundary, stripe edges on device pixels, and the
// grid's rows striped as the tracks' bands.

import { rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import {
  DrawLists,
  type EmphasisHighlight,
  drawView,
  type GridLay,
  recordDrawLists,
  type ViewLayout,
} from '../index.js';
import { renderSvg } from '../render-svg.js';
import { p53GeneSpec, p53TracksSpec } from './fixtures.js';
import { settle } from './helpers.js';
import { gridSpec } from './pairs-fixtures.js';
import { ringSpec } from './rings-fixtures.js';

const Protein = '/data/p53';
const Gene = '/data/tp53';
const Red = 0xff0000ff;
const Green = 0x00ff00ff;
const Blue = 0x0000ffff;

const on = (system: string, ranges: number[], color: number, style?: EmphasisHighlight['style']) =>
  ({ location: rangesLocation(system, ranges), color, ...(style ? { style } : {}) }) satisfies EmphasisHighlight;

interface R {
  readonly x: number;
  readonly y: number;
  readonly w: number;
  readonly h: number;
  readonly paint: string;
}

const rectsOf = (lists: Parameters<typeof recordDrawLists>[0]): R[] =>
  recordDrawLists(lists, ['highlight'])
    .split('\n')
    .filter((l) => l.startsWith('rect '))
    .map((l) => {
      const [, x, y, w, h, ...paint] = l.split(' ');
      return { x: Number(x), y: Number(y), w: Number(w), h: Number(h), paint: paint.join(' ') };
    });

const byColor = (rects: R[], hex: string) => rects.filter((r) => r.paint.startsWith(hex)).sort((a, b) => a.x - b.x);

const Tracks = { size: { width: 800, height: 240 }, viewport: { start: 0, end: 30 } };

async function tracks(highlights: EmphasisHighlight[], config = {}) {
  const r = await renderSvg(p53TracksSpec({ heights: 'fixed' }), { ...Tracks, config, emphasis: { highlights } });
  return rectsOf(r.lists);
}

describe('nested partial overlaps, adjacency, another system', () => {
  it('A [0,20), B [5,10), C [8,15): the slots and the count change in the middle of A’s range', async () => {
    const rects = await tracks([on(Protein, [0, 20], Red), on(Protein, [5, 10], Green), on(Protein, [8, 15], Blue)]);
    const top = rects.reduce((m, r) => Math.min(m, r.y), Infinity);
    const H = Math.max(...rects.map((r) => r.h));
    const slots = (hex: string) => byColor(rects, hex).map((r) => [(r.y - top) / H, r.h / H]);
    // [0,5) alone; [5,8) A, B; [8,10) A, B, C; [10,15) A, C; [15,20) alone
    expect(slots('#ff0000')).toEqual([
      [0, 1],
      [0, 1 / 2],
      [0, 1 / 3],
      [0, 1 / 2],
      [0, 1],
    ]);
    expect(slots('#00ff00')).toEqual([
      [1 / 2, 1 / 2],
      [1 / 3, 1 / 3],
    ]);
    expect(slots('#0000ff')).toEqual([
      [2 / 3, 1 / 3],
      [1 / 2, 1 / 2],
    ]);
  });

  it('adjacency: A ends where B starts, neither is split', async () => {
    const rects = await tracks([on(Protein, [4, 8], Red), on(Protein, [8, 12], Green)]);
    expect(rects.length).toBe(2);
    expect(rects[0].h).toBe(rects[1].h);
    expect(rects[0].y).toBe(rects[1].y);
    expect(rects[0].x + rects[0].w).toBe(rects[1].x);
  });

  it('a channel with its part on another system is not counted on this axis', async () => {
    const rects = await tracks([on(Protein, [4, 10], Red), on('/data/none', [4, 10], Green)]);
    expect(rects.length).toBe(1);
    expect(rects[0].h).toBe(96);
  });
});

describe('reversed axis, cyclic origin, a break', () => {
  it('a reversed axis: the stripes are the same, the runs in display order', async () => {
    const s = await settle(p53GeneSpec());
    const input = { size: { width: 800, height: 300 }, viewport: { start: 170, end: 180 } };
    await s.frame(input);
    const f = s.renderer.frame({
      ...input,
      emphasis: { highlights: [on(Gene, [826, 836], Red), on(Gene, [830, 840], Green)] },
    });
    const gene = f.layout.sections[1];
    const rects = rectsOf(f.lists).filter((r) => r.y >= gene.y - 1 && r.y < gene.y + gene.height);
    const reds = byColor(rects, '#ff0000');
    const greens = byColor(rects, '#00ff00');
    expect(reds.length).toBe(2);
    expect(greens.length).toBe(2);
    // red alone, red over green (a half), green under it, green alone: whole and half bands, either way round
    expect(new Set([...reds, ...greens].map((r) => r.h)).size).toBe(2);
    const half = Math.min(reds[0].h, reds[1].h);
    expect(reds.filter((r) => r.h === half).length).toBe(1);
    expect(greens.filter((r) => r.h === half).length).toBe(1);
    const a = gene.axis!;
    expect(a.x(831)).toBeLessThan(a.x(830));
  });

  it('a cyclic origin: a channel across it is two runs, the other’s stripes where they overlap', async () => {
    const hs = [on('/data/construct', [0, 40, 560, 600], Red), on('/data/construct', [20, 80], Green)];
    const r = await renderSvg(ringSpec(), { size: { width: 600, height: 600 }, emphasis: { highlights: hs } });
    const paths = recordDrawLists(r.lists, ['highlight'])
      .split('\n')
      .filter((l) => l.startsWith('path '))
      .map((l) => l.slice(l.lastIndexOf(' ') + 1));
    // red: [0,20) alone, [20,40) striped, [560,600) alone; green: [20,40) striped, [40,80) alone
    expect(paths.filter((c) => c === '#ff00004d').length).toBe(3);
    expect(paths.filter((c) => c === '#00ff004d').length).toBe(2);
  });

  it('a range wholly in a hidden break lights it once a channel, however its runs split', async () => {
    const s = await settle(p53GeneSpec({ coding: true }));
    const input = { size: { width: 800, height: 300 }, viewport: { start: 0, end: 393 } };
    await s.frame(input);
    const hidden = async (hs: EmphasisHighlight[]) => {
      const f = s.renderer.frame({ ...input, emphasis: { highlights: hs } });
      const gene = f.layout.sections[1];
      return rectsOf(f.lists).filter((r) => r.y >= gene.y - 1 && r.y < gene.y + gene.height);
    };
    // the intron [400, 480) of the gene is hidden: A and B lie in it and overlap, so A’s and B’s runs are split
    const rects = await hidden([on(Gene, [410, 440], Red), on(Gene, [420, 450], Green)]);
    expect(byColor(rects, '#ff0000').length).toBe(1);
    expect(byColor(rects, '#00ff00').length).toBe(1);
  });
});

describe('an outline is whole; stripe edges are on device pixels', () => {
  it('A outlined [0,100) with B [40,60): A is one outline, B keeps the whole band', async () => {
    const rects = await renderSvg(p53TracksSpec({ heights: 'fixed' }), {
      size: { width: 1200, height: 240 },
      viewport: { start: 0, end: 120 },
      emphasis: { highlights: [on(Protein, [0, 100], Red, 'outline'), on(Protein, [40, 60], Green)] },
    }).then((r) => rectsOf(r.lists));
    const red = rects.filter((r) => r.paint.includes('#ff0000'));
    const green = rects.filter((r) => r.paint.includes('#00ff00'));
    expect(red.length).toBe(1);
    expect(green.length).toBe(1);
    // the outline is inside its box by half its width; the tint keeps the whole band
    expect(green[0].h).toBe(96);
    expect(red[0].h).toBeGreaterThan(90);
  });

  it.each([1, 2])('five stripes of a band not divisible by five meet edge to edge at dpr %i', async (dpr) => {
    const hs = [Red, Green, Blue, 0xffcc00ff, 0xff66ccff].map((c) => on(Protein, [4, 10], c));
    const r = await renderSvg(p53TracksSpec({ heights: 'fixed' }), {
      ...Tracks,
      dpr,
      config: { emphasis: { highlightStripes: 5 } },
      emphasis: { highlights: hs },
    });
    const rects = rectsOf(r.lists);
    expect(rects.length).toBe(5);
    for (let i = 0; i + 1 < 5; i++) expect(rects[i].y + rects[i].h).toBe(rects[i + 1].y);
    // each edge on a device pixel
    for (const q of rects) {
      expect(Math.round(q.y * dpr)).toBeCloseTo(q.y * dpr, 6);
      expect(Math.round((q.y + q.h) * dpr)).toBeCloseTo((q.y + q.h) * dpr, 6);
    }
    expect(rects[4].y + rects[4].h - rects[0].y).toBe(96);
  });
});

describe('the grid: rows striped across the plot as the tracks’ bands', () => {
  it('overlapping channels split a row into stripes across the plot; an outline stays whole; a band is a bar', async () => {
    const s = await settle(gridSpec({ n: 20 }));
    const f = await s.renderer.frame({ size: { width: 720, height: 900 }, settled: true });
    const sec = f.layout.sections[1];
    const lay = sec.lays[0][0].lay as GridLay;
    const W = lay.axis!.width;
    const draw = (hs: EmphasisHighlight[]) => {
      const lists = new DrawLists();
      drawView(f.layout as ViewLayout, lists, { highlights: hs }, ['highlight']);
      const r = lists.passes.highlight.rects;
      const out: { x: number; y: number; w: number; h: number; color: number }[] = [];
      for (let k = 0; k < r.n; k++) out.push({ x: r.x[k], y: r.y[k], w: r.w[k], h: r.h[k], color: r.fill[k] });
      // the rows: across the whole plot or its stripes
      return out;
    };
    const S = '/data/s';
    const rowsOnly = (q: { w: number }) => q.w >= W / 5 - 1;
    const rows = draw([on(S, [4, 10], Red), on(S, [6, 12], Green)]).filter(rowsOnly);
    const full = rows.filter((q) => Math.abs(q.w - W) < 1);
    const half = rows.filter((q) => Math.abs(q.w - W / 2) < 1);
    expect(full.length).toBe(2);
    expect(half.length).toBe(2);
    // in channel order across the plot: red’s stripe first, green’s after it
    const x0 = lay.axis!.x0;
    expect(half.find((q) => q.color >>> 8 === 0xff0000)!.x).toBeCloseTo(x0, 3);
    expect(half.find((q) => q.color >>> 8 === 0x00ff00)!.x).toBeCloseTo(x0 + W / 2, 3);
    // an outline takes no stripe: the green alone keeps the whole row
    const outlined = draw([on(S, [4, 10], Red, 'outline'), on(S, [6, 12], Green)]).filter(rowsOnly);
    expect(outlined.filter((q) => Math.abs(q.w - W) < 1).length).toBe(2);
    expect(outlined.filter((q) => Math.abs(q.w - W / 2) < 1).length).toBe(0);
    // a band is a bar of its width
    const band = draw([on(S, [4, 10], Red, 'band')]);
    expect(band.some((q) => Math.abs(q.w - 3) < 1e-6)).toBe(true);
  });
});
