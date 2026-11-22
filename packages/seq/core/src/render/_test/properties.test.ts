/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The draw lists' properties: bounded by the screen at every density step (a 3M-position system and p53 alike), hits
// where the marks are, pure (same inputs, same lists), keyed by what they read (the token table's version among it).

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import {
  colorHex,
  createTokens,
  DarkTokens,
  type DensityStep,
  densityStep,
  DefaultRenderConfig,
  hitView,
  LaneIndex,
  LightTokens,
  type LettersLay,
  type RangesLay,
  recordDrawLists,
  type SeqTokenName,
} from '../index.js';
import { renderSvg } from '../render-svg.js';
import { rangesLocation } from '@molstar/seqx-core/query';
import { BIG, BigRuntime, bigSpec, P53, p53TracksSpec } from './fixtures.js';
import { blocks, settle } from './helpers.js';

const DataWidth = 600;
/** p53's domains (the fixture's). */
const Domains = [
  [0, 61],
  [63, 94],
  [93, 292],
  [318, 360],
  [323, 355],
  [363, 393],
];
const { inset, labelColumn, labelGap } = DefaultRenderConfig;
// the label column measured: its least width (every label here is shorter)
const Width = DataWidth + inset + labelColumn.min + labelGap + inset;
const Size = { width: Width, height: 120 };
/** Lanes a track's box holds: letters 1, ranges 3 (36 px of 12), values 1. */
const Lanes = [1, 3, 1];

/** Pixels per position (CSS = device at DPR 1) that land on each step. */
const Steps: [DensityStep, number][] = [
  ['letters', 20],
  ['cells', 5],
  ['dots', 1.5],
  ['bars', 0.2],
  ['density', 600 / BIG],
];

describe('draw lists: bounded by the screen', () => {
  it('at every density step, at most a primitive a device column a lane, for 3M positions and for p53', async () => {
    const big = bigSpec();
    const p53 = p53TracksSpec();
    for (const [step, px] of Steps) {
      expect(densityStep(px, 1, DefaultRenderConfig)).toBe(step);
      const n = DataWidth / px;
      const mid = Math.floor(BIG / 2);
      const viewport = step === 'density' ? { start: 0, end: BIG } : { start: mid, end: mid + n };
      const b = await renderSvg(big, { size: Size, viewport, runtimes: [BigRuntime] });
      expect(b.settled).toBe(true);
      const counts = blocks(b.lists);
      expect(counts.length).toBe(3);
      counts.forEach((c, t) => {
        expect(c.rects + c.images + c.texts).toBeLessThanOrEqual(DataWidth * Lanes[t]);
        expect(c.letters).toBeLessThanOrEqual(DataWidth);
      });
      const [letters, items, signal] = counts;
      if (step === 'letters' || step === 'cells' || step === 'dots') {
        // a cell a position: one rect each, and a letter each when letters show
        expect(letters.rects).toBe(Math.ceil(n));
        expect(letters.letters).toBe(step === 'letters' ? Math.ceil(n) : 0);
        // a bar a position from 4 px a position (a profile below), and the baseline's hairline
        expect(signal.rects).toBeLessThanOrEqual(Math.ceil(n) + 1);
        expect(signal.rects).toBeGreaterThanOrEqual(1);
      } else {
        // one bar for the letters, a bar a device column for the values
        expect(letters.rects).toBe(1);
        expect(signal.rects).toBeLessThanOrEqual(DataWidth);
      }
      if (step === 'density') expect(items).toMatchObject({ rects: 0, images: 1 });

      // p53 shows the same positions where it has them: the same counts a device column
      if (n > P53.length) continue;
      const p = await renderSvg(p53, { size: Size, viewport: { start: 100, end: 100 + n } });
      const [pl, pr, pv] = blocks(p.lists);
      expect(pl).toEqual(letters);
      expect(pv.rects).toBeLessThanOrEqual(signal.rects);
      // ranges: the same bound, and a mark (with at most its label) a domain in view
      const inView = Domains.filter(([a, b]) => a < 100 + n && b > 100).length;
      expect(pr.rects).toBe(inView);
      expect(pr.rects + pr.images + pr.texts).toBeLessThanOrEqual(DataWidth * Lanes[1]);
    }
  });

  it('p53 whole is bounded too: a cell a position at its own step, the values a profile', async () => {
    const r = await renderSvg(p53TracksSpec(), { size: Size });
    const [letters, items, values] = blocks(r.lists);
    expect(letters.rects).toBe(P53.length);
    expect(items.rects).toBeLessThanOrEqual(6);
    // below 4 px a position the values are a profile: its paths and the baseline's hairline
    expect(values.rects).toBe(1);
  });
});

describe('draw lists: hits where the marks are', () => {
  it('every range mark answers its item, every cell its position, a gap no item', async () => {
    const s = await settle(p53TracksSpec());
    const f = await s.frame({ size: Size, viewport: { start: 40, end: 140 } });
    const track = (i: number) => f.layout.sections[0].tracks[i];
    const ranges = track(1).layers[0].lay as RangesLay;
    expect(ranges.m).toBeGreaterThan(1);
    for (let j = 0; j < ranges.m; j++) {
      const x = (ranges.mx0[j] + ranges.mx1[j]) / 2;
      const y = ranges.laneY + (ranges.mlane[j] + 0.5) * ranges.laneHeight;
      const hit = hitView(f.layout, x, y);
      const item = ranges.mitem[j];
      expect(hit?.element).toEqual({ owner: '/data/domains', key: ranges.keys[item] });
      expect(hit?.extent).toEqual(rangesLocation('/data/p53', [ranges.starts[item], ranges.ends[item]]));
      expect(hit?.system).toBe('/data/p53');
      expect(hit?.section).toBe('protein');
    }
    // lane 0 between TAD's end (61) and PRR's start (63): no item
    const gap = s.renderer.lists;
    expect(gap.size).toBeGreaterThan(0);
    const x = f.layout.sections[0].axis!.x(62) + 1;
    const between = hitView(f.layout, x, ranges.laneY + 0.5 * ranges.laneHeight);
    expect(between).toEqual({
      layer: '/views/main/protein/domains/ranges',
      section: 'protein',
      position: 62,
      system: '/data/p53',
    });

    const letters = track(0).layers[0].lay as LettersLay;
    for (let j = 0; j < letters.n; j++) {
      const hit = hitView(f.layout, (letters.x0[j] + letters.x1[j]) / 2, letters.box.y + 1);
      expect(hit).toEqual({
        layer: '/views/main/protein/letters/letters',
        section: 'protein',
        position: letters.pos[j],
        system: '/data/p53',
      });
    }
    // the label column and the page's margin hit nothing
    expect(hitView(f.layout, 4, letters.box.y + 1)).toBeUndefined();
    expect(hitView(f.layout, 200, 2)).toBeUndefined();
  });

  it('seeded random items in lanes: every drawn mark is hit at its centre by an item it draws, merged runs too', async () => {
    const rng = new PCG(5);
    let merged = 0;
    for (let round = 0; round < 12; round++) {
      const length = 200 + (rng.int() % 4000);
      const n = 1 + (rng.int() % 40);
      const start: number[] = [];
      const end: number[] = [];
      const id: string[] = [];
      for (let i = 0; i < n; i++) {
        const s = rng.int() % length;
        start.push(s);
        end.push(Math.min(length, s + 1 + (rng.int() % Math.max(1, length >> (2 + (rng.int() % 6))))));
        // a few ids repeated: pieces of one item
        id.push(`i${rng.int() % Math.max(1, n - (rng.int() % 4))}`);
      }
      const b = createSeqSpec({ title: 'random ranges' });
      const s = b.data.sequence('s', { alphabet: 'protein', letters: 'M'.repeat(length) });
      const r = b.data.ranges('r', {
        target: s,
        rows: { start: Int32Array.from(start), end: Int32Array.from(end), id: categoryColumn(id) },
      });
      b.views.view('v', (v) =>
        v.section('s', { of: s }).track('r', { height: 96 }).layer('ranges', { type: 'ranges', data: r }),
      );
      const settled = await settle(b.build());
      const a = rng.int() % length;
      const viewports = [
        { start: 0, end: length },
        { start: a, end: Math.min(length, a + 20 + (rng.int() % 200)) },
      ];
      for (const viewport of viewports) {
        const f = await settled.frame({ size: { width: 728, height: 140 }, viewport });
        const lay = f.layout.sections[0].tracks[0].layers[0].lay as RangesLay;
        if (lay.mode === 'density') continue;
        for (let j = 0; j < lay.m; j++) {
          const drawn = new Set<string | number>();
          if (lay.mitem[j] >= 0) drawn.add(lay.keys[lay.mitem[j]]);
          else {
            merged++;
            for (let k = 0; k < lay.hn; k++)
              if (lay.hlane[k] === lay.mlane[j] && lay.hx1[k] >= lay.mx0[j] && lay.hx0[k] <= lay.mx1[j])
                drawn.add(lay.keys[lay.hitem[k]]);
          }
          const x = (lay.mx0[j] + lay.mx1[j]) / 2;
          const y = lay.laneY + (lay.mlane[j] + 0.5) * lay.laneHeight;
          const key = hitView(f.layout, x, y)?.element?.key;
          expect(key !== undefined && drawn.has(key), `round ${round} ${lay.mode} mark ${j}`).toBe(true);
        }
      }
    }
    // the seed reaches runs of several items
    expect(merged).toBeGreaterThan(0);
  });
});

describe('the lane index', () => {
  it('finds the mark under x, else the nearest within the slack (a merged run’s gaps), never past it', () => {
    const index = new LaneIndex();
    // lane 0: [0, 10] and [10.5, 20]; lane 1: [5, 6]
    index.build(3, [0, 0, 1], [0, 10.5, 5], [10, 20, 6], 2);
    expect(index.at(0, 5)).toBe(0);
    expect(index.at(0, 10.25)).toBe(-1);
    expect(index.near(0, 10.25, 1)).toBe(0);
    expect(index.near(0, 10.4, 1)).toBe(1);
    expect(index.near(0, 21.5, 1)).toBe(-1);
    expect(index.near(0, 20.5, 1)).toBe(1);
    expect(index.near(1, 7, 0)).toBe(-1);
    expect(index.near(1, 6.5, 1)).toBe(2);
    expect(index.near(2, 5, 1)).toBe(-1);
    expect(index.near(-1, 5, 1)).toBe(-1);
  });

  it('agrees with a scan over seeded random marks', () => {
    const rng = new PCG(9);
    for (let round = 0; round < 50; round++) {
      const n = rng.int() % 60;
      const lanes = 1 + (rng.int() % 4);
      const lane = Array.from({ length: n }, () => rng.int() % lanes);
      const x0 = Array.from({ length: n }, () => (rng.int() % 4000) / 10);
      const x1 = x0.map((x) => x + (rng.int() % 200) / 10);
      const index = new LaneIndex();
      index.build(n, lane, x0, x1, lanes);
      for (let q = 0; q < 40; q++) {
        const l = rng.int() % lanes;
        const x = (rng.int() % 4400) / 10;
        const slack = (rng.int() % 3) / 2;
        const j = index.near(l, x, slack);
        const dist = (i: number) => (x < x0[i] ? x0[i] - x : x > x1[i] ? x - x1[i] : 0);
        let best = Number.POSITIVE_INFINITY;
        for (let i = 0; i < n; i++) if (lane[i] === l) best = Math.min(best, dist(i));
        if (best > slack) expect(j).toBe(-1);
        else {
          expect(j).toBeGreaterThanOrEqual(0);
          expect(lane[j]).toBe(l);
          expect(dist(j)).toBe(best);
        }
      }
    }
  });
});

describe('draw lists: emphasis', () => {
  const emphasis = {
    selection: rangesLocation('/data/p53', [3, 6, 10, 12]),
    cursor: { system: '/data/p53', position: 20 },
  };

  it('column (the default): a band a range over the section’s tracks and the gaps between them, one outline', async () => {
    const s = await settle(p53TracksSpec());
    const input = { size: Size, viewport: { start: 0, end: 30 } };
    await s.frame(input);
    const f = s.renderer.frame({ ...input, emphasis });
    expect(f.rebuilt).toBe(false);
    // the tracks run from 8 to 118 (20, 50 and 32 px with two gaps of 4: the domains' two lanes and a label lane): two
    // ranges, a tint and an outline each, the tints first so no tint covers an outline; positions 3 to 6 start at
    // 96 + 3 · 20; squircle ends, a thin outline
    // inside the tint at 0.6 of the accent
    expect(recordDrawLists(f.lists, ['selection']).split('\n').slice(2)).toEqual([
      'rect 156 8 60 110 #8f4c062e squircle 8',
      'rect 296 8 40 110 #8f4c062e squircle 8',
      'rect 156.5 8.5 59 109 none stroke #8f4c0699 1 squircle 7.5',
      'rect 296.5 8.5 39 109 none stroke #8f4c0699 1 squircle 7.5',
      '',
    ]);
    // the cursor: one thin outline of its cell over the band, in the neutral ink; a squircle by the cell over the band
    // (capped at 8), inside it by half its width
    expect(recordDrawLists(f.lists, ['cursor']).split('\n').slice(2)).toEqual([
      'rect 496.5 8.5 19 109 none stroke #2d271eff 1 squircle 7.5',
      '',
    ]);
  });

  it('tracks: a box a track, as before', async () => {
    const s = await settle(p53TracksSpec(), [], { emphasis: { selectionShape: 'tracks' } });
    const input = { size: Size, viewport: { start: 0, end: 30 } };
    await s.frame(input);
    const f = s.renderer.frame({ ...input, emphasis });
    // three tracks, two ranges each: a tint and an outline
    expect(f.lists.passes.selection.rects.n).toBe(12);
    expect(f.lists.passes.cursor.rects.n).toBe(3);
    // each track's tints inside its well's clip: the well's shape (_spec/README.md "Render", Shapes)
    const lines = recordDrawLists(f.lists, ['selection']).split('\n');
    expect(lines[2]).toBe('clip 96 8 600 20 squircle 6.67');
    expect(lines[3]).toBe('rect 156 8 60 20 #8f4c062e squircle 6.67');
  });
});

describe('draw lists: pure and keyed', () => {
  it('the same inputs give the same lists', async () => {
    const a = await renderSvg(p53TracksSpec(), { size: Size, viewport: { start: 10, end: 60 } });
    const b = await renderSvg(p53TracksSpec(), { size: Size, viewport: { start: 10, end: 60 } });
    expect(recordDrawLists(a.lists)).toBe(recordDrawLists(b.lists));
    expect(a.svg).toBe(b.svg);
  });

  it('an equal key reuses the layout (a hover redraws the emphasis alone); the token table’s version invalidates', async () => {
    const s = await settle(p53TracksSpec());
    const input = { size: Size, viewport: { start: 0, end: 50 } };
    const first = await s.frame(input);
    const key = first.key;
    const background = recordDrawLists(first.lists, ['background']);

    const again = s.renderer.frame(input);
    expect(again.rebuilt).toBe(false);
    expect(again.key).toBe(key);
    const hover = {
      layer: '/views/main/protein/letters/letters',
      section: 'protein',
      position: 3,
      system: '/data/p53',
    };
    const hovered = s.renderer.frame({ ...input, emphasis: { hover } });
    expect(hovered.rebuilt).toBe(false);
    expect(hovered.lists.passes.hover.size).toBe(1);

    s.renderer.tokens = DarkTokens;
    const dark = s.renderer.frame(input);
    expect(dark.rebuilt).toBe(true);
    expect(dark.key).not.toBe(key);
    expect(recordDrawLists(dark.lists, ['background'])).not.toBe(background);

    // an equal table is the same version: back to the first key
    const colors = Object.fromEntries(Object.entries(LightTokens.colors).map(([k, v]) => [k, colorHex(v)])) as Record<
      SeqTokenName,
      string
    >;
    const light = createTokens('light', colors);
    expect(light.version).toBe(LightTokens.version);
    expect(DarkTokens.version).not.toBe(LightTokens.version);
    s.renderer.tokens = light;
    const back = s.renderer.frame(input);
    expect(back.rebuilt).toBe(true);
    expect(back.key).toBe(key);
    expect(recordDrawLists(back.lists, ['background'])).toBe(background);
  });

  it('an entry arriving changes the key: a frame after it lays out again', async () => {
    const s = await settle(p53TracksSpec());
    const input = { size: Size, viewport: { start: 0, end: 30 } };
    const before = s.renderer.frame(input);
    // nothing was wanted yet: the layers wait, drawn as pending
    expect(blocks(before.lists)[0].letters).toBe(0);
    s.store.want(before.wants);
    await s.store.whenSettled();
    const after = s.renderer.frame(input);
    expect(after.rebuilt).toBe(true);
    expect(after.key).not.toBe(before.key);
    expect(blocks(after.lists)[0].letters).toBe(30);
  });
});
