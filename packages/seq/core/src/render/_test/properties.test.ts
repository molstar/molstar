/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The draw lists' properties: bounded by the screen at every density step (a 3M-position system and p53 alike), hits
// where the marks are, pure (same inputs, same lists), keyed by what they read (the token table's version among it).

import { describe, expect, it } from 'vitest';
import {
  colorHex,
  createTokens,
  DarkTokens,
  type DensityStep,
  densityStep,
  DefaultRenderConfig,
  hitView,
  LightTokens,
  type LettersLay,
  type RangesLay,
  recordDrawLists,
  renderSvg,
  type SeqTokenName,
} from '../index.js';
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
const { inset, labelWidth, labelGap } = DefaultRenderConfig;
const Width = DataWidth + inset + labelWidth + labelGap + inset;
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
        expect(signal.rects).toBe(Math.ceil(n));
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
      expect(pv.rects).toBe(signal.rects);
      // ranges: the same bound, and a mark (with at most its label) a domain in view
      const inView = Domains.filter(([a, b]) => a < 100 + n && b > 100).length;
      expect(pr.rects).toBe(inView);
      expect(pr.rects + pr.images + pr.texts).toBeLessThanOrEqual(DataWidth * Lanes[1]);
    }
  });

  it('p53 whole is bounded too: a cell a position at its own step', async () => {
    const r = await renderSvg(p53TracksSpec(), { size: Size });
    const [letters, items, values] = blocks(r.lists);
    expect(letters.rects).toBe(P53.length);
    expect(items.rects).toBeLessThanOrEqual(6);
    expect(values.rects).toBe(P53.length);
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
      const y = ranges.box.y + (ranges.mlane[j] + 0.5) * ranges.laneHeight;
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
    const between = hitView(f.layout, x, ranges.box.y + 0.5 * ranges.laneHeight);
    expect(between).toEqual({
      layer: '/views/main/protein/domains:0',
      section: 'protein',
      position: 62,
      system: '/data/p53',
    });

    const letters = track(0).layers[0].lay as LettersLay;
    for (let j = 0; j < letters.n; j++) {
      const hit = hitView(f.layout, (letters.x0[j] + letters.x1[j]) / 2, letters.box.y + 1);
      expect(hit).toEqual({
        layer: '/views/main/protein/letters:0',
        section: 'protein',
        position: letters.pos[j],
        system: '/data/p53',
      });
    }
    // the label column and the page's margin hit nothing
    expect(hitView(f.layout, 4, letters.box.y + 1)).toBeUndefined();
    expect(hitView(f.layout, 200, 2)).toBeUndefined();
  });
});

describe('draw lists: emphasis', () => {
  it('the selection tints and outlines its ranges on every track of its system; the cursor outlines a cell', async () => {
    const s = await settle(p53TracksSpec());
    const input = { size: Size, viewport: { start: 0, end: 30 } };
    await s.frame(input);
    const f = s.renderer.frame({
      ...input,
      emphasis: {
        selection: rangesLocation('/data/p53', [3, 6, 10, 12]),
        cursor: { system: '/data/p53', position: 20 },
      },
    });
    expect(f.rebuilt).toBe(false);
    // three tracks, two ranges each: a tint and an outline; positions 3 to 6 start at 128 + 3 · 20 (the cursor's
    // cell, 20, at 128 + 20 · 20)
    expect(f.lists.passes.selection.rects.n).toBe(12);
    expect(recordDrawLists(f.lists, ['selection']).split('\n').slice(2, 4)).toEqual([
      'rect 188 8 60 20 #8f4c062e',
      'rect 188.75 8.75 58.5 18.5 none stroke #8f4c06ff 1.5',
    ]);
    expect(f.lists.passes.cursor.rects.n).toBe(3);
    const lines = recordDrawLists(f.lists, ['cursor']).split('\n');
    expect(lines[2]).toBe('rect 529.25 9.25 17.5 17.5 none stroke #8f4c06ff 2.5');
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
    const hover = { layer: '/views/main/protein/letters:0', section: 'protein', position: 3, system: '/data/p53' };
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
