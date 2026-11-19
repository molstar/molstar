/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A one-position block is a chip in every layout: a fixed slot its label's width (5 px either side), never a
// position's share, at the fit view and zoomed; its letters' band draws its label in a pill, never clipped (canvas and
// SVG, both themes); the row's x skips it as it skips a gap; hit, hover and a selection made elsewhere read it as its
// one position; a zoom about it keeps it whole.

import { rangesLocation } from '@molstar/seqx-core/query';
import { describe, expect, it } from 'vitest';
import {
  type BlockPlacement,
  bundledTextMetrics,
  DarkTokens,
  drawListsToSvg,
  type Frame,
  hitView,
  LightTokens,
  PlacePoint,
  type RingPlacement,
  recordDrawLists,
  TracksLayout,
} from '../index.js';
import { settle } from './helpers.js';
import { blocksSpec } from './rings-fixtures.js';

const Size = { width: 800, height: 600 };
const Atp = '/data/atp';
const LettersKey = '/views/v/complex/letters/letters';
/** The chip's width the rule gives: its label in the label font, 5 px either side. */
const chipWidth = (label: string) => bundledTextMetrics().measure('12px sans', label) + 10;

/** The chip's strip (ATP's). */
const chipOf = (p: BlockPlacement | RingPlacement) => p.strips.axes.findIndex((a) => a.system === Atp);

/** The data pass's texts as records (`text "ATP" x y middle font color width w max m`). */
const dataTexts = (f: Frame) =>
  recordDrawLists(f.lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('text '));

/** The chip's pill and its label drawn whole: the label's measured width within its max width, the pill's. */
function checkPill(f: Frame, p: BlockPlacement, label: string) {
  const i = chipOf(p);
  const a = p.strips.axes[i];
  const line = dataTexts(f).find((l) => l.startsWith(`text "${label}"`));
  expect(line).toBeDefined();
  const width = Number(/ width ([\d.]+)/.exec(line!)![1]);
  const max = Number(/ max ([\d.]+)/.exec(line!)![1]);
  expect(width).toBeGreaterThan(0);
  expect(width).toBeLessThanOrEqual(max);
  expect(max).toBeCloseTo(a.width, 0);
  // the pill: a rect with round ends over the chip's letters band
  const k = i * p.strips.bands + 1;
  const pill = recordDrawLists(f.lists, ['data'])
    .split('\n')
    .filter((l) => l.startsWith('rect ') && l.includes('squircle') && !/ (left|right)$/.test(l))
    .map((l) => l.split(' ').slice(1, 5).map(Number))
    .find(([x, y]) => Math.abs(x - a.x0 - 1) < 1 && Math.abs(y - p.strips.bandY[k] - 1) < 1);
  expect(pill).toBeDefined();
  // the label is not said again above the chip
  const labels = recordDrawLists(f.lists, ['labels'])
    .split('\n')
    .filter((l) => l.startsWith(`text "${label}"`));
  expect(labels).toEqual([]);
}

describe('a chip in the tracks layout', () => {
  for (const tokens of [LightTokens, DarkTokens])
    it(`a row (${tokens === LightTokens ? 'light' : 'dark'}): the chip its label's width at the fit view and zoomed, drawn whole in a pill`, async () => {
      const s = await settle(blocksSpec());
      s.renderer.tokens = tokens;
      for (const viewport of [undefined, { start: 150, end: 201 }, { start: 180, end: 201 }]) {
        const f = await s.frame({ size: Size, ...(viewport ? { viewport } : {}) });
        const p = f.layout.sections[0].placement as BlockPlacement;
        const i = chipOf(p);
        expect(p.strips.chip[i]).toBe(1);
        expect(p.strips.axes[i].width).toBeCloseTo(chipWidth('ATP'), 6);
        expect([p.strips.axes[i].u0, p.strips.axes[i].u1]).toEqual([0, 1]);
        checkPill(f, p, 'ATP');
        // the polymers fill the rest of the width: the last strip ends at the data box's end
        const axis = f.layout.sections[0].axis!;
        expect(p.strips.axes[i].x0 + p.strips.axes[i].width).toBeCloseTo(axis.x0 + axis.width, 6);
        // the SVG draws the label too
        expect(drawListsToSvg(f.lists)).toContain('>ATP<');
      }
    });

  it('a row: the x skips the chip as it skips a gap; placement the inverse of hit-testing over it', async () => {
    const s = await settle(blocksSpec());
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement as BlockPlacement;
    const axis = f.layout.sections[0].axis!;
    // the polymers' scale: the width less the gaps and the chip, over the units less the chip's
    const gaps = p.strips.axes[1].x0 - (p.strips.axes[0].x0 + p.strips.axes[0].width);
    const scale = (axis.width - 2 * gaps - chipWidth('ATP')) / 200;
    expect(p.strips.axes[0].scale).toBeCloseTo(scale, 6);
    expect(p.strips.axes[1].scale).toBeCloseTo(scale, 6);
    const pt = new PlacePoint();
    for (const u of [0, 119, 120, 199, 200]) {
      expect(p.place(u + 0.5, pt)).toBe(true);
      expect(p.at([pt.x, pt.y], pt) && Math.floor(pt.offset)).toBe(u);
    }
    // anywhere across the chip is its one position
    const a = p.strips.axes[2];
    for (const x of [a.x0 + 1, a.x0 + a.width / 2, a.x0 + a.width - 1]) {
      expect(p.at([x, p.strips.bandY[2 * p.strips.bands + 1] + 2], pt)).toBe(true);
      expect(Math.floor(pt.offset)).toBe(200);
    }
  });

  it('a row: the chip hit, hovered and selected as its one position', async () => {
    const s = await settle(blocksSpec());
    const input = { size: Size };
    const f = await s.frame(input);
    const p = f.layout.sections[0].placement as BlockPlacement;
    const i = chipOf(p);
    const a = p.strips.axes[i];
    const y = p.strips.bandY[i * p.strips.bands + 1] + 2;
    for (const x of [a.x0 + 2, a.x0 + a.width - 2])
      expect(hitView(f.layout, x, y)).toMatchObject({ block: i, system: Atp, position: 0 });
    // the hover outlines the whole chip
    const hover = { layer: LettersKey, section: 'complex', system: Atp, position: 0 };
    const h = s.renderer.frame({ ...input, emphasis: { hover } });
    const outline = recordDrawLists(h.lists, ['hover'])
      .split('\n')
      .filter((l) => l.startsWith('rect '))
      .map((l) => l.split(' ').slice(1, 4).map(Number));
    expect(outline.length).toBe(1);
    expect(outline[0][0]).toBeLessThan(a.x0);
    expect(outline[0][0] + outline[0][2]).toBeGreaterThan(a.x0 + a.width);
    // a selection made elsewhere (a link, a structure) on ATP lights the chip
    const sel = s.renderer.frame({ ...input, emphasis: { selection: rangesLocation(Atp, [0, 1]) } });
    const lit = recordDrawLists(sel.lists, ['selection'])
      .split('\n')
      .filter((l) => l.startsWith('rect '))
      .map((l) => l.split(' ').slice(1, 4).map(Number));
    expect(lit.some(([x, , w]) => x <= a.x0 + 1 && x + w >= a.x0 + a.width - 1)).toBe(true);
  });

  it('a zoom about the chip keeps it whole, its width the same', async () => {
    const s = await settle(blocksSpec());
    const f0 = await s.frame({ size: Size });
    const p0 = f0.layout.sections[0].placement as BlockPlacement;
    const a0 = p0.strips.axes[chipOf(p0)];
    const at = p0.offsetAt(a0.x0 + a0.width / 2, p0.strips.bandY[1] + 2);
    expect(Math.floor(at)).toBe(200);
    const axis = f0.layout.sections[0].axis!;
    const c = { total: 201, width: axis.width, height: 600, cyclic: false, maxScale: 40, placement: p0 };
    const v = TracksLayout.nav.zoom(c, { start: 0, end: 201 }, 4, at);
    expect(v.end - v.start).toBeCloseTo(201 / 4, 6);
    const f = await s.frame({ size: Size, viewport: v });
    const p = f.layout.sections[0].placement as BlockPlacement;
    const a = p.strips.axes[chipOf(p)];
    expect(a.width).toBeCloseTo(a0.width, 6);
    expect(p.offsetAt(a.x0 + a.width / 2, p.strips.bandY[1] + 2)).toBeCloseTo(at, 6);
    checkPill(f, p, 'ATP');
  });

  it('split copies: each copy’s blocks, the chip its label’s width', async () => {
    const s = await settle(blocksSpec({ copies: 'each' }));
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement as BlockPlacement;
    expect(p.strips.n).toBe(4);
    expect(p.strips.axes[chipOf(p)].width).toBeCloseTo(chipWidth('ATP'), 6);
    checkPill(f, p, 'ATP');
  });

  it('a stack: the chips on a row of their own, their label’s width, whatever the window', async () => {
    const s = await settle(blocksSpec({ arrange: 'stack' }));
    for (const viewport of [undefined, { start: 60, end: 100 }]) {
      const f = await s.frame({ size: Size, ...(viewport ? { viewport } : {}) });
      const p = f.layout.sections[0].placement as BlockPlacement;
      const i = chipOf(p);
      expect(i).toBe(2);
      // below the polymers' rows, at the first position's x
      expect(p.strips.top[i]).toBeGreaterThan(p.strips.top[1]);
      expect(p.strips.axes[i].x0).toBeCloseTo(p.strips.axes[0].x0, 6);
      expect(p.strips.axes[i].width).toBeCloseTo(chipWidth('ATP'), 6);
      checkPill(f, p, 'ATP');
      const a = p.strips.axes[i];
      expect(hitView(f.layout, a.x0 + a.width / 2, p.strips.bandY[i * p.strips.bands + 1] + 2)).toMatchObject({
        block: i,
        system: Atp,
        position: 0,
      });
    }
  });
});

describe('a chip on a ring', () => {
  it('a fixed-angle arc its label’s width on the outer ring, after the polymers’ arcs', async () => {
    const s = await settle(blocksSpec({ layout: 'rings' }));
    const f = await s.frame({ size: Size });
    const p = f.layout.sections[0].placement as RingPlacement;
    const i = chipOf(p);
    expect(i).toBe(2);
    expect(p.strips.chip[i]).toBe(1);
    expect(p.strips.axes[i].width).toBeCloseTo(chipWidth('ATP'), 6);
    // the polymers share the rest of the circle by their lengths
    const [a, b] = p.strips.axes;
    expect(a.scale).toBeCloseTo(b.scale, 6);
    expect(p.strips.axes[i].x0).toBeGreaterThan(b.x0 + b.width);
    const pt = new PlacePoint();
    p.place(200.5, pt);
    expect(p.at([pt.x, pt.y], pt) && Math.floor(pt.offset)).toBe(200);
  });
});
