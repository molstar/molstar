/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A value's color at full strength, a line's dots and the `points` encoding (_spec/README.md "Render", Layer types,
// `values`): a color from the value is drawn as the scheme gives it (the legend's swatches), a uniform series keeps the
// tints; a line and an area show a dot a value from 6 px a value; `points` is dots only, a thin mark a device column
// below a pixel a value; a missing value has neither; the export draws the canvas's dots.

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { hitView, PlainDatum, type ValuesLay } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';
import { layerLay, type ReadPath, readPaths, readRects } from './small-helpers.js';

const Size = { width: 400, height: 56 };

function spec(values: readonly number[], layer: Record<string, unknown> = {}) {
  const b = createSeqSpec({ title: 'points' });
  const s = b.data.sequence('s', { alphabet: 'dna', letters: 'A'.repeat(values.length) });
  const v = b.data.values('v', { target: s, values: Float32Array.from(values) });
  b.views.view('main', (view) =>
    view
      .section('s', { of: s })
      .track('v')
      .layer('values', { type: 'values', data: v, ...layer }),
  );
  return b.build();
}

async function frame(values: readonly number[], layer: Record<string, unknown> = {}, width = Size.width) {
  const s = await settle(spec(values, layer));
  const f = await s.frame({ size: { ...Size, width } });
  return { s, f, lay: layerLay<ValuesLay>(f, 0) };
}

/** The dots: closed paths of the polygon's 16 corners, a ring (the well's ground) then the dot, a value. */
const dots = (r: ReadPath[]) => r.filter((p) => p.closed && p.xy.length === 32);
const width = (p: ReadPath) => {
  const xs = p.xy.filter((_, k) => k % 2 === 0);
  return Math.max(...xs) - Math.min(...xs);
};

const Ramp = theme.field('value', { domain: [0, 10] });
const Values = [1, 3, 5, 7, 9, 2];
const Wave = Array.from({ length: 1200 }, (_, i) => 5 + 4 * Math.sin(i / 90));

describe('a color from the value', () => {
  it('bars take the scheme’s colors exactly, at full strength (the legend’s swatches)', async () => {
    const { f, lay } = await frame(Values, { color: Ramp });
    expect(lay.perBar).toBe(true);
    const d = new PlainDatum();
    const rects = readRects(f.lists).filter((x) => x.w > 20 && x.h > 3);
    expect(rects.length).toBe(Values.length);
    rects.forEach((x, i) => {
      d.value = Values[i];
      expect(x.fill >>> 0).toBe((lay.color!.at(d) as number) >>> 0);
      expect(x.fill & 255).toBe(255);
    });
  });

  it('a uniform series keeps the tints: bars at 60%, the profile’s area at 35%', async () => {
    const bars = await frame(Values);
    const alphas = readRects(bars.f.lists)
      .filter((x) => x.w > 20 && x.h > 3)
      .map((x) => x.fill & 255);
    expect(alphas).toEqual(Values.map(() => 0x99));
    const profile = readPaths((await frame(Wave)).f.lists).filter((p) => p.closed);
    expect(profile.length).toBeGreaterThan(0);
    for (const p of profile) expect(p.fill & 255).toBe(0x59);
  });

  it('a profile’s area from the value is nearly opaque, its line full', async () => {
    const paths = readPaths((await frame(Wave, { color: Ramp })).f.lists);
    const areas = paths.filter((p) => p.closed);
    expect(areas.length).toBeGreaterThan(0);
    for (const p of areas) expect(p.fill & 255).toBe(230);
    for (const p of paths.filter((q) => !q.closed)) expect(p.stroke & 255).toBe(255);
  });
});

describe('a line’s dots', () => {
  // 520 px, 104 of them the labels': 66 positions are 6.3 px each, 74 are 5.6
  it('appear at 6 px a value, not below', async () => {
    const at6 = await frame(
      Array.from({ length: 66 }, (_, i) => 1 + (i % 5)),
      { encoding: 'line' },
      520,
    );
    expect(at6.lay.axis.scale).toBeGreaterThanOrEqual(6);
    expect(dots(readPaths(at6.f.lists)).length).toBe(2 * 66);
    const under = await frame(
      Array.from({ length: 74 }, (_, i) => 1 + (i % 5)),
      { encoding: 'line' },
      520,
    );
    expect(under.lay.axis.scale).toBeLessThan(6);
    expect(dots(readPaths(under.f.lists))).toEqual([]);
  });

  it('are the value’s color in a ring of the well’s ground, over an area too', async () => {
    const { f, lay } = await frame(Values, { encoding: 'area', color: Ramp });
    const d = dots(readPaths(f.lists));
    expect(d.length).toBe(2 * Values.length);
    const ground = f.layout.tokens.colors['seq-track'];
    const x = new PlainDatum();
    for (let i = 0; i < Values.length; i++) {
      x.value = Values[i];
      expect(d[2 * i].fill >>> 0).toBe(ground >>> 0);
      expect(d[2 * i + 1].fill >>> 0).toBe((lay.color!.at(x) as number) >>> 0);
      // 3 px across, its ring 1 px beyond it
      expect(width(d[2 * i + 1])).toBeCloseTo(3, 3);
      expect(width(d[2 * i])).toBeCloseTo(5, 3);
    }
  });

  it('a missing value has no dot', async () => {
    const { f } = await frame([1, 2, Number.NaN, 3, 2, 1], { encoding: 'line' });
    expect(dots(readPaths(f.lists)).length).toBe(2 * 5);
  });
});

describe('the points encoding', () => {
  it('is dots only: no line, no area; a missing value has none', async () => {
    const { f } = await frame([1, 2, Number.NaN, 3, 2, 1], { encoding: 'points' });
    const paths = readPaths(f.lists);
    expect(paths.length).toBe(2 * 5);
    expect(dots(paths).length).toBe(2 * 5);
  });

  it('below a pixel a value is a thin vertical mark of the column’s extent, min to max, no dots', async () => {
    const values = Array.from({ length: 2000 }, (_, i) => (i % 2 ? 8 : 2));
    const { f, lay } = await frame(values, { encoding: 'points', domain: [0, 10] });
    expect(readPaths(f.lists)).toEqual([]);
    const marks = readRects(f.lists).filter((x) => x.w < 1.5 && x.h > 2);
    expect(marks.length).toBe(lay.n);
    // each column holds both values: the mark covers [2, 8]
    const [lo, hi] = lay.domain!;
    const y = (v: number) => lay.box.y + lay.box.height - ((v - lo) / ((hi - lo) / 0.9)) * lay.box.height;
    for (const m of marks.slice(5, 50)) {
      expect(m.y).toBeCloseTo(y(8), 3);
      expect(m.y + m.h).toBeCloseTo(y(2), 3);
    }
  });

  it('a hit is the nearest value of the column: the position under the pointer', async () => {
    const { f, lay } = await frame(Values, { encoding: 'points' });
    const w = lay.box.width / Values.length;
    for (let i = 0; i < Values.length; i++) {
      const h = hitView(f.layout, lay.box.x + (i + 0.5) * w, lay.box.y + lay.box.height / 2);
      expect(h?.position).toBe(i);
      expect(h?.layer).toBe(lay.key);
    }
  });

  it('hover rings the hovered dot', async () => {
    const { s, f, lay } = await frame(Values, { encoding: 'points' });
    const h = hitView(f.layout, lay.box.x + 1.5 * (lay.box.width / Values.length), lay.box.y + 20)!;
    const g = await s.frame({ size: Size, emphasis: { hover: h } });
    const ring = readPaths(g.lists, 'hover');
    expect(ring.length).toBe(1);
    expect(ring[0].fill === 0 && ring[0].strokeWidth > 0).toBe(true);
  });
});

describe('the export draws the canvas’s dots', () => {
  it.each(['line', 'area', 'points'])('%s: a path a dot ring and a dot, as the lists have them', async (encoding) => {
    const r = await renderSvg(spec([1, 2, Number.NaN, 3, 2, 1], { encoding }), { size: Size });
    const lists = readPaths(r.lists);
    expect(dots(lists).length).toBe(2 * 5);
    // the wells' and clips' paths aside, the data pass's are the lists' (a path an element)
    expect((r.svg.match(/<path d="/g) ?? []).length).toBeGreaterThanOrEqual(lists.length);
    for (const p of dots(lists)) expect(r.svg).toContain(`fill="#${(p.fill >>> 8).toString(16).padStart(6, '0')}`);
  });
});
