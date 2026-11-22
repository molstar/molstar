/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Arcs on one track while the view moves (_spec/README.md "Render", The pairs layer, Arcs), a junction track's shape:
// short pairs between long ones, spaced wider than the window. The arcs' scale is the widest span among the pairs
// within reach of the window (an end in it, or both within `reach` past one edge), so (2) panned a pixel at a time, no
// arc's points jump as a long pair enters: before the fix the scale counted only the pairs drawn, and a long pair's
// first frame halved every arc at once. (1) A pair across the whole window is the widest arc over it, its feet `reach`
// past both edges; an end past the window is never a caret. (3) A pan's frame equals a fresh renderer's.

import { describe, expect, it } from 'vitest';
import { drawListsToSvg, type Frame, type PairsLay } from '../index.js';
import { arcsSpec } from './pairs-fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 900, height: 300 };
const N = 2000;
const Span = 120;
/** Every 200 positions a long pair (40) or a short one (10), a tiny one (5) halfway; and one across them all. */
const Pairs: [number, number][] = [];
for (let i = 0; i < 9; i++) {
  Pairs.push([200 * i, 200 * i + (i % 2 ? 10 : 40)]);
  Pairs.push([200 * i + 100, 200 * i + 105]);
}
const Across: [number, number] = [20, 1990];

const settled = (pairs: [number, number][]) => settle(arcsSpec({ n: N, pairs, height: 60 }));
const Junctions = settled(Pairs);
const WithAcross = settled([...Pairs, Across]);

function pairsLay(f: Frame): PairsLay {
  for (const l of f.layout.sections[0].sectionLays) if (l.type.name === 'pairs') return l.lay as PairsLay;
  throw new Error('no pairs layer');
}

function curves(f: Frame): Map<number, number[]> {
  const lay = pairsLay(f);
  const out = new Map<number, number[]>();
  for (let i = 0; i < lay.n; i++) out.set(lay.pair[i], Array.from(lay.curve.subarray(8 * i, 8 * i + 8)));
  return out;
}

/** Samples along a curve (a cubic's eight coordinates): x and y at t = i / 20. */
function samples(c: readonly number[]): number[] {
  const out: number[] = [];
  for (let i = 0; i <= 20; i++) {
    const t = i / 20;
    const u = 1 - t;
    for (let d = 0; d < 2; d++)
      out.push(u * u * u * c[d] + 3 * u * u * t * c[2 + d] + 3 * u * t * t * c[4 + d] + t * t * t * c[6 + d]);
  }
  return out;
}

/** The most any sample of a curve moved from another's. */
function moved(c: readonly number[], p: readonly number[]): number {
  const a = samples(c);
  const b = samples(p);
  return Math.max(...a.map((v, j) => Math.abs(v - b[j])));
}

describe('arcs on a track while the view moves', () => {
  it('(1) a pair across the window: the widest arc over it, its feet `reach` past both edges; no caret', async () => {
    const s = await WithAcross;
    for (const start of [300, 750, 1210]) {
      const f = await s.frame({ size: Size, viewport: { start, end: start + Span } });
      const lay = pairsLay(f);
      const ax = f.layout.sections[0].axis!;
      expect(lay.carets).toEqual([]);
      const c = curves(f).get(Pairs.length)!;
      expect(c[0]).toBeCloseTo(ax.x0 - 0.25 * ax.width, 6);
      expect(c[6]).toBeCloseTo(ax.x0 + ax.width + 0.25 * ax.width, 6);
      // the widest: its control the highest of all
      for (const [, d] of curves(f)) expect(c[3]).toBeLessThanOrEqual(d[3] + 1e-9);
      // canvas and SVG replay one list: the SVG has its curve
      const drawn = [...drawListsToSvg(f.lists).matchAll(/d="M ([^"C]*) C ([^"C]*)"/g)].map((m) =>
        `${m[1]} ${m[2]}`.split(' ').map(Number),
      );
      expect(drawn.some((q) => q.every((v, j) => Math.abs(v - c[j]) < 0.02))).toBe(true);
    }
  });

  it('(2) panned a pixel at a time as long pairs enter and leave: no arc’s points jump', async () => {
    const s = await Junctions;
    const ax = (await s.frame({ size: Size, viewport: { start: 0, end: Span } })).layout.sections[0].axis!;
    const step = Span / ax.width;
    let prev: Map<number, number[]> | undefined;
    for (let start = 240; start <= 330; start += step) {
      const next = curves(await s.frame({ size: Size, viewport: { start, end: start + Span } }));
      if (prev)
        for (const [k, c] of next) {
          const p = prev.get(k);
          if (p) expect(moved(c, p), `${k} at ${start}`).toBeLessThanOrEqual(2);
        }
      prev = next;
    }
  });

  it('(3) a pan’s frame equals a fresh renderer’s', async () => {
    const s = await Junctions;
    const fresh = await settled(Pairs);
    for (let start = 150; start <= 600; start += 41) {
      const viewport = { start, end: start + Span };
      expect(curves(await s.frame({ size: Size, viewport }))).toEqual(
        curves(await fresh.frame({ size: Size, viewport })),
      );
    }
  });

  it('(4) held pairs past the gate: the gate counts what the arcs reach, so a window across the threshold never bins by device column', async () => {
    // 40 pairs in one stretch, 30 across the whole sequence: an end in view is under the gate (50), what is drawn is not
    const pairs: [number, number][] = [];
    for (let i = 0; i < 40; i++) pairs.push([5000 + i, 5050 + i]);
    for (let i = 0; i < 30; i++) pairs.push([100 + i, 19000 + i]);
    const s = await settle(arcsSpec({ n: 20_000, pairs, height: 60 }), [], { pairs: { maxArcs: 50 } });
    const modes = new Set<string>();
    for (let start = 7500; start >= 2000; start -= 250) {
      const f = await s.frame({ size: Size, viewport: { start, end: start + 9000 } });
      const lay = pairsLay(f);
      const where = `at ${start} n=${lay.n} w=${!!lay.prep!.weight} d=${lay.prep!.data!.n}`;
      // items, or bins at a level: never the device-column bins past `maxArcs` (a density strip, a limit said)
      expect(lay.densityN, where).toBe(0);
      expect(
        f.layout.issues.map((i) => i.code),
        where,
      ).not.toContain('seq.render.pairs.limit');
      expect(lay.n, where).toBeLessThanOrEqual(50);
      modes.add(lay.prep!.weight ? 'bins' : 'items');
    }
    expect([...modes].sort()).toEqual(['bins', 'items']);
  });
});
