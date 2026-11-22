/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Arcs on a row of blocks while the view moves (_spec/README.md "Render", The pairs layer, Arcs): four sequences as
// blocks, the second scaled ×2, the third reversed, a gap between each two, pairs between every two of them (the last
// block's wholly out of view at the start). (1) A partner off the window, on a block partly in view or on one out of
// it, is an arc toward the edge it lies beyond, at most `reach` past it, never a caret; the line's strips share the
// pairs' band and `dy` (an end's strip changes nothing). (2) Panned a pixel at a time across each joint (entering
// and leaving), no arc's samples move more than the blocks themselves do (2 px at least). (3) An arc with a virtual end
// is hit, hovered and selected as drawn (the emphasis passes draw its curve). (4) An export at the panned window draws
// what the canvas path draws.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { describe, expect, it } from 'vitest';
import { type DrawLists, type Frame, hitView, type PairsLay, PlacePoint, recordDrawLists } from '../index.js';
import { renderSvg } from '../render-svg.js';
import { settle } from './helpers.js';
import { letters } from './pairs-fixtures.js';

const Size = { width: 900, height: 260 };
/** The blocks: their lengths and options; section units: 150, then 200 (×2), 150, 150. */
const Blocks = [
  { id: 'b1', n: 150 },
  { id: 'b2', n: 100, scale: 2 },
  { id: 'b3', n: 150, reversed: true },
  { id: 'b4', n: 150 },
] as const;
/** Each block's first section unit and its units a position. */
const First = [0, 150, 350, 500];
const Scale = [1, 2, 1, 1];
/** The joints between the blocks, in section units. */
const Joints = [150, 350, 500];
const Span = 100;
const Reach = 0.25;
/** The pairs: between each two blocks, by positions. */
const Links: [number, number, [number, number][]][] = [
  [
    0,
    1,
    [
      [140, 5],
      [100, 40],
      [20, 90],
      [149, 0],
    ],
  ],
  [
    1,
    2,
    [
      [95, 10],
      [60, 140],
      [10, 75],
    ],
  ],
  [
    2,
    3,
    [
      [5, 10],
      [140, 140],
      [75, 30],
    ],
  ],
  [
    0,
    3,
    [
      [10, 140],
      [130, 5],
    ],
  ],
  [0, 2, [[75, 75]]],
];

function spec() {
  const b = createSeqSpec();
  const seqs = Blocks.map((x, i) => b.data.sequence(`s${i}`, { alphabet: 'rna', letters: letters(x.n, 11 + i) }));
  const nodes = Links.map(([a, c, pairs], j) =>
    b.data.pairs(`p${j}`, {
      a: seqs[a],
      b: seqs[c],
      rows: {
        a_start: pairs.map((p) => p[0]),
        a_end: pairs.map((p) => p[0] + 1),
        b_start: pairs.map((p) => p[1]),
        b_end: pairs.map((p) => p[1] + 1),
        id: pairs.map((p) => `${a}:${p[0]}-${c}:${p[1]}`),
      },
    }),
  );
  b.views.view('v', (v) => {
    const s = v.section('row', {});
    Blocks.forEach((x, i) =>
      s.block(x.id, {
        of: seqs[i],
        ...('scale' in x ? { scale: x.scale } : {}),
        ...('reversed' in x ? { direction: 'reversed' as const } : {}),
      }),
    );
    s.track('letters').layer('letters', { type: 'letters' });
    s.track('pairs', { height: 60 }, (t) => {
      nodes.forEach((p, j) => t.layer(`p${j}`, { type: 'pairs', data: p }));
    });
  });
  return b.build();
}

const Spec = spec();
const Settled = settle(Spec);
const frameAt = async (start: number) =>
  (await Settled).frame({ size: Size, settled: true, viewport: { start, end: start + Span } });

function pairsLays(f: Frame): PairsLay[] {
  const out: PairsLay[] = [];
  for (const l of f.layout.sections[0].sectionLays) if (l.type.name === 'pairs') out.push(l.lay as PairsLay);
  return out;
}

/** A pair's ends' section units: its block's first unit, then its position (reversed: from the end) times the scale. */
function unitsOf(lay: PairsLay, k: number): [number, number] {
  const j = Number(lay.key.slice(lay.key.lastIndexOf('/p') + 2));
  const [a, b] = Links[j];
  const d = lay.prep!.data!;
  const unit = (blk: number, p: number) =>
    First[blk] + Scale[blk] * ('reversed' in Blocks[blk] ? Blocks[blk].n - p - 0.5 : p + 0.5);
  return [unit(a, d.aStart[k]), unit(b, d.bStart[k])];
}

function curves(f: Frame): Map<string, number[]> {
  const out = new Map<string, number[]>();
  for (const lay of pairsLays(f))
    for (let i = 0; i < lay.n; i++)
      out.set(`${lay.key}:${lay.pair[i]}`, Array.from(lay.curve.subarray(8 * i, 8 * i + 8)));
  return out;
}

/** A cubic's point at t. */
function pointAt(c: ArrayLike<number>, t: number): [number, number] {
  const u = 1 - t;
  const at = (d: number) =>
    u * u * u * c[d] + 3 * u * u * t * c[2 + d] + 3 * u * t * t * c[4 + d] + t * t * t * c[6 + d];
  return [at(0), at(1)];
}

function moved(c: readonly number[], p: readonly number[]): number {
  let m = 0;
  for (let i = 0; i <= 20; i++) {
    const a = pointAt(c, i / 20);
    const b = pointAt(p, i / 20);
    m = Math.max(m, Math.abs(a[0] - b[0]), Math.abs(a[1] - b[1]));
  }
  return m;
}

/** The curves a pass draws: each path's numbers. */
function passCurves(lists: DrawLists, pass: 'hover' | 'selection' | 'data'): number[][] {
  return recordDrawLists(lists, [pass])
    .split('\n')
    .filter((l) => l.startsWith('path M'))
    .map((l) => l.split(' ').map(Number).filter(Number.isFinite));
}

const has = (drawn: number[][], c: ArrayLike<number>) =>
  drawn.some((q) => q.length >= 8 && Array.from(c).every((v, j) => Math.abs(v - q[j]) < 0.02));

const Pt = new PlacePoint();

describe('arcs on a row of blocks (scaled, reversed, gapped) while the view moves', () => {
  it('(1) a partner off the window, on a block partly in view or out of it: an arc toward its edge, never a caret', async () => {
    let virtual = 0;
    for (let start = 0; start + Span <= 650; start += 10) {
      const f = await frameAt(start);
      const sec = f.layout.sections[0];
      const ax = sec.axis!;
      const s = sec.placement.strips;
      // one line: every strip's pairs band and `dy` alike
      for (let i = 1; i < s.n; i++) {
        expect(s.dy[i]).toBe(s.dy[0]);
        expect(s.bandY[i * s.bands + 1]).toBe(s.bandY[1]);
        expect(s.bandH[i * s.bands + 1]).toBe(s.bandH[1]);
      }
      for (const lay of pairsLays(f)) {
        const where = `${lay.key} at ${start}`;
        expect(lay.carets, where).toEqual([]);
        let reached = 0;
        for (let k = 0; k < lay.prep!.data!.n; k++) {
          const [a, b] = unitsOf(lay, k);
          if ((a >= ax.u0 && a <= ax.u1) || (b >= ax.u0 && b <= ax.u1) || a < ax.u0 !== b < ax.u0) reached++;
        }
        expect(lay.n, where).toBe(reached);
        for (let i = 0; i < lay.n; i++) {
          const c = lay.curve.subarray(8 * i, 8 * i + 8);
          const [a, b] = unitsOf(lay, lay.pair[i]);
          for (const [u, x] of [
            [Math.min(a, b), c[0]],
            [Math.max(a, b), c[6]],
          ]) {
            if (u >= ax.u0 && u <= ax.u1) {
              expect(sec.placement.place(u, Pt), where).toBe(true);
              expect(x, where).toBeCloseTo(Pt.x, 6);
            } else {
              virtual++;
              const edge = u < ax.u0 ? ax.x0 : ax.x0 + ax.width;
              const past = u < ax.u0 ? edge - x : x - edge;
              expect(past, where).toBeGreaterThanOrEqual(0);
              expect(past, where).toBeLessThanOrEqual(Reach * ax.width + 1e-6);
            }
          }
          expect(c[1] - pointAt(c, 0.5)[1], where).toBeGreaterThanOrEqual(2);
        }
      }
    }
    expect(virtual).toBeGreaterThan(50);
  });

  it('(2) panned a pixel at a time across each joint, entering and leaving: no arc jumps', async () => {
    const ax = (await frameAt(0)).layout.sections[0].axis!;
    const step = Span / ax.width;
    for (const joint of Joints)
      for (const from of [joint - Span - 3, joint - 3]) {
        let prev: Map<string, number[]> | undefined;
        let prevUnits: number[] = [];
        for (let start = Math.max(0, from); start <= from + 6; start += step) {
          const f = await frameAt(start);
          const next = curves(f);
          // the blocks' own move under the step (a gap growing as its boundary leaves the edge): arcs follow it
          const units: number[] = [];
          for (let u = Math.ceil(start) + 0.5; u < start + Span; u += 2)
            units.push(f.layout.sections[0].placement.place(u, Pt) ? Pt.x : Number.NaN);
          let shift = 0;
          for (let j = 0; j < units.length; j++)
            if (Number.isFinite(units[j]) && Number.isFinite(prevUnits[j]))
              shift = Math.max(shift, Math.abs(units[j] - prevUnits[j]));
          const bound = Math.max(2, shift + 0.25);
          if (prev)
            for (const [key, c] of next) {
              const p = prev.get(key);
              if (p) expect(moved(c, p), `${key} at ${start.toFixed(2)}`).toBeLessThanOrEqual(bound);
            }
          prev = next;
          prevUnits = units;
        }
      }
  });

  it('(3) an arc with a virtual end: hit, hovered and selected as drawn', async () => {
    const s = await Settled;
    const viewport = { start: 120, end: 120 + Span };
    const f = await s.frame({ size: Size, settled: true, viewport });
    const ax = f.layout.sections[0].axis!;
    let checked = 0;
    for (const lay of pairsLays(f))
      for (let i = 0; i < lay.n; i++) {
        const c = lay.curve.subarray(8 * i, 8 * i + 8);
        const off = c[0] < ax.x0 || c[6] > ax.x0 + ax.width;
        if (!off) continue;
        // a point of the curve inside the data box, off its feet
        let p: [number, number] | undefined;
        for (let t = 0.05; t < 1 && !p; t += 0.05) {
          const q = pointAt(c, t);
          if (q[0] > ax.x0 + 2 && q[0] < ax.x0 + ax.width - 2 && c[1] - q[1] > 3) p = q;
        }
        if (!p) continue;
        const hit = hitView(f.layout, p[0], p[1]);
        // an arc of another layer over it there: the nearest wins
        if (hit?.element?.owner !== lay.prep!.owner || hit?.element?.key !== lay.prep!.data!.keys[lay.pair[i]])
          continue;
        const hover = await s.frame({ size: Size, settled: true, viewport, emphasis: { hover: hit } });
        expect(has(passCurves(hover.lists, 'hover'), c)).toBe(true);
        const selection = await s.frame({
          size: Size,
          settled: true,
          viewport,
          emphasis: { elements: { keys: [hit.element!] } },
        });
        expect(has(passCurves(selection.lists, 'selection'), c)).toBe(true);
        checked++;
      }
    expect(checked).toBeGreaterThan(2);
  });

  it('(4) an export at a panned window draws what the canvas path draws', async () => {
    for (const start of [100, 330, 480]) {
      const viewport = { start, end: start + Span };
      const canvas = await frameAt(start);
      const svg = await renderSvg(Spec, { size: Size, viewport });
      expect(svg.settled).toBe(true);
      expect(passCurves(svg.lists, 'data')).toEqual(passCurves(canvas.lists, 'data'));
      expect(passCurves(canvas.lists, 'data').length).toBeGreaterThan(3);
    }
  });
});
