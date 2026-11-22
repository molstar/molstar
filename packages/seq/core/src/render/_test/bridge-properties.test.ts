/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge's properties (_spec/README.md "Render", Bridges): the twist by geometry equals GA's three-direction rule over
// every combination of block and section direction (EC-206); every ribbon's ends at its data's positions on reversed
// and scaled blocks; the closed-form hit is the drawn polygon's topmost cover (random ribbons, both directions, off-window
// ends, stubs).

import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { blockUnitOf } from '@molstar/seqx-core/view';
import { hitView } from '../index.js';
import { BlockPlacement } from '../layout/blocks.js';
import type { RibbonsLay } from '../pairs-ribbons.js';
import { type BridgePair, bridgeSpec, type BridgeSpecOptions } from './bridge-fixtures.js';
import { settle } from './helpers.js';

const Size = { width: 900, height: 600 };

async function frame(o: BridgeSpecOptions, config = {}, input: Record<string, unknown> = {}) {
  const s = await settle(bridgeSpec(o), [], config);
  const f = await s.renderer.frame({ size: Size, settled: true, ...input });
  return { f, lay: f.layout.bridges!.list[0].lays[0].lay as RibbonsLay };
}

/** The page x of position `p` on side `i`'s block `b` (a row of blocks, as laid). */
function blockX(f: Awaited<ReturnType<typeof frame>>['f'], i: number, b: number, p: number): number {
  const s = f.layout.sections[i];
  const u = blockUnitOf(s.section.blocks![b], p)!;
  return (s.placement as BlockPlacement).row.xOf(u);
}

describe('a ribbon’s twist (EC-206)', () => {
  it('is crossed exactly when run × display direction of A × of B is −1, over every combination', async () => {
    for (const run of [1, -1])
      for (const topSection of [false, true])
        for (const topBlock of [false, true])
          for (const bottomSection of [false, true])
            for (const bottomBlock of [false, true]) {
              const { lay } = await frame({
                pairs: [{ a: [40, 120], b: [60, 140], strand: run > 0 ? '+' : '-' }],
                top: [{ range: [0, 400], reversed: topBlock }],
                bottom: [{ range: [0, 300], reversed: bottomBlock }],
                topReversed: topSection,
                bottomReversed: bottomSection,
              });
              const dirA = (topSection ? -1 : 1) * (topBlock ? -1 : 1);
              const dirB = (bottomSection ? -1 : 1) * (bottomBlock ? -1 : 1);
              const r = lay.rib;
              const crossed = Math.sign(r[1] - r[0]) !== Math.sign(r[4] - r[3]);
              expect(crossed, `run ${run}, A ${dirA}, B ${dirB}`).toBe(run * dirA * dirB === -1);
            }
  });
});

describe('a ribbon’s ends against its data', () => {
  it('are at its positions on reversed and scaled blocks, every place of a block showing it', async () => {
    const top = [{ range: [0, 150] }, { range: [150, 400], reversed: true }] as const;
    const bottom = [
      { range: [0, 100], scale: 3 },
      { range: [100, 300], reversed: true, scale: 0.5 },
    ] as const;
    const pairs: BridgePair[] = [
      { a: [20, 60], b: [30, 70] },
      { a: [200, 260], b: [150, 250], strand: '-' },
      { a: [100, 200], b: [50, 150] },
    ];
    const { f, lay } = await frame({ pairs, top, bottom });
    const blockOf = (bs: readonly { range: readonly number[] }[], p: number) =>
      bs.findIndex((b) => p >= b.range[0] && p < b.range[1]);
    const seen = new Set<string>();
    for (let i = 0; i < lay.n; i++) {
      const k = lay.pair[i];
      const p = pairs[k];
      const r = lay.rib;
      // the a side's corners hold its first and its end on the block it is drawn on: each a block clipping of it
      const xs = [r[6 * i], r[6 * i + 1]];
      const candidates: number[][] = [];
      for (let b = 0; b < top.length; b++) {
        const s = Math.max(p.a[0], top[b].range[0]);
        const e = Math.min(p.a[1], top[b].range[1]);
        if (e > s) candidates.push([blockX(f, 0, b, s), blockX(f, 0, b, e)]);
      }
      expect(candidates.some((c) => Math.abs(c[0] - xs[0]) < 1e-6 && Math.abs(c[1] - xs[1]) < 1e-6)).toBe(true);
      seen.add(`${k}:${xs[0].toFixed(3)}:${r[6 * i + 3].toFixed(3)}`);
      expect(blockOf(top, p.a[0])).toBeGreaterThanOrEqual(0);
    }
    // the pair across two blocks on each side: a ribbon each place pair (2 × 2), the others one each
    expect(lay.n).toBe(1 + 1 + 4);
    expect(seen.size).toBe(lay.n);
  });
});

/** Whether a point is inside a ribbon's drawn polygon (its two cubics sampled finely, even-odd). */
function inside(r: Float64Array, o: number, x: number, y: number, line: boolean): boolean {
  if (line) return false;
  const pts: number[] = [];
  const side = (x0: number, x1: number, y0: number, y1: number) => {
    for (let s = 0; s <= 256; s++) {
      const t = s / 256;
      const S = t * t * (3 - 2 * t);
      pts.push(x0 + (x1 - x0) * S, y0 + (y1 - y0) * t);
    }
  };
  side(r[o], r[o + 3], r[o + 2], r[o + 5]);
  side(r[o + 4], r[o + 1], r[o + 5], r[o + 2]);
  let c = false;
  for (let i = 0, j = pts.length / 2 - 1; i < pts.length / 2; j = i++) {
    const xi = pts[2 * i];
    const yi = pts[2 * i + 1];
    const xj = pts[2 * j];
    const yj = pts[2 * j + 1];
    if (yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) c = !c;
  }
  return c;
}

/** The smallest distance in x from a point to any ribbon's sides at its height (a point on a boundary is skipped). */
function nearSide(lay: RibbonsLay, x: number, y: number): number {
  let d = Number.POSITIVE_INFINITY;
  for (let i = 0; i < lay.n; i++) {
    const o = 6 * i;
    const u = (y - lay.rib[o + 2]) / (lay.rib[o + 5] - lay.rib[o + 2]);
    if (u < 0 || u > 1) continue;
    const S = u * u * (3 - 2 * u);
    for (const [a, b] of [
      [o, o + 3],
      [o + 1, o + 4],
    ])
      d = Math.min(d, Math.abs(x - (lay.rib[a] + (lay.rib[b] - lay.rib[a]) * S)));
  }
  return d;
}

describe('the closed-form hit', () => {
  it('is the drawn polygon’s topmost cover: random ribbons, both directions, reversed blocks, off-window ends, stubs', async () => {
    const rng = new PCG(17);
    for (let round = 0; round < 6; round++) {
      const pairs: BridgePair[] = [];
      for (let i = 0; i < 24; i++) {
        const a0 = rng.int() % 380;
        const b0 = rng.int() % 380;
        pairs.push({
          a: [a0, a0 + 1 + (rng.int() % 20)],
          b: [b0, Math.min(400, b0 + 1 + (rng.int() % 40))],
          strand: rng.int() % 2 ? '-' : '+',
        });
      }
      const o: BridgeSpecOptions = {
        pairs,
        lengths: [400, 400],
        top: [{ range: [0, 200] }, { range: [200, 400], reversed: round % 2 === 1 }],
        // a block leaving positions out: stubs
        bottom: [{ range: [0, 320], reversed: round % 3 === 0 }],
      };
      const viewports = round > 2 ? { '/views/v/bottom': { start: 50, end: 200 } } : {};
      const { f, lay } = await frame(o, { pairs: { hitPx: 0 } }, { viewports });
      const box = f.layout.bridges!.list[0].box;
      let checked = 0;
      for (let s = 0; s < 400; s++) {
        const x = box.x + rng.float() * box.width;
        const y = box.y + 1 + rng.float() * (box.height - 2);
        if (nearSide(lay, x, y) < 0.05) continue;
        let top = -1;
        for (let i = lay.n - 1; i >= 0 && top < 0; i--)
          if (inside(lay.rib, 6 * i, x, y, (lay.flags[i] & 1) === 1)) top = i;
        const h = hitView(f.layout, x, y);
        const want = top < 0 ? undefined : `p${lay.pair[top]}`;
        expect(h?.element?.key, `round ${round} at ${x}, ${y}`).toBe(want);
        checked++;
      }
      expect(checked).toBeGreaterThan(300);
    }
  });
});
