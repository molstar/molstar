/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The pairs' and grids' time targets (the design's "Performance and size"): a 2,000 × 2,000 PAE's pan frame under 8 ms
// once its tiles are built, a frame while moving building within `grid.buildBudget`; 4,096 arcs a frame under 8 ms and a
// hover over nested arcs under 2 ms. The bounds are the heavy lane's (`itWithin`, best of three); the quick lane runs
// each unbounded, and the budget's count is exact in both.

import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect, it } from 'vitest';
import { type FrameInput, type GridLay, hitView, type PairsLay, type Renderer } from '../index.js';
import { ImageTile } from '../images.js';
import { settle } from './helpers.js';
import { arcsSpec, gridSpec, randomPairs } from './pairs-fixtures.js';

const Size = { width: 1440, height: 1600 };
const Leader = '/views/v/s';

const at = (start: number, end: number, settled = true): FrameInput => ({
  size: Size,
  dpr: 2,
  settled,
  viewports: { [Leader]: { start, end } },
});

/** A pan frame: the next of the inputs (every frame lays out again), timed by the bound around it. */
function panner(renderer: Renderer, inputs: readonly FrameInput[]): () => void {
  let k = 0;
  return () => {
    const f = renderer.frame(inputs[++k % inputs.length]);
    if (!f.rebuilt) throw new Error('a pan frame lays out again');
  };
}

describe('a dot plot’s frame, 2,000 × 2,000', () => {
  itWithin('a pan frame once its tiles are built', async (within) => {
    const s = await settle(gridSpec({ n: 2000 }));
    const inputs = [0, 1, 2, 3].map((k) => at(200 + 3 * k, 1400 + 3 * k));
    for (const i of inputs) await s.frame(i);
    const lay = (await s.frame(inputs[0])).layout.sections[1].lays[0][0].lay as GridLay;
    expect(lay.images.length).toBeGreaterThan(0);
    const pan = panner(s.renderer, inputs);
    for (let k = 0; k < 8; k++) pan();
    within(8, pan);
  });

  it('a frame while moving builds at most a tile’s samples (`grid.buildBudget`), a coarser tile meanwhile', async () => {
    const s = await settle(gridSpec({ n: 2000 }));
    const images = s.renderer.images;
    for (const [start, end] of [
      [0, 2000],
      [300, 1100],
      [500, 700],
    ]) {
      const made = images.made;
      const size = images.size;
      const f = s.renderer.frame(at(start, end, false));
      const lay = f.layout.sections[1].lays[0][0].lay as GridLay;
      expect(lay.images.length).toBeGreaterThan(0);
      // a full tile and the coarse tiles drawn meanwhile (each a sample's worth at its level)
      expect(images.size - size).toBeLessThanOrEqual(1 + 6);
      expect(images.made - made).toBeLessThanOrEqual(1 + 6);
      for (const img of lay.images) expect(img.width * img.height).toBeLessThanOrEqual(ImageTile * ImageTile);
    }
  });
});

describe('4,096 arcs', () => {
  const pairs = randomPairs(20_000, 4096, 17);

  itWithin('a frame of 4,096 arcs', async (within) => {
    const s = await settle(arcsSpec({ n: 20_000, pairs, height: 120 }));
    const inputs = [0, 1, 2, 3].map((k) => at(k, 20_000 - 3 + k));
    for (const i of inputs) await s.frame(i);
    const lay = (await s.frame(inputs[0])).layout.sections[0].sectionLays[0].lay as PairsLay;
    expect(lay.n).toBeGreaterThan(4000);
    const pan = panner(s.renderer, inputs);
    for (let k = 0; k < 8; k++) pan();
    within(8, pan);
  });

  itWithin('a hover over nested arcs', async (within) => {
    const s = await settle(arcsSpec({ n: 20_000, pairs, height: 120 }));
    const f = await s.frame(at(0, 20_000));
    const lay = f.layout.sections[0].sectionLays[0].lay as PairsLay;
    // the points: apexes and the band's foot, where arcs crowd
    const points: [number, number][] = [];
    for (let i = 0; i < lay.n; i += 37) {
      const o = 8 * i;
      points.push([(lay.curve[o] + lay.curve[o + 6]) / 2, (lay.curve[o + 3] + lay.curve[o + 1]) / 2]);
      points.push([lay.curve[o] + 1, lay.curve[o + 1] - 1]);
    }
    // the slowest point bounded: each hover once warm
    for (const p of points) hitView(f.layout, p[0], p[1]);
    for (const p of points) within(2, () => hitView(f.layout, p[0], p[1]));
  });
});
