/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rings layout's targets: the plasmid-sized ring (four rings, about 30 labels) and
// its linear detail, a pan of the detail (the ring's window moved: the detail laid out again, the ring drawn from its
// kept lays, its labels not placed again) under 8 ms a frame; a hover over the ring (the placement's `at` through the
// warp's inverse, the layer's hit, the emphasis redraw) under 2 ms. Each bound is the heavy lane's (`itWithin`, best of
// three); the quick lane runs them unbounded.

import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import { type FrameInput, hitView, type RingPlacement } from '../index.js';
import { settle } from './helpers.js';
import { plasmidRingSpec } from './rings-fixtures.js';

const Size = { width: 1200, height: 1200 };
const Frames = 40;

describe('the rings layout’s time', () => {
  itWithin('a pan of the ring’s linear detail: the ring drawn from its kept lays', async (within) => {
    const s = await settle(plasmidRingSpec());
    const inputs: FrameInput[] = [];
    for (let k = 0; k < 8; k++) {
      const start = 2400 + k * 37;
      inputs.push({ size: Size, dpr: 2, viewports: { '/views/v/ring': { start, end: start + 300 } } });
    }
    for (const i of inputs) await s.frame(i);
    const ring = s.renderer.frame(inputs[0]).layout.sections[0];
    expect((ring.placement as RingPlacement).labels.m).toBeGreaterThanOrEqual(25);
    const time = () => {
      const t: number[] = [];
      for (let k = 0; k < Frames; k++) {
        const t0 = performance.now();
        const f = s.renderer.frame(inputs[(k + 1) % inputs.length]);
        t.push(performance.now() - t0);
        if (!f.rebuilt || f.layout.sections[0] !== ring) throw new Error('a pan lays the detail out, not the ring');
      }
      t.sort((a, b) => a - b);
      return t[t.length >> 1];
    };
    time();
    const median = within(8 * Frames, time);
    console.info(`ring, a pan of its detail: median ${median.toFixed(3)} ms`);
  });

  itWithin('a hover over the ring: the placement through the warp, the hit, the emphasis redraw', async (within) => {
    const s = await settle(plasmidRingSpec());
    const input: FrameInput = { size: Size, dpr: 2, viewports: { '/views/v/ring': { start: 100, end: 400 } } };
    const f = await s.frame(input);
    const layout = f.layout;
    const p = layout.sections[0].placement as RingPlacement;
    const points: { x: number; y: number }[] = [];
    for (let u = 0; u < p.total; u += 23) {
      p.warp.map(p.frameX(u + 0.5), p.strips.bandY[0] + 3);
      points.push({ x: p.warp.px, y: p.warp.py });
    }
    const emphasis = points.map((q) => ({ ...input, emphasis: { hover: hitView(layout, q.x, q.y) } }));
    const hover = () => {
      const t: number[] = [];
      for (let k = 0; k < 200; k++) {
        const q = points[k % points.length];
        const t0 = performance.now();
        hitView(layout, q.x, q.y);
        s.renderer.frame(emphasis[k % emphasis.length]);
        t.push(performance.now() - t0);
      }
      t.sort((a, b) => a - b);
      return t[t.length >> 1];
    };
    hover();
    const median = within(2 * 200, hover);
    console.info(`ring hover: median ${median.toFixed(3)} ms`);
    expect(points.length).toBeGreaterThan(100);
  });
});
