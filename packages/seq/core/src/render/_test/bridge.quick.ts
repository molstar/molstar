/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A bridge's allocation targets (_spec/README.md "Render", Bridges): a warm draw of 4,096 ribbons and a warm hover (the
// hovered ribbon at full opacity, the ribbons at a hovered position lit) allocate nothing; a ribbon's hit allocates no
// hit object. By hand, never in CI: `pnpm -C packages/seq/core test:perf`.

import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, type FrameInput, type Hit, hitView } from '../index.js';
import type { RibbonsLay } from '../pairs-ribbons.js';
import { bridgeSpec, randomBridgePairs } from './bridge-fixtures.js';
import { settle } from './helpers.js';
import { reusedHover } from './scenarios.js';

const N = 20_000;
const input: FrameInput = { size: { width: 1440, height: 900 }, dpr: 2, settled: true };

/** A hit's bytes: the two numbers the test's own call boxes, nothing past them (no hit object). */
const HitBound = 32;

describe('a bridge draws and hovers with nothing allocated once warm', () => {
  it('4,096 ribbons: draw, a hover lighting a ribbon, a hit', async ({ task }) => {
    const s = await settle(bridgeSpec({ pairs: randomBridgePairs(4096, 23), lengths: [N, N] }));
    const f = await s.frame(input);
    const d = quickCase(task, 'bridge draw', () => drawView(f.layout, f.lists), {});
    console.info(`bridge draw: ${d.r.bytesPerStep.toFixed(0)} bytes a frame (${d.r.note})`);
    expect(d.r.bytesPerStep, d.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const lay = f.layout.bridges!.list[0].lays[0].lay as RibbonsLay;
    const points: [number, number][] = [];
    for (let i = 0; i < lay.n; i += 97) {
      const o = 6 * i;
      points.push([
        (lay.rib[o] + lay.rib[o + 1] + lay.rib[o + 3] + lay.rib[o + 4]) / 4,
        (lay.rib[o + 2] + lay.rib[o + 5]) / 2,
      ]);
    }
    const hits: Hit[] = [];
    for (const p of points) {
      const h = hitView(f.layout, p[0], p[1]);
      if (h) hits.push({ ...h });
    }
    expect(hits.length).toBeGreaterThan(20);
    const hover = reusedHover(hits);
    const hoverAt = { ...input, emphasis: hover.emphasis };
    const h = quickCase(
      task,
      'bridge hover',
      (i) => {
        hover.set(i);
        s.renderer.frame(hoverAt);
      },
      {},
    );
    console.info(`bridge hover: ${h.r.bytesPerStep.toFixed(0)} bytes a frame (${h.r.note})`);
    expect(h.r.bytesPerStep, h.r.note).toBeLessThanOrEqual(inject('perfThreshold'));
    const flat = Float64Array.from(points.flat());
    const m = points.length;
    const hit = quickCase(
      task,
      'bridge hit',
      (i) => void hitView(f.layout, flat[2 * (i % m)], flat[2 * (i % m) + 1]),
      {},
    );
    console.info(`bridge hit: ${hit.r.bytesPerStep.toFixed(0)} bytes a hit (${hit.r.note})`);
    expect(hit.r.bytesPerStep, hit.r.note).toBeLessThanOrEqual(HitBound);
  });
});
