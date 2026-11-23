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
import { SeqRuntime } from '@molstar/seqx-core/data';
import { bridgeSpec, randomBridgePairs } from './bridge-fixtures.js';
import { AlnLength, randomAlignments, regionBridgeSpec } from './bridge-region-fixtures.js';
import { settle } from './helpers.js';
import { reusedHover } from './scenarios.js';

const N = 20_000;
const input: FrameInput = { size: { width: 1440, height: 900 }, dpr: 2, settled: true };

/** A hit's bytes: the two numbers the test's own call boxes, nothing past them (no hit object). */
const HitBound = 32;

/**
 * A pan's layout of both sides and the bridge over a region source (its reads kept, `R` read again and fitted to the
 * layer's bounds, the ribbons laid out), bytes a layout: measured 486 kB (440 before `R`'s fit check), held, not a
 * target (the sides' layouts and the reads' regions are most of it).
 */
const RelayoutBound = 500_000;

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

describe('a bridge laid out again', () => {
  it('a pan of a side over a region source relays the bridge each frame, a bounded size a layout', async ({ task }) => {
    const s = await settle(regionBridgeSpec(randomAlignments(500, 3)), [SeqRuntime], {
      pairs: { maxArcs: 1024, binPx: 40 },
    });
    // the top zoomed in (the near side, its `R` read), the bottom at the fit: two windows inside the same tiles
    const at = (x: number) => ({
      ...input,
      viewports: {
        '/views/v/top': { start: 400_000 + x, end: 450_000 + x },
        '/views/v/bottom': { start: 0, end: AlnLength },
      },
    });
    const frames = [at(0), at(37)];
    for (const f of frames) await s.frame(f);
    const step = (i: number) => void s.renderer.frame(frames[i & 1]);
    const r = quickCase(task, 'bridge relayout', step);
    console.info(`bridge relayout: ${r.r.bytesPerStep.toFixed(0)} bytes a layout (${r.r.note})`);
    expect(r.r.bytesPerStep, r.r.note).toBeLessThanOrEqual(RelayoutBound);
  });
});
