/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The rings layout's allocation targets (performance.md, Memory: draw and hover 0 bytes warm through the warp): a warm
// draw of the plasmid-sized ring and its detail (every pass, the ring through the polar warp), a warm hover over the
// ring (the emphasis redraw through the warp, the window lit), the selection and cursor walking on it. By hand, never
// in CI: `pnpm -C packages/seq/core test:perf`.

import { rangesLocation } from '@molstar/seqx-core/query';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { drawView, type FrameInput, type Hit, hitView, type RingPlacement } from '../index.js';
import { settle } from './helpers.js';
import { plasmidRingSpec } from './rings-fixtures.js';
import { reusedHover } from './scenarios.js';

/**
 * What a warm ring frame still allocates: one number boxed a frame where the window's grips go through the warp (16 B,
 * measured 2026-11-19; TODO: found and removed, then this goes).
 */
const WindowBoxed = 16;

const Input: FrameInput = {
  size: { width: 1200, height: 1200 },
  dpr: 2,
  viewports: { '/views/v/ring': { start: 2500, end: 2800 } },
};

async function ring() {
  const s = await settle(plasmidRingSpec());
  const f = await s.frame(Input);
  return { s, f, p: f.layout.sections[0].placement as RingPlacement };
}

describe('a ring draws and hovers with nothing allocated once warm', () => {
  it('the plasmid ring: draw', async ({ task }) => {
    const { f } = await ring();
    const { r } = quickCase(task, 'ring draw', () => drawView(f.layout, f.lists), {});
    console.info(`ring draw: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold') + WindowBoxed);
  });

  it('the plasmid ring: hover', async ({ task }) => {
    const { s, f, p } = await ring();
    const hits: Hit[] = [];
    for (let u = 0; u < p.total; u += 11) {
      p.warp.map(p.frameX(u + 0.5), p.strips.bandY[0] + 3);
      const h = hitView(f.layout, p.warp.px, p.warp.py);
      if (h) hits.push(h);
    }
    const hover = reusedHover(hits);
    const at = { ...Input, emphasis: hover.emphasis };
    const { r } = quickCase(
      task,
      'ring hover',
      (i) => {
        hover.set(i);
        s.renderer.frame(at);
      },
      {},
    );
    console.info(`ring hover: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold') + WindowBoxed);
  });

  it('the plasmid ring: selection and cursor', async ({ task }) => {
    const { s, p } = await ring();
    const system = p.strips.axes[0].system;
    const selections = [0, 1, 2, 3].map((k) => rangesLocation(system, [2600 - 100 * k, 2650 + 10 * k]));
    const cursor = { system, position: 0 };
    const emphasis: { selection?: (typeof selections)[number]; cursor?: typeof cursor } = { cursor };
    const at = { ...Input, emphasis };
    const { r } = quickCase(
      task,
      'ring selection',
      (i) => {
        emphasis.selection = selections[i & 3];
        cursor.position = i % 2686;
        s.renderer.frame(at);
      },
      {},
    );
    console.info(`ring selection and cursor: ${r.bytesPerStep.toFixed(0)} bytes a frame (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold') + WindowBoxed);
  });
});
