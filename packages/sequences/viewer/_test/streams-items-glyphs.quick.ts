/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-items-glyphs.heap.ts (`test:perf`): a pan's pass and composer frame over the streamed
// glyphs at 5 kb (stacked and as points) and over the placeholder bands with every region pending
// (streams-items-glyphs-cases.ts, the guard's cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { glyphsCase, pendingBands } from './streams-items-glyphs-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe('the streamed glyphs and placeholders, quick', () => {
  it('a pass and a frame, 5 kb', async ({ task }) => {
    const g = await glyphsCase('glyphs');
    g.frame();
    expect(g.backend.glyphItems.n).toBeGreaterThan(0);
    const asked = g.asked();
    quickCase(task, 'streams-items-glyphs: 5 kb', g.frame, FRAME);
    expect(g.asked()).toBe(asked);
    g.dispose();
  });

  it('a pass and a frame, 5 kb, the variants as points', async ({ task }) => {
    const g = await glyphsCase('points');
    g.frame();
    expect(g.backend.glyphItems.n).toBeGreaterThan(0);
    const asked = g.asked();
    quickCase(task, 'streams-items-glyphs: 5 kb points', g.frame, FRAME);
    expect(g.asked()).toBe(asked);
    g.dispose();
  });

  it('a frame with every region pending (placeholder bands)', async ({ task }) => {
    const g = await glyphsCase('pending');
    g.frame();
    expect(await pendingBands(g)).toBeGreaterThan(2);
    quickCase(task, 'streams-items-glyphs: pending', g.frame, FRAME);
    g.dispose();
  });
});
