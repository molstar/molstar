/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of text-rows.heap.ts (`test:perf`): an ESPript-like figure (text-rows-cases.ts) drawn as its text
// scrolls, hovered over its rows and a frame's title, and in tracks panned and row-zoomed, each with its positive
// control (seqspec:viewer#rows-in-text, seqspec:viewer#frames).

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { textRows } from './text-rows-cases.js';

/** A frame redraws every line in view: windows of 250 steps, tier-up done by 1,000. */
const DRAW: QuickOptions = { warmSteps: 1_000, steps: 250, windows: 3, minWarmSteps: 1_500 };

describe('an alignment as text with frames, quick', () => {
  it('a frame of the rows as text, a box, frames, cylinders and markers, scrolled', async ({ task }) => {
    const d = await textRows('draw');
    d.resetCounts();
    quickCase(task, 'text-rows: draw', d.step, DRAW);
    expect(d.state().paths).toBeGreaterThan(100);
    expect(d.state().texts).toBeGreaterThan(1000);
    // The frames' outlines were drawn in the steps measured.
    expect(d.state().frames).toBeGreaterThan(100);
  });

  it("the pointer over the rows' cells, a gap and a frame's title", async ({ task }) => {
    const d = await textRows('hover');
    quickCase(task, 'text-rows: hover', d.step);
    expect(d.state().hovered).toBeGreaterThan(0);
  });

  it('frames in tracks: a pan, the strip and its frames', async ({ task }) => {
    const d = await textRows('tracks');
    quickCase(task, 'frames: tracks draw', d.step, DRAW);
    expect(d.state().paths).toBeGreaterThan(0);
  });

  it("frames in tracks: a row zoom (a rescale), the frames over the rows' new heights", async ({ task }) => {
    const d = await textRows('zoom');
    const before = d.stats();
    quickCase(task, 'frames: row zoom', d.step, DRAW);
    expect(d.state().paths).toBeGreaterThan(0);
    // Each step a rescale, no build.
    expect(d.stats().rescales - before.rescales).toBeGreaterThan(1000);
    expect(d.stats().builds).toBe(before.builds);
  });
});
