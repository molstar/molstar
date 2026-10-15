/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of wrapped-tracks.heap.ts (`test:perf`): the same case (wrapped-tracks-cases.ts), the text view's tracks
// drawn in full each step as the rows scroll down and back, with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { wrappedTracks } from './wrapped-tracks-cases.js';

/** A frame redraws every line's tracks (about 1 ms on an M4): windows of 250 steps, tier-up done by 1,000. */
const DRAW: QuickOptions = { warmSteps: 1_000, steps: 250, windows: 3, minWarmSteps: 1_500 };

describe('tracks drawn in full in the text, quick', () => {
  it('a frame of a values line, an area, glyphs, labeled ranges and secondary structure, axes and labels on', async ({
    task,
  }) => {
    const d = await wrappedTracks('all');
    quickCase(task, 'wrapped-tracks: all', d.step, DRAW);
    const { paths, glyphs, texts } = d.state();
    expect(paths).toBeGreaterThan(500);
    expect(glyphs).toBeGreaterThan(0);
    expect(texts).toBeGreaterThan(0);
  });
});
