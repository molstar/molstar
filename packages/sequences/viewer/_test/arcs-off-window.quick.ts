/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of arcs-off-window.heap.ts (`test:perf`): the same case (arcs-off-window-cases.ts), a pan, a frame, a
// hover on an end and its overlay a step, with its positive control.

import { describe, expect, it } from 'vitest';
import { arcsOffWindow } from './arcs-off-window-cases.js';
import { quickCase } from './heap.js';

describe('frames of arcs whose partners are beyond the window, quick', () => {
  it('pans, a hover on an end and its overlay', async ({ task }) => {
    const { frame, glyphs } = await arcsOffWindow();
    expect(glyphs.length).toBeGreaterThan(5);
    quickCase(task, 'arcs-off-window: pan and hover', frame);
  });
});
