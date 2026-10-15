/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of ideogram.heap.ts (`test:perf`): the same case (ideogram-cases.ts), the strip and the window over it
// redrawn a step, with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { ideogramFrame } from './ideogram-cases.js';

describe('the ideogram navigator, quick', () => {
  it('its strip and the window over it, three chromosomes', async ({ task }) => {
    const g = await ideogramFrame();
    g.frame();
    expect(g.backend.rectItems.n).toBeGreaterThan(80);
    expect(g.backend.pathItems.nv).toBeGreaterThanOrEqual(24);
    quickCase(task, 'ideogram: strip and window', g.frame);
    g.dispose();
  });
});
