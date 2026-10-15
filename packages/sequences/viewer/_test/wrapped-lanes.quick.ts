/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of wrapped-lanes.heap.ts (`test:perf`): the same cases (wrapped-lanes-cases.ts), the draw first, each
// with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { DRAW, wrappedLanes } from './wrapped-lanes-cases.js';

describe("the text view's lanes, quick", () => {
  it('a frame of the lanes, scrolled', async ({ task }) => {
    const d = await wrappedLanes('draw');
    quickCase(task, 'wrapped-lanes: draw', d.step, DRAW);
    expect(d.state().rects).toBeGreaterThan(500);
  });

  it('the pointer over the strips and the underlines', async ({ task }) => {
    const d = await wrappedLanes('hover');
    quickCase(task, 'wrapped-lanes: hover', d.step);
    expect(d.state().hovered).toBeGreaterThan(0);
  });
});
