/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of overlay.heap.ts (`test:perf`): the pointer moving over the scale and over the alignment's rows, an
// overlay frame a move (the guard's cases, overlay-cases.ts), each with its control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { hovering, KINDS } from './overlay-cases.js';

describe('overlay frames, quick', () => {
  it.for(KINDS)('hover frames over %s with a selection and a highlight', async (kind, { task }) => {
    const { frame, hovers } = await hovering(kind);
    quickCase(task, `overlay: hover frames over ${kind}`, frame);
    expect(hovers()).toBeGreaterThan(1_000);
  });
});
