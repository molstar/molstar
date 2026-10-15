/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of row-hover.heap.ts (`test:perf`): the hover sweep down the tracks, an overlay frame a move, the band
// drawn (row-hover-cases.ts, the guard's case), with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { rowHover } from './row-hover-cases.js';

it('a hover sweep down the tracks, the band drawn', async ({ task }) => {
  const { frame, hovers, drawn } = await rowHover();
  quickCase(task, 'row-hover: sweep down the tracks', frame);
  expect(hovers()).toBeGreaterThan(1_000);
  expect(drawn()).toBeGreaterThan(1_000);
});
