/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of labels-place.heap.ts (`test:perf`): `placeLabel` at a panned view each step, the guard's step
// (labels-heap-fixture.ts `placing`), with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from '../../core/_test/model/quick.js';
import { placing } from './labels-heap-fixture.js';

it('placeLabel at a panned view each step', ({ task }) => {
  const { step, modes } = placing();
  quickCase(task, 'labels-place: placeLabel at a panned view', step);
  expect([modes.full > 0, modes.cut > 0, modes.none > 0]).toEqual([true, true, true]);
});
