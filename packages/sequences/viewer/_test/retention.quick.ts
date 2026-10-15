/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of retention.heap.ts (`test:perf`). The guard's check (weak references to a disposed viewer's set and
// view, collected a task apart) is once-only and async, not a step; what it exercises repeatably is the projection a
// viewer is used through: an item's extents closed over the domain, projected into the hover marks, closed again
// (retention-cases.ts `projecting`). That is the step measured here, with its positive control.

import { expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { projecting } from './retention-cases.js';

it('the projection a viewer is used through', async ({ task }) => {
  const { step, projected, dispose } = await projecting();
  quickCase(task, 'retention: the projection a viewer is used through', step);
  expect(projected()).toBeGreaterThan(0);
  dispose();
});
