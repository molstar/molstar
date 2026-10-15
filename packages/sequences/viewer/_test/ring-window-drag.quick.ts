/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of ring-window-drag.heap.ts (`test:perf`): the ring window's drag, each move setting the linked linear
// section, on a viewer loaded afresh (ring-window-heap-run.ts, the guards' cases), with its control.

import { expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { RING_WINDOW, setup } from './ring-window-heap-run.js';

it("CG3's ring window drag, each move setting the linear section, quick", async ({ task }) => {
  const s = await setup('drag', 600);
  quickCase(task, 'ring-window: drag', s.step, RING_WINDOW);
  expect(s.owner()).toBe('ringWindow');
});
