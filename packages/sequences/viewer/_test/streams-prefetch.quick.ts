/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-prefetch.heap.ts (`test:perf`): a pass with prefetch on as the view pans to and fro and as
// it zooms in and out (streams-prefetch-cases.ts, the guard's cases), each with its positive control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { MOTIONS, motionName, prefetchCase } from './streams-prefetch-cases.js';

describe('a pass with prefetch on, quick', () => {
  for (const motion of MOTIONS)
    it(motionName(motion), async ({ task }) => {
      const g = await prefetchCase(motion);
      expect(g.prefetched()).toBeGreaterThan(motion === 'pan' ? 3 : 1);
      const asked = g.asked();
      quickCase(task, `streams-prefetch: ${motion}`, g.pass);
      expect(g.asked()).toBe(asked);
      g.dispose();
    });
});
