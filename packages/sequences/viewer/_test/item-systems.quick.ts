/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of item-systems.heap.ts (`test:perf`): a pass and a frame over the transcript's own residues and over
// the genome along its exons, the guard's cases (item-systems-cases.ts), each with its control.

import { describe, expect, it } from 'vitest';
import { quickCase } from './heap.js';
import { SYSTEMS, transcriptFrame } from './item-systems-cases.js';

describe("an item coordinate system's tracks, quick", () => {
  for (const [label, section, block, start, stop] of SYSTEMS)
    it(`a pass and a frame: ${label}`, async ({ task }) => {
      const g = await transcriptFrame(section, block, start, stop);
      g.frame();
      expect(g.backend.rectItems.n).toBeGreaterThan(2);
      quickCase(task, `item-systems: ${section}`, g.frame);
      g.dispose();
    });
});
