/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick twin of streams-paths.heap.ts (`test:perf`): a protein's section showing its genome's streamed tracks
// through the codon map (streams-paths-cases.ts, the guard's cases): a pan's pass and composer frame over 1,000 residues
// and at 100 residues across a split codon, the pointer over the genes and variants there, the codon label's reads and
// a pass with letters wanted, each with its positive control.

import { describe, expect, it } from 'vitest';
import { type QuickOptions, quickCase } from './heap.js';
import { codonReads, pointerOverPath, PROTEIN_VIEWS, proteinCase } from './streams-paths-cases.js';

/**
 * A frame's counts: in a cold process a streamed frame sits in its middle tier for up to about 10,000 steps, reading
 * a steady 16 B a step there; the warm-up runs on until it reads quiet (alloc.ts), at most 40,000 steps.
 */
const FRAME: QuickOptions = { warmSteps: 2_000, steps: 500, windows: 3, maxWarmSteps: 40_000 };

describe("a genome's streamed tracks on a protein's axis, quick", () => {
  for (const [label, start, stop] of PROTEIN_VIEWS)
    it(`a pass and a frame, ${label}`, async ({ task }) => {
      const g = await proteinCase(start, stop);
      g.frame();
      expect(g.backend.rectItems.n).toBeGreaterThan(10);
      const asked = g.asked();
      expect(asked).toBeGreaterThan(0);
      quickCase(task, `streams-paths: ${label}`, g.frame, FRAME);
      expect(g.asked()).toBe(asked);
      g.dispose();
    });

  it('the pointer over the genes and variants (a hit through the path)', async ({ task }) => {
    const g = await pointerOverPath();
    expect(g.parts).toBeGreaterThan(10);
    quickCase(task, 'streams-paths: pointer', g.step);
    g.dispose();
  });

  it("the codon label's reads and a pass with letters wanted", async ({ task }) => {
    const g = await codonReads();
    expect(g.held).toBe(300);
    quickCase(task, 'streams-paths: codon reads', g.read);
    quickCase(task, 'streams-paths: codon reads, letters wanted', g.far);
    g.dispose();
  });
});
