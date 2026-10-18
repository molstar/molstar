/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7's heap guard (seqspec:sources#requests; the memory rule): a protein's section showing its genome's streamed
// genes, variants and phyloP through its codon map, regions held. A pan's pass (the view's residues taken back through
// the path to its exons' runs, the regions asked) and its frame (items read back on the residues, values over codons,
// a split codon's two parts, placeholders mapped through the path) allocate nothing, over 1,000 residues (three exons)
// and at 100 residues across a split codon; so does the pointer over the genes and variants there (`hitStreamed`
// through the path, every region over a split codon tried), and so do the codon label's reads (a residue's positions,
// `Paths.fromAt`, and their letters held, `lettersAt`) and a pass with letters wanted (the aside layer in its order). A
// frame allocating an object a step fails the same guard (the positive control). The cases are streams-paths-cases.ts's,
// shared with the quick twin (streams-paths.quick.ts); each loads the viewer afresh (`vi.resetModules`) and is warmed
// and settled by count (`frameScaling`): colored items tier up slowly.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { codonReads, pointerOverPath, PROTEIN_VIEWS, proteinCase } from './streams-paths-cases.js';

const N = 1_000;

describe("a genome's streamed tracks on a protein's axis allocate nothing per frame once held", () => {
  for (const [label, start, stop] of PROTEIN_VIEWS)
    it(`a pass and a frame, ${label}`, async () => {
      const g = await proteinCase(start, stop);
      g.frame();
      // Not vacuous: values and items drawn (rects).
      expect(g.backend.rectItems.n).toBeGreaterThan(10);
      const asked = g.asked();
      expect(asked).toBeGreaterThan(0);
      const r = frameScaling(g.frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(g.asked()).toBe(asked);
      g.dispose();
    }, 300_000);

  it('the pointer over the genes and variants (a hit through the path)', async () => {
    const g = await pointerOverPath();
    // Not vacuous: items hit.
    expect(g.parts).toBeGreaterThan(10);
    const r = frameScaling(g.step, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    g.dispose();
  }, 300_000);

  it("the codon label's reads: a residue's positions through the path and their letters, held", async () => {
    const g = await codonReads();
    // Not vacuous: every base of the 100 residues held.
    expect(g.held).toBe(300);
    const r = frameScaling(g.read, N);
    expect(isFlat(r), describeScaling(r, N)).toBe(true);
    // A pass with the aside layer's letters wanted: pushed into the pass's order, asked, the order cut back.
    const far = frameScaling(g.far, N);
    expect(isFlat(far), describeScaling(far, N)).toBe(true);
    g.dispose();
  }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await proteinCase(450, 550);
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
