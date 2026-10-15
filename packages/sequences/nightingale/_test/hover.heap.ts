/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Heap guards for NG1's hot paths, on RB0b's stepped `heapScaling` (a step is one hover change, 50,000 and 200,000 of
// them): a host's loci through `Extents.addLoci`, `project` per block into axis runs (sorted and merged), with (a) a
// direct block, (b) a residue onto a gene block through the codon route, (c) a link to a second section's manager, (d)
// the `containing` closure first; (e) the pointer's `hoverAt` across positions; (f) `rangesOf` from a feature's
// fragments and axis runs as block runs. The cases are hover-cases.ts's, shared with the quick twin (hover.quick.ts).
// Counted, not timed: a repeated highlight builds no string, and a unit's lanes probe O(n log n). The positive control
// is the template's: (a)'s step with a leak, caught by the same measurement (NG-11).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, heapScaling, isFlat, leaking } from '../../core/_test/model/heap.js';
import { type Case, codon, containing, direct, fragments, linked, N, pointer, repeated, twice } from './hover-cases.js';

describe('the highlight path allocates nothing per change', () => {
  const run = async (c: Case) => {
    const result = heapScaling(N, c.step);
    expect(c.seen()).toBeGreaterThan(0);
    expect(isFlat(result), describeScaling(result, N)).toBe(true);
  };

  it('(a) a nucleotide on its direct block', async () => run(await direct()));

  it('(a) with a leak is caught (the positive control)', async () => {
    const { step } = await direct();
    expect(controlCaught(() => heapScaling(N, leaking(step)))).toBe(true);
  });

  it('(b) a residue on the gene block through the codon route', async () => run(await codon()));

  it('(c) a residue to every manager: the protein, the gene and the chromosome across the link', async () => {
    await run(await linked());
  });

  it('(d) a nucleotide with the containing closure first, then projected', async () => run(await containing()));

  it('a source used twice, at the end of its other section', async () => run(await twice()));

  it("(e) the pointer's hoverAt across positions, deduped before it locates", async () => run(await pointer()));

  it("(f) a feature's fragments as axis runs, and axis runs as block runs", async () => run(await fragments()));

  it('a repeated highlight builds no string and dispatches nothing', async () => {
    const sent = await repeated();
    expect(sent).toHaveLength(1);
    expect(typeof sent[0]).toBe('string');
  });
});
