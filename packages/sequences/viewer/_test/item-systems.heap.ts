/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST7's item coordinate systems' heap guard (the memory rule): once a transcript's coordinate system is held, a pan's
// pass and frame over its own residues (P1's domains through a document mapping onto it, the genome's variants and
// phyloP through its map) and over the genome (P1's domains along its exons) allocate nothing. A frame allocating an
// object a step fails the same guard (the positive control). Each case loads the viewer afresh (`vi.resetModules`),
// is warmed by count (`frameScaling`) and shares its cases with the quick twin (item-systems-cases.ts,
// item-systems.quick.ts).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { SYSTEMS, transcriptFrame } from './item-systems-cases.js';

const N = 1_000;

describe("an item coordinate system's tracks allocate nothing per frame once held", () => {
  for (const [label, section, block, start, stop] of SYSTEMS)
    it(`a pass and a frame: ${label}`, async () => {
      const g = await transcriptFrame(section, block, start, stop);
      g.frame();
      // Not vacuous: features drawn (rects); on T1's residues, the genome's data read through its map.
      expect(g.backend.rectItems.n).toBeGreaterThan(2);
      if (section === 'tx.item') expect(g.throughMap()).toEqual([true, true]);
      const r = frameScaling(g.frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      g.dispose();
    }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await transcriptFrame('tx.item', 0, 450, 650);
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
