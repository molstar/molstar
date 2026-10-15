/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST6's heap guards (seqspec:sources#the-viewers-side; the memory rule): a pass with prefetch on allocates nothing once
// what it asks is held: the lead past the margin as a pan alternates sides, the next level as a zoom alternates in and
// out, and the snapshots' focus regions (two hints, held). A pass that allocates an object a step fails the same guard
// (the positive control). The cases are streams-prefetch-cases.ts's, shared with the quick twin
// (streams-prefetch.quick.ts); each loads the viewer afresh (`vi.resetModules`) and is warmed and settled by count
// (`frameScaling`).

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { MOTIONS, motionName, prefetchCase } from './streams-prefetch-cases.js';

const N = 1_000;

describe('a pass with prefetch on allocates nothing once what it asks is held', () => {
  for (const motion of MOTIONS)
    it(motionName(motion), async () => {
      const g = await prefetchCase(motion);
      // Not vacuous: prefetches were asked (the lead or the zoom's level, and the snapshots' regions).
      expect(g.prefetched()).toBeGreaterThan(motion === 'pan' ? 3 : 1);
      const asked = g.asked();
      const r = frameScaling(g.pass, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(g.asked()).toBe(asked);
      g.dispose();
    }, 300_000);

  it('the same pass allocating an object a step is caught (positive control)', async () => {
    const g = await prefetchCase('pan');
    expect(controlCaught(() => frameScaling(leaking(g.pass), N))).toBe(true);
    g.dispose();
  }, 300_000);
});
