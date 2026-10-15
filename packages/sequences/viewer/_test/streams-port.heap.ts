/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST1's heap guards (seqspec:sources#the-viewers-side; the memory rule): over the streamed hg38-sized genome with its
// regions held, a pan's pass, and a pan's pass and composer frame (the values drawn from the regions,
// `StreamView.values` per device column), allocate nothing, at the whole view and at 10 kb; a frame that allocates
// fails the same guard (the positive control). The cases are streams-port-cases.ts's, shared with the quick twin
// (streams-port.quick.ts); each loads the viewer and its frame afresh (`vi.resetModules`), as rn3.heap.ts does, and is
// warmed and settled by count (`frameScaling`): at 10 kb with the frame drawn the climb to TurboFan runs thousands of
// frames at 8 to 122 bytes a frame.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { PORT_CASES, portFrame } from './streams-port-cases.js';

const N = 1_000;

describe("the streamed genome's pans allocate nothing per frame once its regions are held", () => {
  for (const { name, at, draw, dynamic } of PORT_CASES)
    it(name, async () => {
      const g = await portFrame(at, draw, dynamic);
      if (dynamic) expect(g.manager.nodes.entries.size).toBeGreaterThan(at === 'whole' ? -1 : 0);
      const asked = g.requests();
      const r = frameScaling(g.frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(g.requests()).toBe(asked);
      if (!dynamic) expect(g.manager.stats.current.hits).toBeGreaterThan(0);
    }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await portFrame('10 kb', true);
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
  }, 300_000);
});
