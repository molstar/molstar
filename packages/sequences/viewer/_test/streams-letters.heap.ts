/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// ST3a's heap guards (seqspec:sources#the-viewers-side; the memory rule): over the streamed genome with a sequence
// track and its regions held, a pan's pass and composer frame allocate nothing: at 10 kb (the letters' cells thinner
// than a device px, read raw per device column from `StreamView.letters`) and at 100 positions (the letters drawn); a
// frame allocating an object a step fails the same guard (the positive control). The cases are
// streams-letters-cases.ts's, shared with the quick twin (streams-letters.quick.ts); each loads the viewer and its
// frame afresh (`vi.resetModules`) and is warmed and settled by count (`frameScaling`): at 10 kb the climb to TurboFan
// runs thousands of frames at 55 to 270 bytes a frame.

import { describe, expect, it } from 'vitest';
import { controlCaught, describeScaling, frameScaling, isFlat, leaking } from './heap.js';
import { LETTER_VIEWS, lettersFrame } from './streams-letters-cases.js';

const N = 1_000;

describe("the streamed letters' pans allocate nothing per frame once their regions are held", () => {
  for (const at of LETTER_VIEWS)
    it(`a pass and a frame, ${at}`, async () => {
      const g = await lettersFrame(at);
      // Not vacuous: the letters drawn (100 b), their cells (10 kb).
      g.frame();
      if (at === '100 b') expect(g.backend.textItems.n).toBeGreaterThan(50);
      else expect(g.backend.rectItems.n).toBeGreaterThan(500);
      const asked = g.requests();
      const r = frameScaling(g.frame, N);
      expect(isFlat(r), describeScaling(r, N)).toBe(true);
      expect(g.requests()).toBe(asked);
    }, 300_000);

  it('the same frame allocating an object a step is caught (positive control)', async () => {
    const g = await lettersFrame('100 b');
    expect(controlCaught(() => frameScaling(leaking(g.frame), N))).toBe(true);
  }, 300_000);
});
