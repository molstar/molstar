/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The axis's performance representative (_spec/README.md "Render", Tracks layout): a pan frame over a reversed axis of
// 150,000 single-position pieces (EC-33) under the pan target (16 ms, a frame's maximum), the best of three; quick
// without the bound, heavy with it (`itWithin`). Its walks (median, p95, maximum) are the heavy frame-time scenario
// `reversed, 150,000 pieces`, and its allocation the opt-in tier's.

import { bestOf, itWithin } from '@molstar/test-harness/timing';
import { beforeAll, describe, expect } from 'vitest';
import { PIECES, piecesSpec } from './fixtures.js';
import { type Settled, settle } from './helpers.js';

const Size = { width: 1440, height: 160 };

describe('a reversed axis of 150,000 single-position pieces', () => {
  let s: Settled;
  beforeAll(async () => {
    s = await settle(piecesSpec());
  }, 60_000);

  itWithin('a pan frame lays out and draws under 16 ms', async (within) => {
    const at = (c: number) => ({ size: Size, dpr: 2, viewport: { start: c - 500, end: c + 500 } });
    // warm, with the data of the windows in the store
    for (let c = PIECES / 2; c < PIECES / 2 + 40; c += 4) await s.frame(at(c));
    let k = 0;
    const f = bestOf(within, 16, () => s.renderer.frame(at(PIECES / 2 + (k++ % 10) * 4 + 1)));
    expect(f.rebuilt).toBe(true);
    expect(f.layout.sections[0].breaks.length / 3).toBeGreaterThan(900);
  });
});
