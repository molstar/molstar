/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The values draw's case, for its heap guard and its quick twin (streams-values.heap.ts, streams-values.quick.ts): ST1's
// values draw over held regions beside the genes, variants and letters, a pan's pass and frame at 100 kb and DPR 2.
// Loads the viewer afresh (`vi.resetModules`).

import { vi } from 'vitest';

export async function valuesFrame() {
  vi.resetModules();
  const { progressGenome } = await import('./stream-progress-fixtures.js');
  const g = await progressGenome(100_000_000, 100_100_000, { values: true, items: true, sequence: true });
  await g.settle();
  let i = 0;
  const frame = () => {
    g.vp.panBy((i++ & 1) === 0 ? 9 : -9, 0, 'pointer');
    g.m.pass();
    g.frame();
  };
  return { g, frame };
}
