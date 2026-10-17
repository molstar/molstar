/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { getHeapSpaceStatistics } from 'node:v8';

// Seq*'s guard harness for Prism's tiers (packages/prism/AGENTS.md, Tests): `quickCase` for the quick allocation tier
// (bytes a step against `perf-baseline.json`, its positive control built in) and the heap guards' scaling (n and 4n
// steps, `leaking` for the control). Test code only: the harness is Seq*'s `_test/model` (packages/sequences/AGENTS.md).

export { quickCase } from '../../../../sequences/core/_test/model/quick.js';
export {
  controlCaught,
  describeScaling,
  heapScaling,
  isFlat,
  leaking,
} from '../../../../sequences/core/_test/model/heap.js';

/**
 * The heap in use outside code (V8's code and trusted spaces apart: jsdom's selector engine compiles a function and a
 * regexp a query, 1 KB a cycle of code that is not the view's; Seq*'s quick tier counts code apart likewise).
 */
export function heapUsedOutsideCode(): number {
  let used = 0;
  for (const s of getHeapSpaceStatistics()) if (!/code|trusted/.test(s.space_name)) used += s.space_used_size;
  return used;
}

const gcNow = () => {
  const gc = (globalThis as { gc?: () => void }).gc;
  if (!gc) throw new Error('Run with --expose-gc (pnpm test:heap).');
  gc();
  gc();
};

/**
 * For a step that allocates by design (a keystroke's matches, jsdom's nodes): what it keeps, the heap after a full
 * collection, over `n` steps and over `4n` more. Flat when neither grew past 1 MB and the longer kept at most 64 KB more
 * than the shorter (Seq*'s limits, heap.ts). `heapScaling` counts garbage as kept (it collects nothing while it
 * measures), so it suits only steps that allocate nothing.
 */
export function retainedScaling(n: number, step: (i: number) => void): { short: number; long: number; scaled: number } {
  gcNow();
  const h0 = heapUsedOutsideCode();
  for (let i = 0; i < n; i++) step(i);
  gcNow();
  const h1 = heapUsedOutsideCode();
  for (let i = 0; i < 4 * n; i++) step(i);
  gcNow();
  const h2 = heapUsedOutsideCode();
  const short = h1 - h0;
  const long = h2 - h1;
  return { short, long, scaled: long - short };
}

/**
 * Whether `retainedScaling` over `n` steps kept nothing a step: neither run past 1 MB, and the longer (3n steps more)
 * kept at most 64 KB plus 8 bytes a step more than the shorter. A reading after a full collection still moves with
 * compiled code and inline caches (about 16 bytes a step on the scene tree's pick in jsdom, where a heap snapshot's diff
 * found no object kept, 2026-10-17); a step that keeps 32 bytes or more is caught.
 */
export const isRetainedFlat = (r: { short: number; long: number; scaled: number }, n: number): boolean =>
  r.short < 1 << 20 && r.long < 1 << 22 && r.scaled < (64 << 10) + 8 * 3 * n;
