/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The Node heap guards' measurement is one for every package: core/_test/model/heap.ts.
export {
  controlCaught,
  describeScaling,
  heapGrowth,
  heapScaling,
  isFlat,
  LIMIT,
  leak,
  leaking,
  SCALED_LIMIT,
  type Scaling,
} from '../../core/_test/model/heap.js';

// The quick tier's case (`test:perf`), for the `*.quick.ts` twins of the guards here.
export {
  type Alloc,
  measureBytesPerStep,
  type QuickOptions,
  type QuickTask,
  quickCase,
} from '../../core/_test/model/quick.js';

import { settleByCount } from '../../core/_test/model/alloc.js';
import { heapScaling, type Scaling } from '../../core/_test/model/heap.js';

// The frame guards' shared run (the conformance kit's, the composer's and every draw guard's): a frame is warmed by
// count (`FRAME_WARM`), then settled by count (`settleFrame`: until it stops allocating, read exactly), then measured
// over `FRAME_N` steps and `4 * FRAME_N`. Trimmed on 2026-10-12 from 1,000 steps (the tier took 568 s).
// `HEAP_DEEP=1` (`test:heap-deep`) runs the guards as before the trim: the full steps, which catch about 22 bytes a
// step (the trimmed run about 43), and the cases the trim dropped.
export const HEAP_DEEP =
  (globalThis as { process?: { env: Record<string, string | undefined> } }).process?.env.HEAP_DEEP === '1';
export const FRAME_N = HEAP_DEEP ? 1_000 : 500;
export const FRAME_WARM = 2_000;
/**
 * `settleFrame`'s cap, in steps. Under V8's default tuning (the heap tier's) a frame climbs to TurboFan over 15,000 to
 * 60,000 steps, allocating in its middle tier until then (2026-10-14, Node 24 and 25, M4: the wrapped lanes' draw 130
 * bytes a step for 14,000 steps; a hover 9 to 26 bytes for 20,000 to 60,000, longer after another case's draws). A
 * timed warm-up (the 3 s `FRAME_WARM_MS` before) ended inside that climb whenever the machine was loaded, which is why
 * the guards wanted a quiet machine and why a case could fail only after another (the "shared feedback" was the climb).
 */
export const FRAME_WARM_CAP = 100_000;

/** Below this a window of `FRAME_N` frames counts as settled, in bytes a frame (read exactly, not by heap growth). */
const SETTLED = 2;

/**
 * A frame warmed past V8's tier-up: `FRAME_WARM` frames, then `settleFrame`. For a guard that measures many frames
 * apiece (the conformance kit's, a layer per zoom) and every draw guard. A frame that allocates once optimized never
 * settles: it runs `FRAME_WARM_CAP` frames, then is measured and caught (draws.heap.ts's positive control).
 */
export function warmFrame(frame: (i: number) => void): void {
  for (let i = 0; i < FRAME_WARM; i++) frame(i);
  settleFrame(frame);
}

/**
 * Windows of `FRAME_N` frames until three in a row allocate under `SETTLED` bytes a frame (core/_test/model/alloc.ts
 * `settleByCount`, exact through `v8.GCProfiler`), at most `cap` frames. By count, so a loaded machine takes longer but
 * measures the same. Returns whether it settled.
 */
export function settleFrame(frame: (i: number) => void, cap = FRAME_WARM_CAP): boolean {
  return settleByCount(frame, 0, FRAME_N, SETTLED, cap).settled;
}

/**
 * The frame guards' measurement in one: `warm` steps (`FRAME_WARM`), `settleFrame`, then `heapScaling(n, …)`, `step`
 * getting a running index throughout (a step that alternates by index stays in phase).
 */
export function frameScaling(step: (i: number) => void, n = FRAME_N, warm = FRAME_WARM): Scaling {
  let i = 0;
  for (; i < warm; i++) step(i);
  i = settleByCount(step, i, FRAME_N, SETTLED, FRAME_WARM_CAP).next;
  return heapScaling(n, () => step(i++));
}

/** `cheapScaling`'s first warm-up, in steps. */
export const CHEAP_WARM = 60_000;

/**
 * `frameScaling` for a cheap step (tens of microseconds or less: hover frames, drags, small draws): `CHEAP_WARM` steps
 * first, and one `heapScaling` run discarded. Such a step climbs to TurboFan over 20,000 to 40,000 steps with quiet
 * windows on the way (0 to 5 bytes a step, then 10 to 30 again: 2026-10-15), so a settle from `FRAME_WARM` could end
 * inside the climb (27 to 57 bytes a step in the hover frames, flat afterwards); and the first measured run after a
 * collection rebuilds what it cleared (the discarded run).
 */
export function cheapScaling(step: (i: number) => void, n = FRAME_N): Scaling {
  let k = 0;
  for (; k < CHEAP_WARM; k++) step(k);
  heapScaling(n, () => step(k++));
  return frameScaling((i) => step(k + i), n, 0);
}
