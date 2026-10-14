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
} from "../../core/_test/model/heap.ts";

// The frame guards' shared run (the conformance kit's and the composer's): a frame is warmed by count (`FRAME_WARM`) and,
// where the draws compile concurrently, by time (`FRAME_WARM_MS`), then measured over `FRAME_N` steps and `4 * FRAME_N`.
// Trimmed on 2026-10-12 from 2,000 / 3 s / 1,000 (the tier took 568 s). streams-port.heap.ts and the other stream guards
// keep 1,000 steps, which catch about 22 bytes a step.
// `HEAP_DEEP=1` (`test:heap-deep`) runs the guards as before the trim: the full steps and warm-ups, which catch about
// 22 bytes a step (the trimmed run about 43), and the cases the trim dropped.
export const HEAP_DEEP =
  (globalThis as { process?: { env: Record<string, string | undefined> } }).process?.env.HEAP_DEEP === "1";
export const FRAME_N = HEAP_DEEP ? 1_000 : 500;
// The warm-up stays as before the trim in both modes: at 1 s and 1,000 frames the ring's and the composer's pan guards
// failed most runs (code still compiling), which cost more than the seconds saved (2026-10-12).
export const FRAME_WARM = 2_000;
export const FRAME_WARM_MS = 3_000;
/**
 * `settleFrame`'s bound after a timed warm-up (`FRAME_WARM_MS`): a whole frame's climb to TurboFan ran 7 s at the
 * composer's cell 0.4 and 6,000 frames past the warm-up on the ring of 120,000 positions (2026-10-13, M4).
 */
export const FRAME_SETTLE_MS = 15_000;

/** Below this growth a window of `FRAME_N` frames counts as settled, in bytes a frame (a fifth of what the measure catches). */
const SETTLED = 8;

/**
 * A frame warmed past V8's tier-up, for the guards that measure many frames apiece (the conformance kit's, a layer per
 * zoom): `FRAME_WARM` frames, then windows of `FRAME_N` until two in a row grow the heap by less than `SETTLED` bytes a
 * frame, for at most `FRAME_WARM_MS`. A draw called once a frame reaches TurboFan only after many frames: tens of
 * thousands where its frame is cheap (below the row's minimum) or after a new zoom or document deoptimized it, and until
 * then its middle tier and the module runner's export getters allocate (2026-10-13, Node 25 on an M4: ranges on
 * features' tiny draw 65 bytes a frame for 41,000 frames, then none; letters on all types 97 bytes for 8,000). Those
 * guards failed once `FRAME_N` was trimmed to 500 (6ee8a642), the measured frames then falling inside the climb; a time
 * warm-up for every frame instead would add minutes to the tier. A frame that allocates once optimized never settles:
 * it is warmed the whole `FRAME_WARM_MS`, then measured and caught (draws.heap.ts's positive control).
 */
export function warmFrame(frame: () => void): void {
  for (let i = 0; i < FRAME_WARM; i++) frame();
  settleFrame(frame, FRAME_WARM_MS);
}

/**
 * `warmFrame`'s second part: windows of `FRAME_N` frames until two in a row grow the heap by less than `SETTLED` bytes
 * a frame, for at most `ms`. The composer's guards call it after their time warm-up: a whole composer frame climbs
 * for longer than `FRAME_WARM_MS` on a slower machine (2026-10-13, M4: 17,000 frames, 7 s, over 400,000 positions at
 * cell 0.4; 33,000, 4.5 s, every type at 0.4), and a fixed time is fewer frames there.
 */
export function settleFrame(frame: () => void, ms: number): void {
  const gc = (globalThis as { gc?: () => void }).gc;
  if (!gc) throw new Error("Run with --expose-gc (the heap config: pnpm test:heap).");
  const t0 = performance.now();
  let settled = 0;
  while (settled < 2 && performance.now() - t0 < ms) {
    gc();
    const before = process.memoryUsage().heapUsed;
    for (let i = 0; i < FRAME_N; i++) frame();
    settled = process.memoryUsage().heapUsed - before < SETTLED * FRAME_N ? settled + 1 : 0;
  }
}
