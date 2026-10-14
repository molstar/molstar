/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The heap guards' measurement (RB0), for every package. The heap config (vitest.shared.ts) runs a guard with
// `--expose-gc` and a fixed 64 MB young generation, so a hot path that allocates only a little leaves its garbage in
// the heap, where the growth shows it. A guard whose run is long enough to trigger a collection would hide that
// allocation, so `heapGrowth` fails on any GC during the measured run.
//
// One run against a fixed limit misses a small allocation per step: 210 bytes a frame over 1,000 frames is 210 KB
// (BK1). A guard of a repeated step uses `heapScaling`, as the browser guard does (viewer/_test/canvas/heap-browser.ts):
// runs of `n` and `4n` steps, and the part of the longer run's growth that its length explains. Node reads the heap
// exactly, so its bound is tighter than the browser's. `heapGrowth` stays for a run that is not a repetition (a build).

import { GCProfiler } from "node:v8";

/** The growth a guard allows over one measured run. */
export const LIMIT = 1 << 20;

/** The growth `heapScaling`'s longer run may add to the shorter's: a 1,000-step guard catches 22 bytes a step. */
export const SCALED_LIMIT = 64 << 10;

const gc = (globalThis as { gc?: () => void }).gc;

/**
 * Heap growth over a second run of `run` (the first warms up). Throws when a GC ran during the measured run, after
 * one retry (the first measured run can still allocate in a lower tier, or meet a stray collection request): shorten
 * the run until it does not collect, because the growth it would report is then unknowable.
 */
export function heapGrowth(run: () => void): number {
  if (!gc) throw new Error("Run with --expose-gc (the heap config: pnpm test:heap).");
  let collections = 0;
  for (let attempt = 0; attempt < 2; attempt++) {
    // Garbage a guard before this one left would be collected after the warm-up, and its maps would deoptimize the
    // hot code they share. Collected first, the warm-up optimizes again; each guard is independent of the others.
    gc();
    run();
    gc();
    gc();
    const profiler = new GCProfiler();
    profiler.start();
    const before = process.memoryUsage().heapUsed;
    run();
    const after = process.memoryUsage().heapUsed;
    collections = profiler.stop()?.statistics.length ?? 0;
    if (collections === 0) return after - before;
  }
  throw new Error(
    `${collections} GC(s) ran during the measured run, twice: it hides allocation. Make the run shorter.`,
  );
}

/** The growth over `n` steps after a warm-up of `n`, or `undefined` when a GC ran during the measured steps. */
function stepGrowth(n: number, step: (i: number) => void, gc: () => void): number | undefined {
  // As heapGrowth's: collected before the warm-up and after it, so earlier garbage neither deoptimizes nor counts.
  gc();
  for (let i = 0; i < n; i++) step(i);
  gc();
  gc();
  const profiler = new GCProfiler();
  profiler.start();
  const before = process.memoryUsage().heapUsed;
  for (let i = 0; i < n; i++) step(i);
  const after = process.memoryUsage().heapUsed;
  return (profiler.stop()?.statistics.length ?? 0) === 0 ? after - before : undefined;
}

/** `heapScaling`'s measurement: the growth over `n` steps (`short`), over `4n` (`long`), and `long - short`. */
export interface Scaling {
  short: number;
  long: number;
  scaled: number;
}

/**
 * The heap growth over `n` steps and over `4n`, each after a warm-up of the same length; `step` gets the step's index.
 * Allocation per step shows in `scaled`, three times `n` steps' worth; allocation per run (a result kept, a first-use
 * cache) cancels. Throws when a GC ran during a measured run, after one retry of both: shorten `n`.
 */
export function heapScaling(n: number, step: (i: number) => void): Scaling {
  if (!gc) throw new Error("Run with --expose-gc (the heap config: pnpm test:heap).");
  for (let attempt = 0; attempt < 2; attempt++) {
    const short = stepGrowth(n, step, gc);
    const long = short === undefined ? undefined : stepGrowth(4 * n, step, gc);
    if (short !== undefined && long !== undefined) return { short, long, scaled: long - short };
  }
  throw new Error("GC(s) ran during the measured run, twice: it hides allocation. Make the run shorter.");
}

const leaked: unknown[] = [null];

/**
 * The positive controls' allocation, for a guard's step to call: an object of eight small-integer fields (88 bytes in
 * Node, measured: twice what a 500-step guard catches, so a frame's own noise cannot hide it), garbage at the next
 * call. A measured run collects nothing, so garbage counts as a kept object would, and the bytes a step are the same in
 * every run. A list kept instead grows by reallocation, its backing store copied at 1.5 times its length: one copy
 * inside `heapScaling`'s shorter measured run cancels the longer run's growth, and where it falls depends on how many
 * steps a timed warm-up made (2026-10-13: streams-port's and the TP53 folded guard's controls, an object kept a step,
 * measured flat).
 */
export function leak(k: number): void {
  leaked[0] = { a: k, b: k, c: k, d: k, e: k, f: k, g: k, h: k };
}

/** True when neither run grew past `LIMIT` and the longer grew at most `SCALED_LIMIT` beyond the shorter. */
export const isFlat = (r: Scaling): boolean => r.short < LIMIT && r.long < LIMIT && r.scaled < SCALED_LIMIT;

/** The positive control's step (the guard template): the guard's own step, then `leak(i)`. */
export const leaking =
  (step: (i: number) => void) =>
  (i: number): void => {
    step(i);
    leak(i);
  };

/** The measurements' refusal of a run that collected (`heapGrowth`, `heapScaling`). */
const COLLECTED = /GC\(s\) ran during the measured run/;

/**
 * The guard template's control (lib/seqstar/AGENTS.md, Tests): `measure` runs `leaking(step)` through the file's own
 * warm-up, settle and measurement, exactly as its guards run `step`. Caught when the measurement is not flat (a
 * `Scaling`, or a number of `limit` or more: by default a growth past `LIMIT`; a guard judged per position or per cell
 * passes its own measure and bound, LA-24), or when it refused the run because a collection ran during it: the leak
 * filled the young generation, which is detection too. Any other throw is the guard's own failure, rethrown. Every
 * `*.heap.ts` calls it (the spec check).
 */
export function controlCaught(measure: () => Scaling | number, limit = LIMIT): boolean {
  try {
    const r = measure();
    return typeof r === "number" ? r >= limit : !isFlat(r);
  } catch (e) {
    if (e instanceof Error && COLLECTED.test(e.message)) return true;
    throw e;
  }
}

/** A failed guard's message: the three growths in KB, and the per-step bytes the scaled part suggests. */
export const describeScaling = (r: Scaling, n: number): string =>
  `short ${(r.short / 1024).toFixed(1)} KB, long ${(r.long / 1024).toFixed(1)} KB, scaled ` +
  `${(r.scaled / 1024).toFixed(1)} KB (${(r.scaled / (3 * n)).toFixed(1)} bytes a step)`;
