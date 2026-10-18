/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The browser heap guards' measurement (RB0). A page cannot see a scavenge, so one run could hide its garbage behind a
// collection. Two run lengths, `n` and `4n` steps, each after a warm-up and a forced collection: allocation per step
// shows as growth that scales with the length, and a flat path as growth that does
// not. Run with the browser heap config (`HEAP=1`: `gc` exposed, a fixed young generation, precise memory info).
//
// The two lengths alone alias: a step whose `n` steps fill one young generation and a little more collects in both
// runs, and both end under the limit (1 KB a step over 64,010 steps did, 2026-10-11). A collection frees the steps'
// garbage, so the heap, read 16 times in a run, falls when one ran; a fall refuses the run, as the Node guard does.

/** The growth a guard allows over one measured run. */
export const LIMIT = 1 << 20;

const gc = () => (globalThis as { gc?: () => void }).gc?.();
const heap = () => (performance as unknown as { memory: { usedJSHeapSize: number } }).memory.usedJSHeapSize;

/** The growth over `n` steps after a warm-up, or `undefined` when the heap fell during them (a collection ran). */
function growth(n: number, step: (i: number) => void): number | undefined {
  // Garbage an earlier test or file left in this isolate, collected after the warm-up, would deoptimize the hot code
  // its maps share (core/_test/model/heap.ts); collected first, the warm-up optimizes again.
  gc();
  for (let i = 0; i < n; i++) step(i);
  gc();
  gc();
  // A fall between two reads is missed only when the steps between them allocate more than the collection freed,
  // about a young generation: a run would have to allocate 16 of them.
  const every = Math.max(1, n >> 4);
  const before = heap();
  let last = before;
  for (let i = 0; i < n; i++) {
    step(i);
    if (i % every === 0) {
      const h = heap();
      if (h < last) return undefined;
      last = h;
    }
  }
  const after = heap();
  return after < last ? undefined : after - before;
}

/**
 * The heap growth over `n` steps and over `4n`, and the part of the longer run's that the length explains. Throws when
 * a collection ran during a measured run, after one retry: shorten the run, the growth it would report is unknowable.
 */
export function heapScaling(n: number, step: (i: number) => void): { short: number; long: number; scaled: number } {
  if (!(globalThis as { gc?: unknown }).gc) throw new Error('Run with HEAP=1 (the browser heap config).');
  for (let attempt = 0; attempt < 2; attempt++) {
    const short = growth(n, step);
    const long = short === undefined ? undefined : growth(4 * n, step);
    if (short !== undefined && long !== undefined) return { short, long, scaled: long - short };
  }
  throw new Error('A GC ran during a measured run, twice: it hides allocation. Make the run shorter.');
}

/** True when neither run grew past the limit and the longer run did not grow with its length. */
export const isFlat = (r: { short: number; long: number; scaled: number }): boolean =>
  r.short < LIMIT && r.long < LIMIT && r.scaled < LIMIT / 2;

const leaked: unknown[] = [null];

/**
 * The positive controls' allocation, as the Node guards' (core/_test/model/heap.ts): garbage at the next call. With
 * `slots`, an array of that many slots (4 bytes each with pointer compression, Chromium's default, 8 without) instead
 * of the object of eight fields (about 44 bytes, 88 without).
 */
export function leak(k: number, slots = 0): void {
  leaked[0] = slots > 0 ? new Array(slots) : { a: k, b: k, c: k, d: k, e: k, f: k, g: k, h: k };
}

/**
 * The array slots a step must allocate for a run of `n` steps to grow at least twice what `isFlat` allows: a page
 * cannot count on a collection to catch a small leak (the young generation is 64 MB), so the longer run's extra `3n`
 * steps must pass `LIMIT / 2` by themselves.
 */
export const plantedSlots = (n: number): number => Math.ceil(LIMIT / (3 * n) / 4);

/**
 * The positive control's step (the guard template): the guard's own step, then `leak(i)`. Pass the measured run's
 * length `n`: the planted garbage is then sized to it (`plantedSlots`), so a short run's control is not blind.
 */
export const leaking = (step: (i: number) => void, n?: number) => {
  const slots = n ? plantedSlots(n) : 0;
  return (i: number): void => {
    step(i);
    leak(i, slots);
  };
};

/**
 * The guard template's control (packages/sequences/AGENTS.md, Tests), as the Node guards' `controlCaught`: `measure` runs
 * `leaking(step)` through the file's own warm-up and measurement. Caught when not flat (or a number of `limit` or more,
 * for a guard judged by its own measure: what a retention guard finds still reachable), or when the measurement threw
 * because the heap fell during a run (a page cannot see a scavenge: a collection the leak caused counts as detection).
 */
export function controlCaught(
  measure: () => { short: number; long: number; scaled: number } | number,
  limit = LIMIT,
): boolean {
  try {
    const r = measure();
    return typeof r === 'number' ? r >= limit : !isFlat(r);
  } catch (e) {
    if (e instanceof Error && /GC ran during a measured run/.test(e.message)) return true;
    throw e;
  }
}
