/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The quick allocation tier's measurement (SPIKE; packages/sequences/AGENTS.md, Tests), for every package: bytes
// allocated per step, read directly instead of inferred from heap growth. `v8.GCProfiler` records the heap spaces
// before and after every collection in a window, so the bytes allocated over the window are Σ(before − the previous
// after) plus the final delta, whatever collected in between: a window may collect as often as it likes (the quick
// config runs a 1 MB young generation, `viewer/vitest.perf.config.ts`). Exact to the byte in Node 25: `leak`'s object
// reads 88.0 bytes a step. Objects only: the code and trusted spaces (code, bytecode, deoptimization data) are the
// compiler's, reported apart (`codePerStep`). The window's own statistics objects are measured once (`overhead`) and
// subtracted. `getHeapStatistics().total_allocated_bytes` (Node 25 only) counts by linear allocation buffer and swings
// by tens of bytes a step: not used.
//
// The warm-up is by count, not time: `warmSteps`, then windows of `steps` until two in a row are quiet (under `QUIET`)
// and compiled nothing (the code and trusted spaces grew by under `COMPILED`), at most `maxWarmSteps`. Windows that
// agree on a steady allocation do not end it: a step sits in a lower tier at a steady 16 or 130 bytes for thousands of
// steps, and a lull between tiers compiles nothing (2026-10-15); a step that allocates once optimized runs to the cap
// and is measured. The quick config lowers V8's tier-up thresholds (`--invocation-count-for-turbofan=300
// --minimum-invocations-after-ic-update=50`): with the defaults a step stays in its middle tier (Maglev) for 15,000 to
// 60,000 steps (the wrapped lanes' draw 130 bytes a step, the hover 9 to 26); with them most are optimized in a few
// thousand, and the steady state is the defaults' (2026-10-14, Node 24 and 25, M4). `--single-threaded` compiles on the
// main thread, so tier-up falls on the same step whatever else the machine runs.

import { GCProfiler, getHeapSpaceStatistics, getHeapStatistics } from 'node:v8';

/** Below this a window counts as allocating nothing, in bytes a step. */
export const QUIET = 0.5;
/** A measured window reading below minus this many bytes a step lost a collection's account (re-run). */
export const TOLERANCE = 1;

export interface AllocOptions {
  /** Steps run before the first window, unmeasured. Default 2,000. */
  warmSteps?: number;
  /** Steps per window. Default 1,000. */
  steps?: number;
  /** Measured windows after the warm-up; the result is their median. Default 3. */
  windows?: number;
  /** The warm-up runs at least this many steps. Default 0. */
  minWarmSteps?: number;
  /**
   * The warm-up's cap, in steps. Default 30,000. A step that allocates by design never reads quiet and runs it: give it
   * a smaller cap.
   */
  maxWarmSteps?: number;
}

export interface Alloc {
  /** The median of the measured windows' bytes a step, objects only (code and bytecode apart). */
  bytesPerStep: number;
  /** Each measured window's bytes a step. */
  windows: number[];
  /** The largest minus the smallest window. */
  spread: number;
  /** Code, bytecode and other compiler output a step (the code and trusted spaces), over the measured windows. */
  codePerStep: number;
  /** Steps run before the measured windows. */
  warmed: number;
  /** How the warm-up ended: two quiet windows (compiling nothing), or the cap. */
  settled: 'quiet' | 'capped';
  /** Collections inside the measured windows, all accounted for. */
  collections: number;
  /** A line for a report: the windows, the warm-up and the collections. */
  note: string;
}

/** The compiler's spaces: code, bytecode and their metadata, apart from the objects a step allocates. */
const CODE = /code|trusted|read_only|shared/;

interface Space {
  spaceName?: string;
  spaceUsedSize?: number;
  space_name?: string;
  space_used_size?: number;
}

/** Used bytes in the compiler's spaces (`CODE`): GCProfiler's names are camelCase, v8's snake_case. */
function codeUsed(spaces: readonly Space[]): number {
  let n = 0;
  for (const s of spaces)
    if (CODE.test(s.spaceName ?? s.space_name ?? '')) n += s.spaceUsedSize ?? s.space_used_size ?? 0;
  return n;
}

/**
 * The bytes (objects, code) `n` calls of `step` allocate, and the collections that ran meanwhile. The heap's total is
 * read with `getHeapStatistics` (one object, made after the read) and the code spaces apart: `getHeapSpaceStatistics`
 * makes an object per space between its reads, so the objects' share of a per-space total moved by a kilobyte when
 * V8 pretenured them, a window's constant that a calibration made before could not subtract (2026-10-14).
 */
function window(
  step: (i: number) => void,
  start: number,
  n: number,
): { bytes: number; code: number; collections: number } {
  const profiler = new GCProfiler();
  const code0 = codeUsed(getHeapSpaceStatistics());
  const total0 = getHeapStatistics().used_heap_size;
  profiler.start();
  for (let i = start; i < start + n; i++) step(i);
  const total1 = getHeapStatistics().used_heap_size;
  const stats = profiler.stop()?.statistics ?? [];
  const code1 = codeUsed(getHeapSpaceStatistics());
  let total = 0;
  let code = 0;
  let prevTotal = total0;
  let prevCode = code0;
  for (const s of stats) {
    total += s.beforeGC.heapStatistics.usedHeapSize - prevTotal;
    code += codeUsed(s.beforeGC.heapSpaceStatistics as Space[]) - prevCode;
    prevTotal = s.afterGC.heapStatistics.usedHeapSize;
    prevCode = codeUsed(s.afterGC.heapSpaceStatistics as Space[]);
  }
  total += total1 - prevTotal;
  code += code1 - prevCode;
  return { bytes: total - code, code, collections: stats.length };
}

const noop = () => {};
/**
 * The measurement's own allocation in a window (its statistics' objects), subtracted from every window: measured once
 * per process, so a case and its control subtract the same (a window of 5 steps, a step of milliseconds, turns a
 * calibration's 200 bytes into 40 a step).
 */
let own = Number.NaN;
function overhead(): number {
  if (!Number.isNaN(own)) return own;
  // The first windows run the measurement in its lower tiers, which allocate more: the median of the last 101.
  const b: number[] = [];
  for (let k = 0; k < 300; k++) b.push(window(noop, 0, 0).bytes);
  own = b.slice(-101).sort((x, y) => x - y)[50];
  return own;
}

/**
 * The deep guards' warm-up by count (`viewer/_test/heap.ts` `warmFrame`, `settleFrame`): windows of `n` steps from
 * `start` until they allocate under `quiet` bytes a step and compile nothing over `SETTLE_SPAN` steps (three windows at least), read exactly as `measureBytesPerStep` reads them, at
 * most `cap` steps. Under V8's default tuning a step stays in its middle tier for 15,000 to 60,000 steps, allocating
 * 9 to 130 bytes a step there (2026-10-14): a heap guard measured inside that climb, which a timed warm-up ended in
 * under load. Returns the next step's index and whether it settled (a step that allocates once optimized never does:
 * the guard then measures it, and catches it).
 */
export function settleByCount(
  step: (i: number) => void,
  start: number,
  n: number,
  quiet: number,
  cap: number,
): { next: number; settled: boolean } {
  const own = overhead();
  let i = start;
  let quietRun = 0;
  while (i - start < cap) {
    const w = window(step, i, n);
    const b = (w.bytes - own) / n;
    i += n;
    // Quiet and compiling nothing: a lull between tiers allocates nothing for a while too (`COMPILED`), and the climb's
    // allocation comes in bursts thousands of steps apart (a ring frame of 120,000 positions: quiet 500-step windows
    // between bursts of 300 bytes a step, 2026-10-15), so the quiet run spans `SETTLE_SPAN` steps, three windows at least.
    quietRun = b < quiet && w.code <= COMPILED ? quietRun + n : 0;
    if (quietRun >= Math.max(3 * n, SETTLE_SPAN)) return { next: i, settled: true };
  }
  return { next: i, settled: false };
}

/** `settleByCount`'s least quiet span, in steps. */
const SETTLE_SPAN = 6_000;

/**
 * A window whose code and trusted spaces grew by more than this (bytes) compiled something: V8 is still tiering the
 * step up, and its allocation may change (a compile installs kilobytes; a quiet window's code reads 0).
 */
const COMPILED = 512;

/**
 * Bytes allocated per call of `step` (it gets a running index), after a warm-up by count. Throws only on a step's own
 * throw: a step that keeps allocating is measured, and its `settled` says the warm-up did not end quiet.
 */
export function measureBytesPerStep(step: (i: number) => void, options: AllocOptions = {}): Alloc {
  const { warmSteps = 2_000, steps = 1_000, windows = 3, minWarmSteps = 0, maxWarmSteps = 30_000 } = options;
  const own = overhead();
  let i = 0;
  for (; i < warmSteps; i++) step(i);
  let settled: Alloc['settled'] = 'capped';
  let warmed = 0;
  let per: number[] = [];
  let collections = 0;
  let code = 0;
  let rejected = 0;
  let resumed = 0;
  for (;;) {
    let last = Number.NaN;
    let lastCompiled = true;
    settled = 'capped';
    while (i < maxWarmSteps) {
      const w = window(step, i, steps);
      const b = (w.bytes - own) / steps;
      const compiled = w.code > COMPILED;
      i += steps;
      // Only two quiet windows in which V8 compiled nothing end the warm-up early. Windows that agree on a steady
      // allocation do not: a step sits in a lower tier for thousands of steps at a steady 16 or 130 bytes (a boxed
      // double, 2026-10-15), and a lull between tiers compiles nothing. A step that allocates once optimized runs to the
      // cap.
      if (i >= minWarmSteps && !compiled && !lastCompiled && b < QUIET && last < QUIET) {
        settled = 'quiet';
        break;
      }
      last = b;
      lastCompiled = compiled;
    }
    warmed = i;
    per = [];
    collections = code = rejected = 0;
    let recompiled = false;
    for (let w = 0; w < windows; w++) {
      const r = window(step, i, steps);
      i += steps;
      const b = (r.bytes - own) / steps;
      // A window that reads less than nothing lost a collection's account (seen once in Node 24: -834 bytes a step, a
      // major GC's sweeping finished inside the window): measured again, at most `windows` times in all. No full
      // collection before the windows: the steps after one rebuilt what it cleared and recompiled (7 to 17 bytes a step
      // over the next two windows of a text frame).
      if (b < -TOLERANCE && rejected < windows) {
        rejected++;
        w--;
        continue;
      }
      per.push(b);
      code += r.code;
      collections += r.collections;
      recompiled ||= r.code > COMPILED;
    }
    // A quiet warm-up whose measured windows then compiled ended in a lull: warmed on, and measured again.
    if (!(settled === 'quiet' && recompiled && i < maxWarmSteps && resumed < 4)) break;
    resumed++;
  }
  const codePerStep = code / (windows * steps);
  const sorted = [...per].sort((a, b) => a - b);
  const bytesPerStep = sorted[sorted.length >> 1];
  const spread = sorted[sorted.length - 1] - sorted[0];
  const note =
    `windows ${per.map((b) => b.toFixed(2)).join(' / ')} B (${steps} steps each), warmed ${warmed} steps (${settled}), ` +
    `${collections} GC(s) accounted, code ${codePerStep.toFixed(2)} B${rejected ? `, ${rejected} negative window(s) re-run` : ''}` +
    `${resumed ? `, warmed on ${resumed} time(s) after a lull` : ''}`;
  return { bytesPerStep, windows: per, spread, codePerStep, warmed, settled, collections, note };
}
