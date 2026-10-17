/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { appendFileSync } from 'node:fs';
import { GCProfiler, getHeapStatistics } from 'node:v8';
import { expect } from 'vitest';

// The scaling tier (`test:scaling`, packages/prism/AGENTS.md: bulk operations are one change): a bulk path (expand all,
// a select of many, a preset over many structures, a filter over many commands) run at n and at `factor`·n, its time
// and its allocation compared. Linear work grows by about the factor (n log n a little more); quadratic by its square.
// The ratio must stay under `band` (twice the factor by default: 8 for 4×), robust to a machine's noise; the positive
// control, a deliberately quadratic path, reads about 16 and is caught. Time is the fastest of `reps` runs (the noise is
// all upward); allocation is the bytes allocated over the run, read from the collector's own accounts (`GCProfiler`:
// what each collection freed counted back, as Seq*'s quick tier does, sequences/core/_test/model/alloc.ts).

/** A bulk path at a size: `setup` makes its state (untimed), `run` does the bulk operation on it (timed). */
export interface BulkCase<S> {
  setup(n: number): S;
  run(state: S): void;
  /** Called after each run (untimed): disposes what `setup` made. */
  teardown?(state: S): void;
}

export interface Scaling {
  n: number;
  factor: number;
  /** ms at n and at factor·n (the fastest run of each). */
  ms: [number, number];
  /** Bytes allocated at n and at factor·n (the run's median). */
  bytes: [number, number];
  timeRatio: number;
  allocRatio: number;
}

const used = () => getHeapStatistics().used_heap_size;

/** Bytes allocated by `fn`: the heap's growth, with what collections freed during it counted back. */
export function allocated(fn: () => void): number {
  const profiler = new GCProfiler();
  profiler.start();
  const before = used();
  fn();
  const after = used();
  const stats = profiler.stop()?.statistics ?? [];
  let bytes = 0;
  let from = before;
  for (const s of stats) {
    bytes += s.beforeGC.heapStatistics.usedHeapSize - from;
    from = s.afterGC.heapStatistics.usedHeapSize;
  }
  return bytes + after - from;
}

const median = (xs: number[]) => [...xs].sort((a, b) => a - b)[Math.floor(xs.length / 2)];

function measure<S>(c: BulkCase<S>, n: number, reps: number): { ms: number; bytes: number } {
  const times: number[] = [];
  const bytes: number[] = [];
  for (let r = 0; r < reps; r++) {
    const s = c.setup(n);
    let ms = 0;
    bytes.push(
      allocated(() => {
        const t0 = performance.now();
        c.run(s);
        ms = performance.now() - t0;
      }),
    );
    times.push(ms);
    c.teardown?.(s);
  }
  return { ms: Math.min(...times), bytes: median(bytes) };
}

/** Runs the case at n and factor·n (`reps` each, after one warm run of each). */
export function scaling<S>(c: BulkCase<S>, n: number, options: { factor?: number; reps?: number } = {}): Scaling {
  const factor = options.factor ?? 4;
  const reps = options.reps ?? 5;
  // Warm: both sizes once, so neither is measured unoptimized.
  measure(c, n, 1);
  measure(c, factor * n, 1);
  const a = measure(c, n, reps);
  const b = measure(c, factor * n, reps);
  return {
    n,
    factor,
    ms: [a.ms, b.ms],
    bytes: [a.bytes, b.bytes],
    timeRatio: b.ms / Math.max(a.ms, 1e-3),
    allocRatio: b.bytes / Math.max(a.bytes, 1),
  };
}

/** An async bulk path (a load, a preset, a select over a real plugin): as `BulkCase`. */
export interface AsyncBulkCase<S> {
  setup(n: number): Promise<S>;
  run(state: S): Promise<void>;
  teardown?(state: S): Promise<void> | void;
}

/**
 * `scaling` for an async path: time from the run's start to its end (the plugin's own work in it); allocation is not
 * measured (other tasks run meanwhile), its ratio reads 1.
 */
export async function scalingAsync<S>(
  c: AsyncBulkCase<S>,
  n: number,
  options: { factor?: number; reps?: number } = {},
): Promise<Scaling> {
  const factor = options.factor ?? 4;
  const reps = options.reps ?? 3;
  const at = async (size: number, k: number) => {
    let best = Number.POSITIVE_INFINITY;
    for (let r = 0; r < k; r++) {
      const s = await c.setup(size);
      const t0 = performance.now();
      await c.run(s);
      best = Math.min(best, performance.now() - t0);
      await c.teardown?.(s);
    }
    return best;
  };
  await at(n, 1);
  await at(factor * n, 1);
  const a = await at(n, reps);
  const b = await at(factor * n, reps);
  return { n, factor, ms: [a, b], bytes: [1, 1], timeRatio: b / Math.max(a, 1e-3), allocRatio: 1 };
}

/** Grows about linearly: both ratios under `band` (2·factor by default). */
export function isLinear(s: Scaling, band = 2 * s.factor): boolean {
  return s.timeRatio < band && s.allocRatio < band;
}

export const describeScaling = (s: Scaling): string =>
  `n ${s.n} → ${s.factor * s.n}: ${s.ms[0].toFixed(2)} → ${s.ms[1].toFixed(2)} ms (×${s.timeRatio.toFixed(1)}), ` +
  `${(s.bytes[0] / 1024).toFixed(0)} → ${(s.bytes[1] / 1024).toFixed(0)} KB (×${s.allocRatio.toFixed(1)})`;

/**
 * The positive control: a bulk path that copies the whole collection per item (the pattern this tier is for), n items.
 * Reads about factor² and must be caught.
 */
export const quadraticControl: BulkCase<number> = {
  setup: (n) => n,
  run: (n) => {
    let list: number[] = [];
    for (let i = 0; i < n; i++) list = [...list, i];
    if (list.length !== n) throw new Error('control');
  },
};

/** Asserts the case grows about linearly (`isLinear`); `SCALING_REPORT=<file>` appends its numbers there. */
export function checkLinear(name: string, s: Scaling, linear = true): void {
  if (process.env.SCALING_REPORT) appendFileSync(process.env.SCALING_REPORT, `${name} | ${describeScaling(s)}\n`);
  expect(isLinear(s), `${name}: ${describeScaling(s)}`).toBe(linear);
}
