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

// The scaling tier of this package (`test:scaling`, opt in, never in CI), as prism-ui's `src/_test/scaling.ts` has it
// (packages/prism/AGENTS.md, bulk operations): an async path run at n and at 4n, the fastest time and the bytes
// allocated compared, both under 8x (linear reads about 4, a quadratic about 16); a deliberately quadratic control must
// be caught. Self-contained: this package imports nothing of Prism.

export interface AsyncCase<S> {
  /** Makes the state (untimed). */
  setup(n: number): Promise<S> | S;
  /** The timed path. */
  run(state: S): Promise<void> | void;
}

export interface Scaling {
  n: number;
  factor: number;
  ms: [number, number];
  bytes: [number, number];
  timeRatio: number;
  allocRatio: number;
}

const used = () => getHeapStatistics().used_heap_size;

/** Bytes allocated while `fn` runs: the heap's growth with what collections freed counted back. */
async function allocated(fn: () => Promise<void> | void): Promise<number> {
  const profiler = new GCProfiler();
  profiler.start();
  const before = used();
  await fn();
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

const median = (xs: number[]) => [...xs].sort((a, b) => a - b)[Math.floor(xs.length / 2)] as number;

async function measure<S>(c: AsyncCase<S>, n: number, reps: number) {
  const times: number[] = [];
  const bytes: number[] = [];
  for (let r = 0; r < reps; r++) {
    const s = await c.setup(n);
    let ms = 0;
    bytes.push(
      await allocated(async () => {
        const t0 = performance.now();
        await c.run(s);
        ms = performance.now() - t0;
      }),
    );
    times.push(ms);
  }
  return { ms: Math.min(...times), bytes: median(bytes) };
}

/** Runs the case at n and 4n (`reps` each, after one warm run of each). */
export async function scaling<S>(c: AsyncCase<S>, n: number, options: { factor?: number; reps?: number } = {}) {
  const factor = options.factor ?? 4;
  const reps = options.reps ?? 5;
  await measure(c, n, 1);
  await measure(c, factor * n, 1);
  const a = await measure(c, n, reps);
  const b = await measure(c, factor * n, reps);
  return {
    n,
    factor,
    ms: [a.ms, b.ms],
    bytes: [a.bytes, b.bytes],
    timeRatio: b.ms / Math.max(a.ms, 1e-3),
    allocRatio: b.bytes / Math.max(a.bytes, 1),
  } satisfies Scaling;
}

export const isLinear = (s: Scaling, band = 2 * s.factor) => s.timeRatio < band && s.allocRatio < band;

export const describeScaling = (s: Scaling): string =>
  `n ${s.n} → ${s.factor * s.n}: ${s.ms[0].toFixed(2)} → ${s.ms[1].toFixed(2)} ms (×${s.timeRatio.toFixed(1)}), ` +
  `${(s.bytes[0] / 1024).toFixed(0)} → ${(s.bytes[1] / 1024).toFixed(0)} KB (×${s.allocRatio.toFixed(1)})`;

/** The positive control: a copy of the whole list per item, n items (about factor² when run at n and 4n). */
export const quadraticControl: AsyncCase<number> = {
  setup: (n) => n,
  run: (n) => {
    let list: number[] = [];
    for (let i = 0; i < n; i++) list = [...list, i];
    if (list.length !== n) throw new Error('control');
  },
};

/** Asserts the case grows about linearly; `SCALING_REPORT=<file>` appends its numbers there. */
export function checkLinear(name: string, s: Scaling, linear = true): void {
  if (process.env.SCALING_REPORT) appendFileSync(process.env.SCALING_REPORT, `${name} | ${describeScaling(s)}\n`);
  expect(isLinear(s), `${name}: ${describeScaling(s)}`).toBe(linear);
}
