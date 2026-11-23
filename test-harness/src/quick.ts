/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The quick tier's case (`test:perf`; packages/sequences/AGENTS.md, Tests), for every package: a step measured in bytes a
// step (`measureBytesPerStep`, alloc.ts) against the package's baseline (`_test/perf-baseline.json`), then its positive
// control, `leaking(step)`, on the same state. A case over its baseline by more than the threshold fails (`PERF_THRESHOLD`,
// 4 bytes a step, or a tenth of a baseline above 40: `allowance`); a smaller move is listed by the report (perf-reporter.ts), which also rewrites the
// baseline under `--update` (or `PERF_UPDATE=1`). The config provides the baseline and the mode (vitest.ts
// `perfTests`). Every `*.quick.ts` calls `quickCase` (the spec check); a case without its control says why.
// A case that allocates kilobytes a step has its control scaled to it (`scaledBytes`).

import { expect, inject } from 'vitest';
import { type Alloc, type AllocOptions, measureBytesPerStep } from './alloc.js';
import { leak, leaking } from './heap.js';

export { type Alloc, type AllocOptions, measureBytesPerStep } from './alloc.js';

/** What a case leaves on its test's meta for the report. */
export interface AllocRecord {
  id: string;
  bytesPerStep: number;
  /** The baseline, absent for a new case or a positive control. */
  baseline?: number;
  /** The positive control's bytes a step over the case's (88 for `leak`). */
  control?: number;
  note: string;
}

declare module 'vitest' {
  interface ProvidedContext {
    perfBaseline: Record<string, number>;
    perfUpdate: boolean;
    perfThreshold: number;
  }
  interface TaskMeta {
    alloc?: AllocRecord[];
  }
}

/**
 * How far over its baseline a case may read: the threshold (bytes a step), or a tenth of the baseline where that is
 * more. A step that allocates kilobytes by design (a query run, a file read) reads a few percent apart run to run: its
 * major collections' sweeping finishes inside later windows (2026-10-14: a SeqQL run of about 10 KB read 9.3 to
 * 10.6 KB). The report (perf-reporter.ts) applies the same rule.
 */
export const allowance = (baseline: number, threshold: number): number => Math.max(threshold, 0.1 * baseline);

// `leak` optimized once per process: its object read 56 bytes in a lower tier, and a step of milliseconds gives the
// control only a few dozen calls (the alignment's and the molecule's checks, 2026-10-15).
for (let k = 0; k < 20_000; k++) leak(k);

/** `leak`'s object, in bytes (heap.ts): a control reads this much more a step than its case. */
export const LEAK_BYTES = 88;

/** A case reading above this (bytes a step) has its control scaled (`scaledBytes`). */
export const SCALED_ABOVE = 1_024;
/**
 * So has a case whose windows spread above this (bytes a step): half of what the 88-byte control must clear (66). The
 * quiet cases of the seq perf runs spread by 6 B at most (2026-11-23, under load).
 */
export const SPREAD_ABOVE = 32;
/** A scaled control is at least this many bytes. */
export const SCALED_MIN = 4_096;

/**
 * The control's bytes a step for a measured case: `LEAK_BYTES` for one that allocates little and steadily; for one that
 * reads above `SCALED_ABOVE` or spreads above `SPREAD_ABOVE`, a hundredth of its reading, at least `SCALED_MIN`, and
 * eight times its windows' spread where that is more. Over a step of 100 to 500 KB the windows wander by 20 to 250 bytes
 * a step (300 to 1,200 collections a window): an 88-byte object read from -13 to 105 over the case there (the bridge's
 * relayout, the categorical grid's pan, 2OCJ's hover), a 4 KB array is well above it; over a megabyte a frame the windows
 * are 25 KB apart (a layout's sweeping), and the control is eight times that.
 */
export function scaledBytes(r: Pick<Alloc, 'bytesPerStep' | 'spread'>): number {
  if (r.bytesPerStep <= SCALED_ABOVE && r.spread <= SPREAD_ABOVE) return LEAK_BYTES;
  return Math.max(SCALED_MIN, Math.ceil(0.01 * r.bytesPerStep), Math.ceil(8 * r.spread));
}

/** The array's header (`JSArray` and its `FixedArray`), in bytes, and a slot's. */
const ARRAY_HEADER = 16;
const SLOT = 8;
const garbage: unknown[] = [null];

/**
 * `step`, then an array of `bytes` (a step's garbage at the next one): the control of a case that allocates kilobytes. The
 * array's size is what V8 reports for it, header and slots, so the control reads as `bytes`.
 */
export const leakingBytes =
  (step: (i: number) => void, bytes: number) =>
  (i: number): void => {
    step(i);
    garbage[0] = new Array(Math.max(0, Math.round((bytes - ARRAY_HEADER) / SLOT))).fill(i);
  };

/** A case's options: the measurement's counts, and the control (on unless a reason says why not, and scaled to a case that allocates kilobytes). */
export type QuickOptions = AllocOptions & ({ control?: true } | { control: false; why: string });

/** The test's own task (`it(name, ({ task }) => …)`; `it.for` for a table, `it.each` passes no context). */
export interface QuickTask {
  meta: { alloc?: AllocRecord[] };
}

/**
 * Measures `step` (it gets a running index) and, unless `control: false`, `leaking(step)` right after it on the same
 * state; records both for the report. Fails when the case is over its baseline by more than the threshold (not under
 * `--update`), and when the control does not read at least three quarters of `leak`'s bytes over the case. `id` names
 * the case in the baseline: unique in the package, stable across runs.
 */
export function quickCase(
  task: QuickTask,
  id: string,
  step: (i: number) => void,
  options: QuickOptions = {},
): { r: Alloc; c?: Alloc } {
  const r = measureBytesPerStep(step, options);
  const baseline = inject('perfBaseline')[id];
  const records = (task.meta.alloc ??= []);
  records.push({ id, bytesPerStep: r.bytesPerStep, baseline, note: r.note });
  let c: Alloc | undefined;
  let controlBytes = LEAK_BYTES;
  if (options.control !== false) {
    // The step is warm by now, and the control never reads quiet: a window's warm-up for its wrapper (and `leak`), then
    // three measured (a lower tier's extra bytes only raise it, and the check is a floor).
    // At least 40 steps, for a step of milliseconds (5 a window).
    const n = options.steps ?? 1_000;
    const w = Math.max(n, 40);
    // A case that reads kilobytes a step, or whose windows wander, is measured to a few tens of bytes at best (hundreds
    // of collections a window): its control is scaled to it (`controlBytes`).
    controlBytes = scaledBytes(r);
    c = measureBytesPerStep(controlBytes === LEAK_BYTES ? leaking(step) : leakingBytes(step, controlBytes), {
      warmSteps: w,
      steps: n,
      windows: 3,
      maxWarmSteps: w,
    });
    records.push({
      id: `${id} +leak`,
      bytesPerStep: c.bytesPerStep,
      control: c.bytesPerStep - r.bytesPerStep,
      note: `${controlBytes === LEAK_BYTES ? '' : `${controlBytes} B control; `}${c.note}`,
    });
  }
  if (!inject('perfUpdate') && baseline !== undefined) {
    expect(r.bytesPerStep, `${id}: ${r.note}; baseline ${baseline}`).toBeLessThanOrEqual(
      baseline + allowance(baseline, inject('perfThreshold')),
    );
  }
  if (c) expect(c.bytesPerStep - r.bytesPerStep, `${id} +leak: ${c.note}`).toBeGreaterThan(0.75 * controlBytes);
  return { r, c };
}
