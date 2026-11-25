/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Task } from '../task.js';
import { isTimingMode } from '@molstar/core/util/debug';

const hasPerformance = typeof performance !== 'undefined' && !!performance.mark && performance.measure;
const timingEnabled = hasPerformance && isTimingMode;

export namespace UserTiming {
  function startMarkName(task: Task<any>) {
    return `startTask${task.id}`;
  }
  function endMarkName(task: Task<any>) {
    return `endTask${task.id}`;
  }
  export function markStart(task: Task<any>) {
    if (timingEnabled) performance.mark(startMarkName(task));
  }
  export function markEnd(task: Task<any>) {
    if (timingEnabled) performance.mark(endMarkName(task));
  }
  export function measure(task: Task<any>) {
    if (timingEnabled) performance.measure(`✳️ ${task.name}`, startMarkName(task), endMarkName(task));
  }

  /**
   * `f` run as a measure `name` in timing mode (`setTimingMode`, read each call: set after this module loads): a stage
   * of work in a browser's profile (its Timings track) and in `performance.getEntriesByType('measure')` for benchmarks;
   * `f` alone otherwise
   */
  export function span<T>(name: string, f: () => T): T {
    if (!hasPerformance || !isTimingMode) return f();
    const start = performance.now();
    try {
      return f();
    } finally {
      performance.measure(name, { start, end: performance.now() });
    }
  }
}
