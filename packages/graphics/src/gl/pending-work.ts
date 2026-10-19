/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Subject } from 'rxjs';

/** Asynchronous work a frame may depend on, e.g. images loading into textures */
export interface PendingWork {
  /** Tracks `promise` until it settles; a rejection is reported with `label` */
  track<T>(promise: Promise<T>, label: string): Promise<T>;
  /** Whether nothing tracked is pending */
  isReady(): boolean;
  /**
   * Number of work items tracked so far; when it changes during a render, the render started work it may have skipped
   * draws for (e.g. a pipeline compiling), even if that work has settled since
   */
  readonly trackedCount: number;
  /** Labels of the work tracked since `trackedCount` was `count` (the most recent ones, for diagnostics) */
  trackedSince(count: number): string[];
  /** Resolves when nothing tracked is pending, including work tracked while waiting */
  whenReady(): Promise<void>;
  /** Emits when tracked work settles, e.g. to request a redraw */
  readonly settled: Subject<void>;
  /** Forgets the tracked work (context loss, destroy) and releases waiters */
  clear(): void;
}

const RecentLabelCount = 64;

export function createPendingWork(): PendingWork {
  const pending = new Map<Promise<unknown>, string>();
  const settled = new Subject<void>();
  let waiters: (() => void)[] = [];
  let trackedCount = 0;
  /** labels of the last tracked work, `recent[i]` tracked as number `trackedCount - recent.length + i` */
  const recent: string[] = [];

  const release = () => {
    const w = waiters;
    waiters = [];
    for (const resolve of w) resolve();
  };

  return {
    track: <T>(promise: Promise<T>, label: string) => {
      pending.set(promise, label);
      trackedCount += 1;
      recent.push(label);
      if (recent.length > RecentLabelCount) recent.shift();
      const done = () => {
        if (!pending.delete(promise)) return;
        settled.next();
        if (pending.size === 0) release();
      };
      promise.then(done, (e) => {
        // not when forgotten (`clear`): failing with the context, e.g. reads of a released device
        if (pending.has(promise)) console.error(`'${label}' failed`, e);
        done();
      });
      return promise;
    },
    isReady: () => pending.size === 0,
    get trackedCount() {
      return trackedCount;
    },
    trackedSince: (count: number) => recent.slice(Math.max(0, recent.length - (trackedCount - count))),
    whenReady: () => (pending.size === 0 ? Promise.resolve() : new Promise<void>((resolve) => waiters.push(resolve))),
    settled,
    clear: () => {
      pending.clear();
      release();
    },
  };
}
