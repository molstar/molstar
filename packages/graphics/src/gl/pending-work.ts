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
  /** Resolves when nothing tracked is pending, including work tracked while waiting */
  whenReady(): Promise<void>;
  /** Emits when tracked work settles, e.g. to request a redraw */
  readonly settled: Subject<void>;
  /** Forgets the tracked work (context loss, destroy) and releases waiters */
  clear(): void;
}

export function createPendingWork(): PendingWork {
  const pending = new Map<Promise<unknown>, string>();
  const settled = new Subject<void>();
  let waiters: (() => void)[] = [];
  let trackedCount = 0;

  const release = () => {
    const w = waiters;
    waiters = [];
    for (const resolve of w) resolve();
  };

  return {
    track: <T>(promise: Promise<T>, label: string) => {
      pending.set(promise, label);
      trackedCount += 1;
      const done = () => {
        if (!pending.delete(promise)) return;
        settled.next();
        if (pending.size === 0) release();
      };
      promise.then(done, (e) => {
        console.error(`'${label}' failed`, e);
        done();
      });
      return promise;
    },
    isReady: () => pending.size === 0,
    get trackedCount() {
      return trackedCount;
    },
    whenReady: () => (pending.size === 0 ? Promise.resolve() : new Promise<void>((resolve) => waiters.push(resolve))),
    settled,
    clear: () => {
      pending.clear();
      release();
    },
  };
}
