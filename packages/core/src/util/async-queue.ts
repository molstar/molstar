/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arrayRemoveInPlace } from './array.js';

export class AsyncQueue<T> {
  private queue: T[] = [];
  /** The waiters of each queued value: a dequeue resolves only the new head's, not a broadcast to all (perf-audit-2 03/P6) */
  private waiters = new Map<T, ((stillPresent: boolean) => void)[]>();

  get length() {
    return this.queue.length;
  }

  enqueue(v: T) {
    this.queue.push(v);
    if (this.queue.length === 1) return true;
    return this.waitFor(v);
  }

  handled(v: T) {
    arrayRemoveInPlace(this.queue, v);
    if (this.queue.length > 0) {
      this.resolve(this.queue[0], true);
    }
  }

  remove(v: T) {
    const rem = arrayRemoveInPlace(this.queue, v);
    if (rem) this.resolve(v, false);
    return rem;
  }

  private resolve(v: T, stillPresent: boolean) {
    const waiting = this.waiters.get(v);
    if (!waiting) return;
    this.waiters.delete(v);
    for (const res of waiting) res(stillPresent);
  }

  private waitFor(t: T): Promise<boolean> {
    return new Promise((res) => {
      const waiting = this.waiters.get(t);
      if (waiting) waiting.push(res);
      else this.waiters.set(t, [res]);
    });
  }
}
