/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subscriber } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import { AsyncQueue } from '../async-queue.js';

describe('AsyncQueue', () => {
  it('runs values in order; a removed value resolves false', async () => {
    const q = new AsyncQueue<string>();
    const order: string[] = [];
    const run = async (v: string) => {
      const go = await q.enqueue(v);
      order.push(`${v}:${go}`);
      if (go) q.handled(v);
    };
    const all = [run('a'), run('b'), run('c'), run('d')];
    expect(q.remove('c')).toBe(true);
    await Promise.all(all);
    expect(order).toEqual(['a:true', 'c:false', 'b:true', 'd:true']);
    expect(q.length).toBe(0);
  });

  it('a dequeue resolves only the new head (perf-audit-2 03/P6)', async () => {
    const next = vi.spyOn(Subscriber.prototype, 'next');
    try {
      const q = new AsyncQueue<number>();
      const k = 50;
      const done: number[] = [];
      const waits = Array.from({ length: k }, (_, i) =>
        Promise.resolve(q.enqueue(i)).then((go) => {
          done.push(i);
          if (go) q.handled(i);
        }),
      );
      await Promise.all(waits);
      expect(done).toEqual(Array.from({ length: k }, (_, i) => i));
      // a broadcast per dequeue called every waiter: k²/2 calls
      expect(next.mock.calls.length).toBeLessThan(k);
    } finally {
      next.mockRestore();
    }
  });
});
