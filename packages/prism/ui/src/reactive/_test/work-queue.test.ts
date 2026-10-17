/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { createRoot, flush, requireScope, watch } from '../../dom/scope.js';
import { workQueue } from '../../dom/scoped.js';
import { flushFrame } from '../frame-state.js';
import { observableState } from '../observable-state.js';
import { createLatest, createWorkQueue } from '../work-queue.js';

// Work wanted by key, run once in its lane (rendering.md §2.7).

describe('workQueue', () => {
  it('the last want of a key wins: a burst runs once, the latest', () => {
    const q = createWorkQueue({ lane: 'frame' });
    const ran: number[] = [];
    for (let i = 0; i < 10; i++) q.want('k', () => ran.push(i));
    expect(ran).toEqual([]);
    flushFrame();
    expect(ran).toEqual([9]);
    expect(q.size).toBe(0);
  });

  it('the task lane runs in the runtime flush, before the watches it feeds render', () => {
    const q = createWorkQueue({ lane: 'task' });
    const rows = observableState(0);
    const seen: number[] = [];
    createRoot(() => watch([rows.asReadonly], () => seen.push(rows.value)));
    for (let i = 1; i <= 5; i++) q.want('flatten', () => rows.set(i));
    flush();
    // One render, of the value the one run computed.
    expect(seen).toEqual([0, 5]);
  });

  it('the task lane runs in its own microtask without a flush', async () => {
    const q = createWorkQueue({ lane: 'task' });
    let n = 0;
    q.want('a', () => n++);
    q.want('a', () => n++);
    await Promise.resolve();
    await Promise.resolve();
    expect(n).toBe(1);
  });

  it('runs by priority, asked at each run; stillWanted drops what is no longer wanted', () => {
    const q = createWorkQueue({ lane: 'frame' });
    const order: string[] = [];
    let far = 1;
    q.want('near', () => order.push('near'), { priority: 5 });
    q.want('far', () => order.push('far'), { priority: () => far });
    q.want('gone', () => order.push('gone'), { priority: 9, stillWanted: () => false });
    far = 10;
    flushFrame();
    expect(order).toEqual(['far', 'near']);
  });

  it('a budget carries the rest to the next run, the highest first', () => {
    const q = createWorkQueue({ lane: 'frame', budgetMs: 1 });
    const ran: number[] = [];
    const spin = () => {
      const t = performance.now();
      while (performance.now() - t < 2) {}
    };
    for (let i = 0; i < 5; i++)
      q.want(
        i,
        () => {
          spin();
          ran.push(i);
        },
        { priority: i },
      );
    flushFrame();
    // At least one a run, however slow.
    expect(ran).toEqual([4]);
    flushFrame();
    expect(ran).toEqual([4, 3]);
    q.flush();
    expect(ran).toEqual([4, 3, 2, 1, 0]);
  });

  it('cancel drops a key; work wanted while running waits for the next run', () => {
    const q = createWorkQueue({ lane: 'frame' });
    const ran: string[] = [];
    q.want('a', () => {
      ran.push('a');
      q.want('a', () => ran.push('a again'));
    });
    q.want('b', () => ran.push('b'));
    q.cancel('b');
    flushFrame();
    expect(ran).toEqual(['a']);
    flushFrame();
    expect(ran).toEqual(['a', 'a again']);
  });

  it('a scoped queue pauses with its scope (what waited runs on resume) and is dropped when it ends', () => {
    let q!: ReturnType<typeof workQueue>;
    let setActive!: (a: boolean) => void;
    const dispose = createRoot((d) => {
      q = workQueue({ lane: 'frame' });
      const scope = requireScope('test');
      setActive = (a) => scope.setActive(a);
      return d;
    });
    const run = vi.fn();
    setActive(false);
    q.want('k', run);
    flushFrame();
    expect(run).not.toHaveBeenCalled();
    setActive(true);
    flushFrame();
    expect(run).toHaveBeenCalledTimes(1);
    q.want('k', run);
    dispose();
    flushFrame();
    expect(run).toHaveBeenCalledTimes(1);
  });
});

describe('latest', () => {
  it('a newer call aborts the one before: the stale answer resolves to undefined', async () => {
    const seen: AbortSignal[] = [];
    const resolvers: ((v: string) => void)[] = [];
    const search = createLatest((signal, q: string) => {
      seen.push(signal);
      return new Promise<string>((r) => resolvers.push(() => r(q)));
    });
    const a = search('a');
    const b = search('ab');
    expect(seen[0].aborted).toBe(true);
    expect(seen[1].aborted).toBe(false);
    resolvers[1]('');
    resolvers[0]('');
    expect(await b).toBe('ab');
    expect(await a).toBeUndefined();
  });

  it('abort() supersedes the call in flight', async () => {
    const run = createLatest(async (signal) => {
      await Promise.resolve();
      if (signal.aborted) throw new Error('aborted');
      return 1;
    });
    const p = run();
    run.abort();
    expect(await p).toBeUndefined();
  });
});
