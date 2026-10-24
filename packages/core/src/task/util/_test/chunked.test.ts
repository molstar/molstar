/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { now } from '../../../util/now.js';
import type { RuntimeContext } from '../../execution/runtime-context.js';
import { chunkedSubtask } from '../chunked.js';

describe('chunkedSubtask', () => {
  afterEach(() => {
    vi.useRealTimers();
  });

  it('keeps the gaps between updates near the update interval when the rate changes (perf-audit-2 03/B8)', async () => {
    vi.useFakeTimers({ toFake: ['hrtime', 'performance', 'Date'] });
    const interval = 100;
    let lastUpdate = now();
    const gaps: number[] = [];
    const ctx: RuntimeContext = {
      isSynchronous: false,
      get shouldUpdate() {
        return now() - lastUpdate >= interval;
      },
      update() {
        gaps.push(now() - lastUpdate);
        lastUpdate = now();
      },
    };

    // an item takes 1 ms, 3 ms in every other block of 1000 items
    const state = { done: 0, count: 6000 };
    await chunkedSubtask(
      ctx,
      10,
      state,
      (size, s) => {
        const k = Math.min(size, s.count - s.done);
        for (let i = s.done; i < s.done + k; i++) vi.advanceTimersByTime(Math.floor(i / 1000) % 2 ? 3 : 1);
        s.done += k;
        return k;
      },
      (c) => c.update(),
    );

    expect(state.done).toBe(6000);
    expect(gaps.length).toBeGreaterThan(30);
    // within 1.25× of the interval, 1.5× right after the rate triples (the old sizing reached ~10×)
    for (const g of gaps) expect(g).toBeLessThanOrEqual(1.5 * interval);
    expect(gaps.filter((g) => g > 1.25 * interval).length).toBeLessThanOrEqual(3);
  });
});
