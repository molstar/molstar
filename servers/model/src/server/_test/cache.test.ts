/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { Cache, type CacheParams } from '../cache.js';

interface Item {
  key: string;
  size: number;
}

const item = (key: string, size = 10): Item => ({ key, size });

function createCache(params: Partial<CacheParams> = {}) {
  const p: CacheParams = { maxApproximateSizeInBytes: 100, entryTimeoutInMs: 1000, ...params };
  return new Cache<Item>(
    (i) => i.key,
    (i) => i.size,
    () => p,
  );
}

describe('model server cache', () => {
  beforeEach(() => {
    vi.useFakeTimers();
    vi.spyOn(console, 'log').mockImplementation(() => {});
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it('drops an entry and its size once it expires', () => {
    const cache = createCache();
    cache.add(item('a', 30));
    expect(cache.stats).toEqual({ count: 1, approximateSize: 30 });

    vi.advanceTimersByTime(1000);
    expect(cache.has('a')).toBe(false);
    expect(cache.get('a')).toBeUndefined();
    expect(cache.stats).toEqual({ count: 0, approximateSize: 0 });
    expect(vi.getTimerCount()).toBe(0);
  });

  it('refreshes the timeout on access without duplicating the entry', () => {
    // molstar#1816: every access appended another list node, and the expiry timer then removed the wrong one.
    const cache = createCache();
    const a = cache.add(item('a', 30));
    for (let i = 0; i < 5; i++) {
      vi.advanceTimersByTime(600);
      expect(cache.get('a')).toBe(a);
    }
    expect(cache.stats).toEqual({ count: 1, approximateSize: 30 });
    expect(vi.getTimerCount()).toBe(1);

    vi.advanceTimersByTime(1000);
    expect(cache.stats).toEqual({ count: 0, approximateSize: 0 });
    expect(vi.getTimerCount()).toBe(0);
  });

  it('keeps one entry when a key is added again', () => {
    const cache = createCache();
    cache.add(item('a', 30));
    const second = cache.add(item('a', 40));
    expect(cache.get('a')).toBe(second);
    expect(cache.stats).toEqual({ count: 1, approximateSize: 40 });
    expect(vi.getTimerCount()).toBe(1);
  });

  it('evicts least recently used entries until the new one fits', () => {
    const cache = createCache();
    cache.add(item('a', 30));
    cache.add(item('b', 30));
    cache.add(item('c', 30));
    cache.get('a'); // b is now the least recently used

    cache.add(item('d', 60));
    expect(cache.has('a')).toBe(true);
    expect(cache.has('b')).toBe(false);
    expect(cache.has('c')).toBe(false);
    expect(cache.has('d')).toBe(true);
    expect(cache.stats).toEqual({ count: 2, approximateSize: 90 });
  });

  it('does not cache an item larger than the whole budget', () => {
    const cache = createCache();
    cache.add(item('a', 30));
    const big = item('big', 101);
    expect(cache.add(big)).toBe(big);
    expect(cache.has('big')).toBe(false);
    expect(cache.has('a')).toBe(true);
    expect(cache.stats).toEqual({ count: 1, approximateSize: 30 });
  });

  it('expires single keys and everything', () => {
    const cache = createCache();
    cache.add(item('a'));
    cache.add(item('b'));
    cache.add(item('c'));

    cache.expire('b');
    expect(cache.has('b')).toBe(false);
    expect(cache.stats).toEqual({ count: 2, approximateSize: 20 });

    cache.expireAll();
    expect(cache.stats).toEqual({ count: 0, approximateSize: 0 });
    expect(vi.getTimerCount()).toBe(0);
  });

  it('reads its parameters at each use', () => {
    const p: CacheParams = { maxApproximateSizeInBytes: 0, entryTimeoutInMs: 1000 };
    const cache = new Cache<Item>(
      (i) => i.key,
      (i) => i.size,
      () => p,
    );
    cache.add(item('a'));
    expect(cache.has('a')).toBe(false);

    p.maxApproximateSizeInBytes = 100;
    cache.add(item('a'));
    expect(cache.has('a')).toBe(true);
  });

  it('does not keep the process alive', () => {
    vi.useRealTimers();
    const cache = createCache({ entryTimeoutInMs: 60_000 });
    cache.add(item('a'));
    const timer = (cache as any).entries.get('a').timeoutId as ReturnType<typeof setTimeout>;
    expect(timer.hasRef()).toBe(false);
    cache.expireAll();
  });
});
