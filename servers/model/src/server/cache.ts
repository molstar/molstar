/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { ConsoleLogger } from '@molstar/core/util/console-logger';

interface CacheEntry<T> {
  key: string;
  approximateSize: number;
  timeoutId: ReturnType<typeof setTimeout> | undefined;
  item: T;
}

export interface CacheParams {
  /** Budget for the summed `approximateSize` of the entries; 0 turns the cache off. */
  maxApproximateSizeInBytes: number;
  /** An entry is dropped when it has not been accessed for this long. */
  entryTimeoutInMs: number;
}

export interface CacheStats {
  count: number;
  approximateSize: number;
}

/**
 * LRU cache with a size budget and an idle timeout per entry.
 *
 * The entries live in a `Map`, whose iteration order is insertion order: an access re-inserts its entry, so the
 * first entry is always the least recently used one and eviction takes from the front.
 */
export class Cache<T> {
  private entries = new Map<string, CacheEntry<T>>();
  private approximateSize = 0;

  private schedule(e: CacheEntry<T>) {
    if (e.timeoutId !== void 0) clearTimeout(e.timeoutId);
    e.timeoutId = setTimeout(() => this.remove(e, 'expired'), this.params().entryTimeoutInMs);
    // cached entries must not keep the process alive (local CLI exits when its jobs are done)
    e.timeoutId.unref?.();
  }

  private remove(e: CacheEntry<T>, reason?: 'expired' | 'evicted') {
    // timer may fire for an entry that was replaced in the meantime
    if (this.entries.get(e.key) !== e) return;
    if (e.timeoutId !== void 0) clearTimeout(e.timeoutId);
    e.timeoutId = void 0;
    this.entries.delete(e.key);
    this.approximateSize -= e.approximateSize;
    if (reason) ConsoleLogger.log('Cache', `${e.key} ${reason} (${this.formatStats()}).`);
  }

  private formatStats() {
    return `${this.entries.size} entries, ~${Math.round(this.approximateSize / (1024 * 1024))} MB`;
  }

  get stats(): CacheStats {
    return { count: this.entries.size, approximateSize: this.approximateSize };
  }

  expireAll() {
    for (const e of Array.from(this.entries.values())) this.remove(e);
  }

  expire(key: string) {
    const e = this.entries.get(key);
    if (e) this.remove(e, 'expired');
  }

  /** Adds `item`, evicting the least recently used entries until it fits. An item larger than the whole budget is not cached. */
  add(item: T) {
    const key = this.keyGetter(item);
    const approximateSize = this.sizeGetter(item);

    const existing = this.entries.get(key);
    if (existing) this.remove(existing);

    const { maxApproximateSizeInBytes } = this.params();
    if (approximateSize > maxApproximateSizeInBytes) {
      ConsoleLogger.log('Cache', `${key} not cached: larger than the cache (${this.formatStats()}).`);
      return item;
    }

    for (const e of Array.from(this.entries.values())) {
      if (this.approximateSize + approximateSize <= maxApproximateSizeInBytes) break;
      this.remove(e, 'evicted');
    }

    const entry: CacheEntry<T> = { key, approximateSize, timeoutId: void 0, item };
    this.entries.set(key, entry);
    this.approximateSize += approximateSize;
    this.schedule(entry);
    ConsoleLogger.log('Cache', `${key} added (${this.formatStats()}).`);

    return item;
  }

  has(key: string) {
    return this.entries.has(key);
  }

  get(key: string) {
    const e = this.entries.get(key);
    if (!e) return void 0;
    // re-insert: the entry becomes the most recent one
    this.entries.delete(key);
    this.entries.set(key, e);
    this.schedule(e);
    ConsoleLogger.log('Cache', `${key} accessed.`);
    return e.item;
  }

  constructor(
    private keyGetter: (i: T) => string,
    private sizeGetter: (i: T) => number,
    /** Read on every use: the server configures itself after this module is loaded. */
    private params: () => CacheParams,
  ) {}
}
