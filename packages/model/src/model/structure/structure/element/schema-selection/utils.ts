/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Mapping } from '@molstar/core/util/multi-map';

/** Implementation of `Map` where keys are integers and the values are non-negative integers (so they can be kept in a typed array)
 * and most keys are expected to be from interval `[offset, offset + limit)`.
 * For the keys within this interval, performance is better than `Map` (implemented by array, 4 bytes per key).
 * For the keys out of this interval, performance is slightly worse than `Map`.
 * Keys are numbers: a key of another type has no value (it would otherwise be converted to a number by the array). */
export class NumberMap<K extends number, V extends number> implements Mapping<K, V> {
  private array: Int32Array;
  private map: Map<K, V>;
  constructor(
    public readonly limit: number,
    public readonly offset: number = 0,
  ) {
    this.array = new Int32Array(limit).fill(-1);
    this.map = new Map();
  }
  get(key: K): V | undefined {
    if (typeof key !== 'number') return undefined;
    const i = key - this.offset;
    if (0 <= i && i < this.limit) {
      const value = this.array[i];
      return value < 0 ? undefined : (value as V);
    }
    return this.map.get(key);
  }
  set(key: K, value: V): void {
    const i = key - this.offset;
    if (0 <= i && i < this.limit) this.array[i] = value;
    else this.map.set(key, value);
  }
  /** The number of keys which are not in the interval of the array */
  get overflowSize() {
    return this.map.size;
  }
}

/** Return `true` if `value` is not `undefined` or `null`.
 * Prefer this over `value !== undefined`
 * (for maybe if we want to allow `null` in `AnnotationRow` in the future) */
export function isDefined<T>(value: T | undefined | null): value is T {
  return value !== undefined && value !== null;
}
/** Return `true` if at least one of `values` is not `undefined` or `null`. */
export function isAnyDefined(...values: any[]): boolean {
  return values.some(isDefined);
}
