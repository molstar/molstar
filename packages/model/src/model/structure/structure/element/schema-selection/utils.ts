/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** A map where values are arrays. Handles missing keys when adding values. */
export class MultiMap<K, V> implements Mapping<K, V[]> {
  private _map = new Map();

  /** Return the array of values assidned to a key (or `undefined` if no such values) */
  get(key: K): V[] | undefined {
    return this._map.get(key);
  }
  /** Append value to a key (handles missing keys) */
  add(key: K, value: V) {
    if (!this._map.has(key)) {
      this._map.set(key, []);
    }
    this._map.get(key)!.push(value);
  }
}

/** Basic subset of `Map<K, V>`, only needs to have `get` method */
export type Mapping<K, V> = Pick<Map<K, V>, 'get'>;

/** Implementation of `Map` where keys are integers
 * and most keys are expected to be from interval `[0, limit)`.
 * For the keys within this interval, performance is better than `Map` (implemented by array).
 * For the keys out of this interval, performance is slightly worse than `Map`. */
export class NumberMap<K extends number, V> implements Mapping<K, V> {
  private array: V[];
  private map: Map<K, V>;
  constructor(public readonly limit: K) {
    this.array = new Array(limit);
    this.map = new Map();
  }
  get(key: K): V | undefined {
    if (0 <= key && key < this.limit) return this.array[key];
    else return this.map.get(key);
  }
  set(key: K, value: V): void {
    if (0 <= key && key < this.limit) this.array[key] = value;
    else this.map.set(key, value);
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
