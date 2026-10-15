/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Basic subset of `Map<K, V>`, only needs to have `get` method */
export type Mapping<K, V> = Pick<Map<K, V>, 'get'>;

/** A map where values are arrays. Handles missing keys when adding values. */
export class MultiMap<K, V> implements Mapping<K, V[]> {
  private _map = new Map<K, V[]>();

  /** Return the array of values assigned to a key (or `undefined` if no such values) */
  get(key: K): V[] | undefined {
    return this._map.get(key);
  }
  /** Append value to a key (handles missing keys) */
  add(key: K, value: V) {
    const values = this._map.get(key);
    if (values) values.push(value);
    else this._map.set(key, [value]);
  }
}
