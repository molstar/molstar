/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { KeyedValue } from './keyed-selector.js';
import { Listened } from './state.js';

/**
 * Values by key, each told apart (view-models.md §3: a row never listens to a whole shared value it cares about one key
 * of): a tree's rows' places by row key, where a reflow far below changes none of the rows on screen. `key(k)` is a
 * state of `k`'s value (undefined when it has none), told only when that value changes (`same`). `replace` takes the
 * whole new map and tells only the listened keys whose value changed: its cost is the keys listened to, not the map.
 */
export class KeyedState<K, V> {
  #values: Map<K, V> = new Map();
  readonly #keys = new Map<K, Set<KeyValueState<K, V>>>();
  readonly #same: (a: V, b: V) => boolean;

  constructor(same: (a: V, b: V) => boolean = Object.is) {
    this.#same = same;
  }

  get(k: K): V | undefined {
    return this.#values.get(k);
  }

  has(k: K): boolean {
    return this.#values.has(k);
  }

  /** All the values now (read-only, and live: change them with `set`, `delete` or `replace`). */
  get values(): ReadonlyMap<K, V> {
    return this.#values;
  }

  /** Sets one key's value; tells it only if the value changed. */
  set(k: K, v: V): void {
    const had = this.#values.has(k);
    if (had && this.#same(this.#values.get(k) as V, v)) return;
    // In place: a set is one key's, nothing copied.
    this.#values.set(k, v);
    this.#tell(k);
  }

  delete(k: K): void {
    if (!this.#values.delete(k)) return;
    this.#tell(k);
  }

  /** Replaces every value at once; tells the listened keys whose value changed (or came, or went). */
  replace(values: ReadonlyMap<K, V>): void {
    const before = this.#values;
    // Taken as it is (the caller's: not changed after).
    this.#values = values as Map<K, V>;
    for (const k of this.#keys.keys()) {
      const had = before.has(k);
      const has = values.has(k);
      if (had !== has || (has && !this.#same(before.get(k) as V, values.get(k) as V))) this.#tell(k);
    }
  }

  /** `k`'s value as a state, told only when it changes (for `watch`, `bind`, `derive`). */
  key(k: K): KeyedValue<V | undefined> {
    return new KeyValueState(this, k);
  }

  #tell(k: K) {
    // `forEach` with a module function: no iterator, in any tier (a `for…of` allocated one before optimization).
    this.#keys.get(k)?.forEach(tellValue);
  }

  /** @internal */
  attach(s: KeyValueState<K, V>) {
    let set = this.#keys.get(s.k);
    if (!set) this.#keys.set(s.k, (set = new Set()));
    set.add(s);
  }

  /** @internal */
  detach(s: KeyValueState<K, V>) {
    const set = this.#keys.get(s.k);
    if (!set?.delete(s)) return;
    if (set.size === 0) this.#keys.delete(s.k);
  }
}

class KeyValueState<K, V> extends Listened<V | undefined> implements KeyedValue<V | undefined> {
  constructor(
    readonly state: KeyedState<K, V>,
    readonly k: K,
  ) {
    super();
  }
  get value(): V | undefined {
    return this.state.get(this.k);
  }
  changed() {
    this.tell(this.value);
  }
  protected override onFirst() {
    this.state.attach(this);
  }
  protected override onLast() {
    this.state.detach(this);
  }
}

const tellValue = (s: { changed(): void }) => s.changed();

/** Values by key, each told apart (`KeyedState`). */
export function keyedState<K, V>(same?: (a: V, b: V) => boolean): KeyedState<K, V> {
  return new KeyedState<K, V>(same);
}
