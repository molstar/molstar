/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { batch, Derived, signalChange } from './state.js';
import { Observable } from './observable.js';
import type { ReadonlyState, Subscribable, Unsubscribable } from './types.js';

// Keyed states (reactive §8): one value of many, told apart. A key's state is a node of the state graph (`derive`
// follows it as one), made by `key(k)`, registered while followed and dropped with its last follower. A change visits
// the nodes of the keys it concerns, never every key listened.

/** The followed nodes of each key: one, or an array when several `key(k)` nodes of one key are followed. */
class KeyRegistry<K, N> {
  private readonly _byKey = new Map<K, N | N[]>();

  get size(): number {
    return this._byKey.size;
  }

  keys(): IterableIterator<K> {
    return this._byKey.keys();
  }

  get(k: K): N | N[] | undefined {
    return this._byKey.get(k);
  }

  add(k: K, node: N): void {
    const at = this._byKey.get(k);
    if (at === undefined) this._byKey.set(k, node);
    else if (Array.isArray(at)) at.push(node);
    else this._byKey.set(k, [at, node]);
  }

  delete(k: K, node: N): void {
    const at = this._byKey.get(k);
    if (at === node) this._byKey.delete(k);
    else if (Array.isArray(at)) {
      const i = at.indexOf(node);
      if (i >= 0) at.splice(i, 1);
      if (at.length === 1) this._byKey.set(k, at[0]);
    }
  }
}

/**
 * Marks a key's followed nodes. The global version moves even with none: a cold value derived from the key reads it
 * again (`signalChange` per call: callers batch).
 */
function signalKey<K>(registry: KeyRegistry<K, Derived<any>>, k: K): void {
  signalChange(registry.get(k) as Derived<any> | Derived<any>[] | undefined);
}

class KeyValue<K, V> extends Derived<V | undefined> {
  constructor(
    private readonly _owner: KeyedState<K, V>,
    /** @internal */
    readonly _key: K,
  ) {
    super([], Object.is);
    // No sources: its value is read from the owner, which tells it of changes; unfollowed, read again each time.
    this._hasForeign = true;
  }

  protected _updateStamps(): boolean {
    return true;
  }

  protected _compute(): V | undefined {
    return this._owner.get(this._key);
  }

  protected override _hot(): void {
    super._hot();
    this._owner._registry.add(this._key, this);
  }

  protected override _cold(): void {
    this._owner._registry.delete(this._key, this);
    super._cold();
  }
}

/**
 * Values by key, each a state of its own (`key(k)`): a set tells only that key's followers, and `replace` tells the
 * followed keys whose value changed, so its cost is the keys followed, not the map. `equals` (`Object.is` by default)
 * decides a change.
 */
export class KeyedState<K, V> {
  private _values: Map<K, V> = new Map();
  /** @internal */
  readonly _registry = new KeyRegistry<K, Derived<any>>();

  constructor(private readonly _equals: (a: V, b: V) => boolean = Object.is) {}

  get(k: K): V | undefined {
    return this._values.get(k);
  }

  has(k: K): boolean {
    return this._values.has(k);
  }

  /** All the values now (read-only and live: change them with `set`, `delete` or `replace`). */
  get values(): ReadonlyMap<K, V> {
    return this._values;
  }

  /** Sets one key's value; tells its followers only if the value changed. */
  set(k: K, v: V): void {
    const values = this._values;
    if (values.has(k) && this._equals(values.get(k) as V, v)) return;
    values.set(k, v);
    signalKey(this._registry, k);
  }

  delete(k: K): void {
    if (!this._values.delete(k)) return;
    signalKey(this._registry, k);
  }

  /** Replaces every value (`values` is taken as is: not changed after); tells the followed keys that changed. */
  replace(values: ReadonlyMap<K, V>): void {
    const before = this._values;
    this._values = values as Map<K, V>;
    const registry = this._registry;
    if (registry.size === 0) {
      signalChange(undefined);
      return;
    }
    batch(() => {
      for (const k of registry.keys()) {
        const had = before.has(k);
        const has = values.has(k);
        if (had !== has || (has && !this._equals(before.get(k) as V, values.get(k) as V))) signalKey(registry, k);
      }
      // Unfollowed keys are read again on their next read (the global version moved).
      signalChange(undefined);
    });
  }

  /** `k`'s value as a state (undefined while it has none), told only when it changes. */
  key(k: K): ReadonlyState<V | undefined> {
    return new KeyValue(this, k);
  }
}

class KeyFlag<K> extends Derived<boolean> {
  constructor(
    private readonly _selector: KeyedSelector<K>,
    /** @internal */
    readonly _key: K,
  ) {
    super([], Object.is);
    this._hasForeign = true;
  }

  protected _updateStamps(): boolean {
    return true;
  }

  protected _compute(): boolean {
    return this._selector.isCurrent(this._key);
  }

  protected override _hot(): void {
    super._hot();
    this._selector._attach(this);
  }

  protected override _cold(): void {
    this._selector._detach(this);
    super._cold();
  }
}

/**
 * Whether a key is the current value of `source` (the picked row of thousands): `key(k)` is a state of whether `k` is
 * current, and a change of the current value tells only the two keys it concerns (reactive §8). It follows `source`
 * only while a key is followed. With a source that has a `value` (a state, a `BehaviorSubject`), a key reads it, so a
 * key and anything else derived from the source agree within one change.
 */
export class KeyedSelector<K> {
  /** @internal */
  readonly _registry = new KeyRegistry<K, Derived<any>>();
  private _followed = 0;
  private _current: K | undefined = undefined;
  private _hasCurrent = false;
  private _stop: Unsubscribable | undefined = undefined;
  private readonly _hasValue: boolean;

  constructor(private readonly _source: ReadonlyState<K> | Observable<K> | Subscribable<K>) {
    this._hasValue = 'value' in _source;
  }

  /** Whether `k` is the current value now. */
  isCurrent(k: K): boolean {
    if (this._hasValue) return Object.is((this._source as ReadonlyState<K>).value, k);
    return this._hasCurrent && Object.is(this._current, k);
  }

  /** Whether `k` is current, as a state told only when that changes. */
  key(k: K): ReadonlyState<boolean> {
    return new KeyFlag(this, k);
  }

  /** The same as `key(k)` (a state is a stream that replays its value). */
  is(k: K): ReadonlyState<boolean> {
    return this.key(k);
  }

  private readonly _changed = (v: K) => {
    const had = this._hasCurrent;
    const previous = this._current;
    if (had && Object.is(previous, v)) return;
    this._current = v;
    this._hasCurrent = true;
    const registry = this._registry;
    batch(() => {
      if (had) signalKey(registry, previous as K);
      signalKey(registry, v);
    });
  };

  /** @internal */
  _attach(flag: KeyFlag<K>): void {
    this._registry.add(flag._key, flag);
    if (this._followed++ > 0) return;
    const source = this._source;
    if (this._hasValue) {
      this._current = (source as ReadonlyState<K>).value;
      this._hasCurrent = true;
    }
    if (typeof (source as ReadonlyState<K>).listen === 'function') {
      const stop = (source as ReadonlyState<K>).listen(this._changed);
      this._stop = { unsubscribe: stop };
    } else if (source instanceof Observable) {
      this._stop = source.subscribe(this._changed);
    } else {
      this._stop = (source as Subscribable<K>).subscribe({ next: this._changed });
    }
  }

  /** @internal */
  _detach(flag: KeyFlag<K>): void {
    this._registry.delete(flag._key, flag);
    if (--this._followed > 0) return;
    const stop = this._stop;
    this._stop = undefined;
    this._hasCurrent = false;
    this._current = undefined;
    stop?.unsubscribe();
  }
}
