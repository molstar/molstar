/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable, Subscription } from 'rxjs';
import type { ReadonlyState } from './readonly-state.js';
import { Listened } from './state.js';

/** A keyed value: a state (`value`, `listen`, an Observable when subscribed) told only when its key's value changes. */
export type KeyedValue<T> = ReadonlyState<T> & { listen(fn: (value: T) => void): () => void };

/**
 * Whether a key is the current one (the picked row of thousands; view-models.md §3: a row never listens to a whole
 * shared value it cares about one key of). `key(k)` is a state of whether `k` is current, told only when `k` becomes
 * or stops being it: a change of the current value tells the two keys it concerns, not every row. Follows its source
 * by `listen` when it is a state (no rxjs on that path), else by one subscription; either only while a key is listened.
 */
export class KeyedSelector<K> {
  readonly #source: Observable<K> | ReadonlyState<K>;
  readonly #flags = new Map<K, Set<KeyFlag<K>>>();
  #listened = 0;
  #current: K | undefined;
  #has = false;
  #stop: (() => void) | Subscription | undefined;

  constructor(source: Observable<K> | ReadonlyState<K>) {
    this.#source = source;
  }

  /** Whether `k` is the current value now. */
  isCurrent(k: K): boolean {
    if (this.#has) return Object.is(this.#current, k);
    const s = this.#source as Partial<ReadonlyState<K>>;
    return 'value' in s ? Object.is(s.value, k) : false;
  }

  /** Whether `k` is current, as a state told only when that changes (for `watch`, `bind`, `derive`). */
  key(k: K): KeyedValue<boolean> {
    return new KeyFlag(this, k);
  }

  /** The same as a stream (it replays whether `k` is current). */
  is(k: K): Observable<boolean> {
    return this.key(k);
  }

  readonly #changed = (v: K) => {
    if (this.#has && Object.is(v, this.#current)) return;
    const previous = this.#current;
    const had = this.#has;
    this.#current = v;
    this.#has = true;
    if (had) this.#tell(previous as K);
    this.#tell(v);
  };

  #tell(k: K) {
    // A Set iterates safely while one is removed (a flag that stops listening in its own call); no copy a change, and
    // `forEach` with a module function: no iterator, in any tier.
    this.#flags.get(k)?.forEach(tellFlag);
  }

  /** @internal */
  attach(flag: KeyFlag<K>) {
    let set = this.#flags.get(flag.k);
    if (!set) this.#flags.set(flag.k, (set = new Set()));
    set.add(flag);
    if (this.#listened++ > 0) return;
    const s = this.#source as Partial<ReadonlyState<K>> & Observable<K>;
    if (typeof s.listen === 'function' && 'value' in s) {
      this.#current = s.value as K;
      this.#has = true;
      this.#stop = s.listen(this.#changed);
    } else this.#stop = s.subscribe(this.#changed);
  }

  /** @internal */
  detach(flag: KeyFlag<K>) {
    const set = this.#flags.get(flag.k);
    if (!set?.delete(flag)) return;
    if (set.size === 0) this.#flags.delete(flag.k);
    if (--this.#listened > 0) return;
    const stop = this.#stop;
    this.#stop = undefined;
    this.#has = false;
    this.#current = undefined;
    if (typeof stop === 'function') stop();
    else stop?.unsubscribe();
  }
}

const tellFlag = (f: { changed(): void }) => f.changed();

class KeyFlag<K> extends Listened<boolean> implements KeyedValue<boolean> {
  #last = false;
  constructor(
    readonly selector: KeyedSelector<K>,
    readonly k: K,
  ) {
    super();
  }
  get value(): boolean {
    return this.selector.isCurrent(this.k);
  }
  changed() {
    const v = this.value;
    if (v === this.#last) return;
    this.#last = v;
    this.tell(v);
  }
  protected override onFirst() {
    this.selector.attach(this);
    this.#last = this.value;
  }
  protected override onLast() {
    this.selector.detach(this);
  }
}

/** A `KeyedSelector` over `source`. */
export function keyedSelector<K>(source: Observable<K> | ReadonlyState<K>): KeyedSelector<K> {
  return new KeyedSelector(source);
}
