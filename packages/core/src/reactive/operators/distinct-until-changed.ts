/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lift, Subscriber } from '../observable.js';
import type { MonoTypeOperatorFunction } from '../types.js';

function strictEquals(a: unknown, b: unknown) {
  return a === b;
}

class DistinctSubscriber<T, K> extends Subscriber<T> {
  private first = true;
  private last: K | undefined = undefined;
  constructor(
    destination: Subscriber<T>,
    private readonly comparator: (previous: K, current: K) => boolean,
    private readonly keySelector: ((value: T) => K) | undefined,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    let key: K;
    try {
      key = this.keySelector ? this.keySelector(value) : (value as unknown as K);
      if (!this.first && this.comparator(this.last as K, key)) return;
    } catch (e) {
      this.error(e);
      return;
    }
    this.first = false;
    this.last = key;
    this.destination!.next(value);
  }
}

/** Drops a value equal (`comparator`, `===` by default; of `keySelector`'s keys when given) to the last told. */
export function distinctUntilChanged<T>(comparator?: (previous: T, current: T) => boolean): MonoTypeOperatorFunction<T>;
export function distinctUntilChanged<T, K>(
  comparator: (previous: K, current: K) => boolean,
  keySelector: (value: T) => K,
): MonoTypeOperatorFunction<T>;
export function distinctUntilChanged<T, K>(
  comparator: (previous: K, current: K) => boolean = strictEquals,
  keySelector?: (value: T) => K,
): MonoTypeOperatorFunction<T> {
  return (source) =>
    lift(source, (destination: Subscriber<T>) => new DistinctSubscriber(destination, comparator, keySelector));
}
