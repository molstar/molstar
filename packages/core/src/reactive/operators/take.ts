/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EMPTY } from '../creation.js';
import { lift, Subscriber } from '../observable.js';
import type { MonoTypeOperatorFunction } from '../types.js';

class TakeSubscriber<T> extends Subscriber<T> {
  constructor(
    destination: Subscriber<T>,
    private remaining: number,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    if (--this.remaining < 0) return;
    const destination = this.destination!;
    destination.next(value);
    if (this.remaining === 0) this.complete();
  }
}

/** The first `count` values, then complete (and unsubscribe from the source). */
export function take<T>(count: number): MonoTypeOperatorFunction<T> {
  return (source) =>
    count <= 0 ? EMPTY : lift(source, (destination: Subscriber<T>) => new TakeSubscriber(destination, count));
}
