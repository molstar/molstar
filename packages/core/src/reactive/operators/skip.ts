/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lift, Subscriber } from '../observable.js';
import type { MonoTypeOperatorFunction } from '../types.js';

class SkipSubscriber<T> extends Subscriber<T> {
  constructor(
    destination: Subscriber<T>,
    private remaining: number,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    if (this.remaining > 0) this.remaining--;
    else this.destination!.next(value);
  }
}

/** Drops the first `count` values. */
export function skip<T>(count: number): MonoTypeOperatorFunction<T> {
  return (source) =>
    count <= 0 ? source : lift(source, (destination: Subscriber<T>) => new SkipSubscriber(destination, count));
}
