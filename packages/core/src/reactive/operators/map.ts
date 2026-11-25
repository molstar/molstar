/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lift, Subscriber } from '../observable.js';
import type { OperatorFunction } from '../types.js';

class MapSubscriber<T, R> extends Subscriber<T> {
  private index = 0;
  constructor(
    destination: Subscriber<R>,
    private readonly project: (value: T, index: number) => R,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    let result: R;
    try {
      result = this.project(value, this.index++);
    } catch (e) {
      this.error(e);
      return;
    }
    this.destination!.next(result);
  }
}

/** Each value through `project`; what it throws is the stream's error. */
export function map<T, R>(project: (value: T, index: number) => R): OperatorFunction<T, R> {
  return (source) => lift(source, (destination: Subscriber<R>) => new MapSubscriber(destination, project));
}
