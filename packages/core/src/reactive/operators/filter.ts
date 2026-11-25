/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lift, Subscriber } from '../observable.js';
import type { MonoTypeOperatorFunction, OperatorFunction } from '../types.js';

type Falsy = false | 0 | -0 | 0n | '' | null | undefined;

class FilterSubscriber<T> extends Subscriber<T> {
  private index = 0;
  constructor(
    destination: Subscriber<T>,
    private readonly predicate: (value: T, index: number) => unknown,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    let pass: unknown;
    try {
      pass = this.predicate(value, this.index++);
    } catch (e) {
      this.error(e);
      return;
    }
    if (pass) this.destination!.next(value);
  }
}

/** The values `predicate` passes. */
export function filter<T, S extends T>(predicate: (value: T, index: number) => value is S): OperatorFunction<T, S>;
export function filter<T>(predicate: BooleanConstructor): OperatorFunction<T, Exclude<T, Falsy>>;
export function filter<T>(predicate: (value: T, index: number) => boolean): MonoTypeOperatorFunction<T>;
export function filter<T>(predicate: (value: T, index: number) => unknown): MonoTypeOperatorFunction<T> {
  return (source) => lift(source, (destination: Subscriber<T>) => new FilterSubscriber(destination, predicate));
}
