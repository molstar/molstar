/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Observable } from '../observable.js';
import type { OperatorFunction } from '../types.js';

/** `values`, then the source's. */
export function startWith<T, D extends readonly unknown[]>(...values: D): OperatorFunction<T, T | D[number]> {
  return (source) =>
    new Observable<T | D[number]>((subscriber) => {
      for (let i = 0; i < values.length; i++) {
        if (subscriber.closed) return;
        subscriber.next(values[i]);
      }
      return source.subscribe(subscriber);
    });
}
