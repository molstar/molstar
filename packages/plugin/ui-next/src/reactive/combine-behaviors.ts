/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { combineLatest, config, type Observable } from 'rxjs';

import type { Disposable } from './disposable-group.js';
import { shallowEqual } from './shallow-equal.js';
import { observableState } from './observable-state.js';
import type { ReadonlyState } from './readonly-state.js';

type BehaviorSourceMap = Record<string, ReadonlyState<unknown>>;

type SourceValues<T extends BehaviorSourceMap> = {
  [K in keyof T]: T[K] extends ReadonlyState<infer Value> ? Value : never;
};

export type CombinedBehavior<T> = ReadonlyState<T> & Disposable;

export function combineBehaviors<T extends BehaviorSourceMap, U extends Record<string, unknown>>(options: {
  readonly sources: T;
  readonly project: (values: SourceValues<T>) => U;
}): CombinedBehavior<U> {
  const entries = Object.entries(options.sources) as [keyof T & string, ReadonlyState<unknown>][];
  const observables = entries.map(([, source]) => source) as Observable<unknown>[];

  function readValues(rawValues?: readonly unknown[]): SourceValues<T> {
    const values = {} as SourceValues<T>;
    for (let index = 0; index < entries.length; index += 1) {
      const key = entries[index]?.[0];
      if (key === undefined) {
        continue;
      }
      const value = rawValues === undefined ? entries[index]?.[1].value : rawValues[index];
      values[key as keyof SourceValues<T>] = value as SourceValues<T>[keyof SourceValues<T>];
    }
    return values;
  }

  const projected = observableState(options.project(readValues()), { compare: shallowEqual });
  // The projection runs in the subscriber, guarded: one that throws is reported
  // (rxjs' onUnhandledError, else the console; never an uncaught exception) and the value keeps its
  // last good state, then follows the sources again, instead of the chain ending for good.
  const subscription = combineLatest(observables).subscribe((rawValues) => {
    let value: U;
    try {
      value = options.project(readValues(rawValues));
    } catch (error) {
      if (config.onUnhandledError) config.onUnhandledError(error);
      else console.error(error);
      return;
    }
    projected.set(value);
  });
  const asReadonly = projected.asReadonly as CombinedBehavior<U>;
  return Object.assign(asReadonly, {
    dispose() {
      subscription.unsubscribe();
    },
  });
}
