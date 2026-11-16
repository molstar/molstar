/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isState, type ReadonlyState, RecordState, State, toState } from '@molstar/core/reactive';
import { isPlainRecord, sameShallow } from '@molstar/core/util/object';
import { type MaybeReactive, onCleanup, toStream } from './dom/index.js';

/**
 * A view's state (core's `State`) compared shallowly: a set of an equal copy (a record or an array built anew) tells
 * nothing. A plain record's state is a `RecordState` (`update(partial)`).
 */
export function state<T extends Record<string, unknown>>(initial: T, equals?: (a: T, b: T) => boolean): RecordState<T>;
export function state<T>(initial: T, equals?: (a: T, b: T) => boolean): State<T>;
export function state<T>(initial: T, equals: (a: T, b: T) => boolean = sameShallow): State<T> {
  return isPlainRecord(initial)
    ? (new RecordState(initial as Record<string, unknown>, equals as never) as unknown as State<T>)
    : new State(initial, equals);
}

/**
 * A prop that may be a plain value, a state or a stream, as a state: a state as it is; a stream fed into a state
 * (ended with the current scope); a plain value (or `fallback` for none) as a state of its own.
 */
export function held<T>(value: MaybeReactive<T> | undefined, fallback: T): ReadonlyState<T> {
  if (isState(value)) return value as ReadonlyState<T>;
  const stream = toStream(value as MaybeReactive<T>);
  if (stream) {
    const s = toState(stream, fallback);
    onCleanup(() => s.dispose());
    return s;
  }
  return new State(value === undefined ? fallback : (value as T));
}
