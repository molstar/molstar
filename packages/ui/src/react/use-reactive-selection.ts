/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { listen, type ReadonlyState, type Subscribable } from '@molstar/core/reactive';
import { sameShallow } from '@molstar/core/util/object';
import { useMemo, useRef, useSyncExternalStore } from 'react';

/**
 * What a hook reads: a state (`State`, `derive`, a `BehaviorSubject`: it has a current `value`), or a stream (an
 * `Observable`, ours or rxjs's), or nothing.
 */
export type ReactiveSource<T> = ReadonlyState<T> | Subscribable<T> | null | undefined;

export type ReactiveOptions<T, U = T> = {
  /** Seeds a new source's lifetime (a stream before its first value); it never replaces a value, undefined included. */
  readonly defaultValue?: T;
  /** Whether two selections are the same (no render): shallow by default (react §2). */
  readonly compare?: (left: U, right: U) => boolean;
};

type Callbacks<T, U> = { selector: (value: T) => U; compare: (left: U, right: U) => boolean };

const empty = Symbol('empty');

/** A source with a current value: read it, never wait for (or count) its replay. */
function hasValue<T>(source: ReactiveSource<T>): source is ReadonlyState<T> {
  return source != null && 'value' in source;
}

/**
 * The store `useSyncExternalStore` reads (react §1). A source with a current value is read through `value` on each
 * snapshot and followed through core's `listen`: changes only, so subscribing to a state or a `BehaviorSubject` (ours
 * or rxjs's) never replays its value into the store, and React is told once per change, never on subscribing. A stream
 * without one is subscribed: what it emits while subscribing (a cold `of(x)`) is kept without telling React, which
 * reads the snapshot again after it subscribes; later values tell it.
 */
function createStore<T, U>(source: ReactiveSource<T>, initial: T | typeof empty, callbacks: Callbacks<T, U>) {
  let snapshot = initial;
  let selection: U | typeof empty = empty;
  // The selector the selection was made with: a new one selects again; an equal result keeps the selection's
  // reference, so an inline selector does not loop.
  let selectedWith: ((value: T) => U) | undefined;
  const select = (value: T): boolean => {
    const { selector, compare } = callbacks;
    selectedWith = selector;
    const next = selector(value);
    if (selection !== empty && compare(selection, next)) return false;
    selection = next;
    return true;
  };
  const getSnapshot = (): U | undefined => {
    let value = snapshot;
    if (hasValue(source)) {
      try {
        value = source.value;
      } catch {
        // A current-value stream that ended in an error throws on `value`: its last snapshot stays.
      }
    }
    if (
      value !== empty &&
      (selection === empty || !Object.is(snapshot, value) || selectedWith !== callbacks.selector)
    ) {
      snapshot = value;
      select(value);
    }
    return selection === empty ? undefined : selection;
  };
  const subscribe = (onStoreChange: () => void): (() => void) => {
    if (!source) return () => {};
    if (hasValue(source)) {
      // The next snapshot reads `value`: the listener only says it changed.
      return listen(source, onStoreChange);
    }
    let subscribing = true;
    const subscription = source.subscribe({
      next: (value) => {
        snapshot = value;
        if (select(value) && !subscribing) onStoreChange();
      },
    });
    subscribing = false;
    return () => subscription.unsubscribe();
  };
  return { getSnapshot, subscribe };
}

/**
 * A selection of a source's value, rendering again only when the selection changes (`compare`). The selector and the
 * comparer may be new on each render (their latest is used); the source's identity is its lifetime, so it must be
 * stable across renders: a new one (an inline `.pipe(...)` or `derive(...)`) starts over every render and shows only
 * the default. Memoize it; a source new on renders in a row warns once.
 */
export function useReactiveSelection<T, U>(
  source: ReactiveSource<T>,
  selector: (value: T) => U,
  options?: ReactiveOptions<T, U>,
): U | undefined {
  const sourceRef = useRef(source);
  const changedInARow = useRef(0);
  if (sourceRef.current !== source) {
    sourceRef.current = source;
    changedInARow.current += 1;
    if (changedInARow.current === 2)
      console.warn(
        'useReactive: the source is new on every render (an inline pipe or derive?). Memoize it: a new source starts a new lifetime and shows only its default.',
      );
  } else changedInARow.current = 0;
  const callbacksRef = useRef<Callbacks<T, U> | undefined>(undefined);
  const compare = options?.compare ?? sameShallow;
  if (!callbacksRef.current) callbacksRef.current = { selector, compare };
  else {
    callbacksRef.current.selector = selector;
    callbacksRef.current.compare = compare;
  }
  const callbacks = callbacksRef.current;
  // `defaultValue` only seeds a new source's lifetime, so it is not a dependency.
  const store = useMemo(
    () => createStore(source, options && 'defaultValue' in options ? (options.defaultValue as T) : empty, callbacks),
    [source, callbacks],
  );
  return useSyncExternalStore(store.subscribe, store.getSnapshot, store.getSnapshot);
}
