/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useMemo, useRef, useSyncExternalStore } from 'react';
import { type Observable, skip } from 'rxjs';

import { sameShallow } from '../reactive/same-shallow.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

export type ObservableOptions<T, U = T> = {
  /** Seeds a new source lifetime; it does not replace an emitted value, including undefined. */
  readonly defaultValue?: T;
  readonly compare?: (left: U, right: U) => boolean;
};

type ObservableCallbacks<T, U> = {
  selector: (value: T) => U;
  compare: (left: U, right: U) => boolean;
};

const empty = Symbol('empty');

function hasValue<T>(source: Observable<T> | null | undefined): source is ReadonlyState<T> {
  return source != null && 'value' in source;
}

function createStore<T, U>(
  source: Observable<T> | null | undefined,
  initial: T | typeof empty,
  callbacks: ObservableCallbacks<T, U>,
) {
  let snapshot = initial;
  let selection: U | typeof empty = empty;
  // The selector the selection was made with; a new one selects again (an
  // equal result keeps the selection's reference, so an inline selector does not loop).
  let selectedWith: ((value: T) => U) | undefined;
  const select = (value: T): boolean => {
    const { selector, compare } = callbacks;
    selectedWith = selector;
    const next = selector(value);
    if (selection !== empty && compare(selection, next)) return false;
    selection = next;
    return true;
  };
  const getSnapshot = () => {
    let value = snapshot;
    try {
      if (hasValue(source)) value = source.value;
    } catch {
      // A terminated current-value source may throw; retain its last snapshot.
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
  return {
    getSnapshot,
    subscribe(onStoreChange: () => void) {
      const subscription = (hasValue(source) ? source.pipe(skip(1)) : source)?.subscribe((value) => {
        snapshot = value;
        if (select(value)) onStoreChange();
      });
      return () => subscription?.unsubscribe();
    },
  };
}

/**
 * Callback identities may change; their selection/equality semantics stay stable per source. The
 * source's identity is its lifetime, so it must be stable across renders: a new observable each
 * render (an inline `.pipe(...)`) starts over each time and shows only the default. Memoize it.
 */
export function useObservableSelection<T, U>(
  source: Observable<T> | null | undefined,
  selector: (value: T) => U,
  options?: ObservableOptions<T, U>,
): U | undefined {
  // A source that changes on renders in a row is almost surely made inline;
  // say so once, instead of silently showing the default.
  const sourceRef = useRef(source);
  const changedInARow = useRef(0);
  if (sourceRef.current !== source) {
    sourceRef.current = source;
    changedInARow.current += 1;
    if (changedInARow.current === 2)
      console.warn(
        'useObservable: the source is a new observable on every render (an inline .pipe?). Memoize it: a new source starts a new lifetime and shows only its default.',
      );
  } else changedInARow.current = 0;
  const callbacksRef = useRef<ObservableCallbacks<T, U> | undefined>(undefined);
  const compare = options?.compare ?? sameShallow;
  if (!callbacksRef.current) {
    callbacksRef.current = { selector, compare };
  } else {
    callbacksRef.current.selector = selector;
    callbacksRef.current.compare = compare;
  }
  const callbacks = callbacksRef.current;
  // defaultValue only seeds a new source lifetime, so it is not a dependency.
  const store = useMemo(
    () => createStore(source, options && 'defaultValue' in options ? (options.defaultValue as T) : empty, callbacks),
    [source, callbacks],
  );
  return useSyncExternalStore(store.subscribe, store.getSnapshot, store.getSnapshot);
}
