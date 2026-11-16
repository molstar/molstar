/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { listen, type ReadonlyState, type Subscribable } from '@molstar/core/reactive';
import { sameShallow } from '@molstar/core/util/object';
import { type MutableRefObject, useEffect, useRef, useSyncExternalStore } from 'react';

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

const empty = Symbol('empty');
type Raw<T> = T | typeof empty;

/** A source with a current value: read it, never wait for (or count) its replay. */
function hasValue<T>(source: ReactiveSource<T>): source is ReadonlyState<T> {
  return source != null && 'value' in source;
}

/** A stream that ended in an error (ours or rxjs's subjects say so): its `value` throws that error from then on. */
function endedInError(source: unknown): boolean {
  return (source as { hasError?: unknown }).hasError === true;
}

/**
 * One source's lifetime (react §1): its raw value, nothing selected. A source with a current value is read through
 * `value` and followed through core's `listen` (changes only: subscribing never replays into the store, so the host is
 * told once per change and never on subscribing). A stream without one keeps what it emits and tells the host of each
 * change, including one emitted while subscribing (a cold `of(x)`): React reads the snapshot again after subscribing,
 * but Preact's `useSyncExternalStore` (and React's shim) checks before and never after, so such a value would be lost.
 */
class Store<T> {
  #raw: Raw<T>;

  constructor(
    readonly source: ReactiveSource<T>,
    initial: Raw<T>,
  ) {
    this.#raw = initial;
  }

  /** The raw value now; `empty` for a stream before its first value with no default. */
  raw = (): Raw<T> => {
    const source = this.source;
    if (hasValue(source)) {
      try {
        this.#raw = source.value;
      } catch (e) {
        // Only a stream that ended in an error keeps its last value; a throw from a derived value (a bug) is the
        // host's to see (an error boundary, the console).
        if (!endedInError(source)) throw e;
      }
    }
    return this.#raw;
  };

  subscribe = (onStoreChange: () => void): (() => void) => {
    const source = this.source;
    if (!source) return () => {};
    // The next snapshot reads `value`: the listener only says it changed.
    if (hasValue(source)) return listen(source, onStoreChange);
    const subscription = source.subscribe({
      next: (value) => {
        if (Object.is(value, this.#raw)) return;
        this.#raw = value;
        onStoreChange();
      },
    });
    return () => subscription.unsubscribe();
  };
}

type Selected<T, U> = {
  store: Store<T> | undefined;
  selector: ((value: T) => U) | undefined;
  raw: Raw<T>;
  value: U | undefined;
};

/** The committed store when it is the source's, else the pending one, made for the source when it is not. */
function storeFor<T>(
  source: ReactiveSource<T>,
  committed: MutableRefObject<Store<T> | undefined>,
  pending: MutableRefObject<Store<T> | undefined>,
  options: { readonly defaultValue?: T } | undefined,
): Store<T> {
  const done = committed.current;
  if (done && done.source === source) return done;
  const next = pending.current;
  if (next && next.source === source) return next;
  // `defaultValue` only seeds a new source's lifetime.
  const store = new Store<T>(source, options && 'defaultValue' in options ? (options.defaultValue as T) : empty);
  pending.current = store;
  return store;
}

/**
 * A selection of a source's value, rendering again only when the selection changes (`compare`, shallow by default).
 * The selector and the comparer may be new on each render; the selection made with the last committed render is the
 * one a new selection is compared with, so an equal result keeps its reference. A new selector selects again on the
 * render it comes with: memoize an expensive one (`useCallback`). The source's identity is its lifetime, so it must be
 * stable across renders: a new one (an inline `.pipe(...)` or `derive(...)`) starts over every render and shows only
 * the default. Memoize it; a source new on committed renders in a row warns once.
 */
export function useReactiveSelection<T, U>(
  source: ReactiveSource<T>,
  selector: (value: T) => U,
  options?: ReactiveOptions<T, U>,
): U | undefined {
  const compare = options?.compare ?? sameShallow;
  // The store of the committed render's source, and one made by a render not committed yet (react §1): a render React
  // discards never replaces the committed source's store (its stream's latest value). Kept by source, not in `useMemo`
  // (which React may drop).
  const committed = useRef<Store<T> | undefined>(undefined);
  const pending = useRef<Store<T> | undefined>(undefined);
  const store = storeFor(source, committed, pending, options);
  // The last selection, keyed by everything it was made from (the store, the selector, the raw value): a read from a
  // render React discarded only misses, never answers for another. A new value is compared with it (`compare`, this
  // render's), so an equal selection keeps its reference; a comparer that changed alone selects nothing again.
  const cell = useRef<Selected<T, U>>({ store: undefined, selector: undefined, raw: empty, value: undefined });
  const getSnapshot = (): U | undefined => {
    const raw = store.raw();
    if (raw === empty) return undefined;
    const c = cell.current;
    const same = c.store === store;
    if (same && c.selector === selector && Object.is(c.raw, raw)) return c.value;
    const next = selector(raw);
    const value = same && c.raw !== empty && compare(c.value as U, next) ? (c.value as U) : next;
    cell.current = { store, selector, raw, value };
    return value;
  };
  const value = useSyncExternalStore(store.subscribe, getSnapshot, getSnapshot);

  const seen = useRef({ source, inARow: 0, warned: false });
  useEffect(() => {
    committed.current = store;
    if (pending.current === store) pending.current = undefined;
    // A source new on committed renders in a row is almost surely made inline: said once.
    const s = seen.current;
    if (s.source === source) {
      s.inARow = 0;
      return;
    }
    s.source = source;
    s.inARow++;
    if (s.inARow === 2 && !s.warned) {
      s.warned = true;
      console.warn(
        'useReactive: the source is new on every render (an inline pipe or derive?). Memoize it: a new source starts a new lifetime and shows only its default.',
      );
    }
  });
  return value;
}
