/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ObservableState, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { merge, Observable, type Subscription } from 'rxjs';
import { cancelFrameCallback, requestFrameCallback } from '@molstar/prism-ui/reactive/frame-state';

// A source that follows the plugin only while someone watches it (prism AGENTS.md: hidden UI holds no subscriptions;
// the pattern of the state tree's `watched()`, once): the events are listened to from the first subscriber to the last,
// bursts of them (a transaction fires several, a playing trajectory every frame) are read once per animation frame, and
// `value` is always the current read (fresh when nobody watches, so a caller that only asks pays for that one read).

export interface WatchedOptions<T> {
  /** Equal reads are not emitted. */
  readonly compare?: (a: T, b: T) => boolean;
  /** How events are coalesced: once per animation frame (default), or each one reads at once. */
  readonly coalesce?: 'frame' | 'none';
  /**
   * The state is computed from other watched states, which may hold a read for their frame: `value` is read again each
   * time, so it is current whatever they hold (the emissions follow them as they come).
   */
  readonly derived?: boolean;
}

type Frame = (cb: () => void) => () => void;

/** The next animation frame where there is one (a page), else a timer of about one. */
const nextFrame: Frame = (cb) => {
  if (typeof requestAnimationFrame === 'function') {
    const id = requestFrameCallback(cb);
    return () => cancelFrameCallback(id);
  }
  const id = setTimeout(cb, 16);
  return () => clearTimeout(id);
};

/**
 * `read()`, re-read when `events()` emit, while subscribed. `events` is called when the first subscriber comes and its
 * subscription released with the last; it may give several streams (any of them re-reads). `value` flushes a read that
 * waits for its frame, so it is never stale.
 */
export function watchedState<T>(
  read: () => T,
  events: () => Observable<unknown> | readonly Observable<unknown>[],
  options: WatchedOptions<T> = {},
): ReadonlyState<T> {
  const coalesce = options.coalesce ?? 'frame';
  // Made at the first read (a comparison never sees a value that is not one).
  let state: ObservableState<T> | undefined;
  let watchers = 0;
  let listening: Subscription | undefined;
  let cancelFrame: (() => void) | undefined;
  let dirty = false;

  const refresh = () => {
    cancelFrame?.();
    cancelFrame = undefined;
    dirty = false;
    const value = read();
    if (state) state.set(value);
    else state = observableState<T>(value, options.compare ? { compare: options.compare } : undefined);
  };
  const mark = () => {
    if (coalesce === 'none') return refresh();
    dirty = true;
    cancelFrame ??= nextFrame(() => {
      cancelFrame = undefined;
      if (dirty) refresh();
    });
  };

  const stream = new Observable<T>((subscriber) => {
    if (watchers++ === 0) {
      // Whatever happened while nobody looked: read again, then follow.
      refresh();
      const e = events();
      listening = (Array.isArray(e) ? merge(...e) : (e as Observable<unknown>)).subscribe(mark);
    }
    const sub = state!.asReadonly.subscribe(subscriber);
    return () => {
      sub.unsubscribe();
      if (--watchers > 0) return;
      listening?.unsubscribe();
      listening = undefined;
      cancelFrame?.();
      cancelFrame = undefined;
      dirty = false;
    };
  }) as ReadonlyState<T>;
  Object.defineProperty(stream, 'value', {
    enumerable: true,
    get: () => {
      if (watchers === 0 || options.derived) return read();
      if (dirty) refresh();
      return state!.value;
    },
  });
  return stream;
}
