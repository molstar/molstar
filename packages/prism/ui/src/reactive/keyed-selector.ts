/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Observable, type Subscription } from 'rxjs';
import type { Listenable } from './listen.js';

/**
 * Whether a key is the current one (the selected row of thousands): `is(key)` emits true or false, and a change of the
 * current value notifies only the two keys it concerns, not every row (`source.pipe(map((s) => s === key))` per row
 * would run every row's map on each change). `key(key)` is the same for `watch` (rendering.md §2): a source that tells
 * that key's row when it becomes or stops being the current one (the row reads the current value itself), so a move
 * renders two rows, not all of them. One upstream subscription, while anything listens.
 */
export function keyedSelector<K>(source: Observable<K>): {
  is(key: K): Observable<boolean>;
  key(key: K): Listenable;
} {
  const listeners = new Map<K, Set<() => void>>();
  let current: K | undefined;
  let has = false;
  let upstream: Subscription | undefined;
  let count = 0;
  const notify = (key: K | undefined) => {
    if (key === undefined && !has) return;
    const set = listeners.get(key as K);
    if (set) for (const fn of set) fn();
  };
  const add = (key: K, fn: () => void) => {
    let set = listeners.get(key);
    if (!set) listeners.set(key, (set = new Set()));
    set.add(fn);
    if (count++ === 0) {
      upstream = source.subscribe((v) => {
        if (has && v === current) return;
        const previous = current;
        const hadPrevious = has;
        current = v;
        has = true;
        if (hadPrevious) notify(previous);
        notify(v);
      });
    }
    return () => {
      set!.delete(fn);
      if (set!.size === 0) listeners.delete(key);
      if (--count === 0) {
        upstream?.unsubscribe();
        upstream = undefined;
        has = false;
      }
    };
  };
  return {
    is(key) {
      return new Observable<boolean>((subscriber) => {
        // A wrapper per subscriber: a set keeps each once.
        const stop = add(key, () => subscriber.next(has && current === key));
        subscriber.next(has && current === key);
        return stop;
      });
    },
    key(key) {
      return { listen: (fn) => add(key, () => fn()) };
    },
  };
}
