/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Observable, type Subscriber, type Subscription } from 'rxjs';

/**
 * Whether a key is the current one (the selected row of thousands): `is(key)` emits true or false, and a change of the
 * current value notifies only the two keys it concerns, not every row (`source.pipe(map((s) => s === key))` per row
 * would run every row's map on each change). One upstream subscription, while anything listens.
 */
export function keyedSelector<K>(source: Observable<K>): { is(key: K): Observable<boolean> } {
  const listeners = new Map<K, Set<Subscriber<boolean>>>();
  let current: K | undefined;
  let has = false;
  let upstream: Subscription | undefined;
  let count = 0;
  const notify = (key: K | undefined, on: boolean) => {
    if (key === undefined && !has) return;
    const set = listeners.get(key as K);
    if (set) for (const s of set) s.next(on);
  };
  return {
    is(key) {
      return new Observable<boolean>((subscriber) => {
        let set = listeners.get(key);
        if (!set) listeners.set(key, (set = new Set()));
        set.add(subscriber);
        if (count++ === 0) {
          upstream = source.subscribe((v) => {
            if (has && v === current) return;
            const previous = current;
            const hadPrevious = has;
            current = v;
            has = true;
            if (hadPrevious) notify(previous, false);
            notify(v, true);
          });
        }
        subscriber.next(has && current === key);
        return () => {
          set!.delete(subscriber);
          if (set!.size === 0) listeners.delete(key);
          if (--count === 0) {
            upstream?.unsubscribe();
            upstream = undefined;
            has = false;
          }
        };
      });
    },
  };
}
