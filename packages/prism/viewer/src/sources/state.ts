/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import type { BehaviorSubject, Observable, Subscription } from 'rxjs';

// One shape for every source (state.md §4): the plugin's patterns (a BehaviorSubject, a mutable field with a void
// `changed`, events to accumulate) turned into `ReadonlyState<T>`: a synchronous `value`, emitting only on change. Each
// helper takes a `SourceScope`, which holds the subscriptions until the viewer goes.

export class SourceScope {
  readonly #subs: Subscription[] = [];
  readonly #disposers: (() => void)[] = [];
  #disposed = false;

  add(sub: Subscription | (() => void)): void {
    if (this.#disposed) {
      if (typeof sub === 'function') sub();
      else sub.unsubscribe();
      return;
    }
    if (typeof sub === 'function') this.#disposers.push(sub);
    else this.#subs.push(sub);
  }

  dispose(): void {
    if (this.#disposed) return;
    this.#disposed = true;
    for (const s of this.#subs) s.unsubscribe();
    for (const d of this.#disposers) d();
    this.#subs.length = 0;
    this.#disposers.length = 0;
  }
}

/** A plugin `BehaviorSubject`, taken as is with equality (state.md §4, row 1). */
export function fromBehavior<T>(
  scope: SourceScope,
  subject: BehaviorSubject<T>,
  compare?: (a: T, b: T) => boolean,
): ReadonlyState<T> {
  return fromChanged(scope, () => subject.value, subject, compare);
}

/**
 * A value re-read on an event, emitted when it differs (state.md §4, row 2): a mutable field with a void `changed`,
 * or a getter whose writes the viewer makes emit.
 */
export function fromChanged<T>(
  scope: SourceScope,
  read: () => T,
  event: Observable<unknown>,
  compare?: (a: T, b: T) => boolean,
): ReadonlyState<T> {
  const state = observableState<T>(read(), compare ? { compare } : undefined);
  scope.add(event.subscribe(() => state.set(read())));
  return state.asReadonly;
}

/** Events accumulated into a list (state.md §4, row 5). */
export function scanList<E, T>(
  scope: SourceScope,
  events: Observable<E>,
  reduce: (list: readonly T[], event: E) => readonly T[],
  initial: readonly T[] = [],
): ReadonlyState<readonly T[]> {
  const state = observableState<readonly T[]>(initial);
  scope.add(events.subscribe((e) => state.set(reduce(state.value, e))));
  return state.asReadonly;
}

/** Equal arrays, item by item. */
export function sameItems<T>(a: readonly T[], b: readonly T[], same: (x: T, y: T) => boolean = Object.is): boolean {
  if (a === b) return true;
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (!same(a[i], b[i])) return false;
  return true;
}
