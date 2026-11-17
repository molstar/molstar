/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { from, Observable, Subscriber } from './observable.js';
import type { ObservableInput, ObservableInputTuple } from './types.js';

/** Each value, then complete; cold. */
export function of<A extends readonly unknown[]>(...values: A): Observable<A[number]> {
  return new Observable((subscriber) => {
    for (let i = 0; i < values.length; i++) {
      if (subscriber.closed) return;
      subscriber.next(values[i]);
    }
    subscriber.complete();
  });
}

/** Completes at once. */
export const EMPTY: Observable<never> = new Observable<never>((subscriber) => subscriber.complete());

/** Never tells anything. */
export const NEVER: Observable<never> = new Observable<never>(() => {});

/** `0` after `dueTime` ms, then complete; with `period`, `1, 2, …` every `period` ms after it. */
export function timer(dueTime: number, period?: number): Observable<number> {
  return new Observable<number>((subscriber) => {
    let n = 0;
    let handle: ReturnType<typeof setTimeout> | undefined = setTimeout(function tick() {
      subscriber.next(n++);
      if (subscriber.closed) return;
      if (period === undefined || period < 0) {
        handle = undefined;
        subscriber.complete();
      } else {
        handle = setTimeout(tick, period);
      }
    }, dueTime);
    return () => {
      if (handle !== undefined) clearTimeout(handle);
      handle = undefined;
    };
  });
}

/** The state one `merge` subscription shares between its inner subscribers. */
interface MergeState {
  active: number;
}

class MergeInner<T> extends Subscriber<T> {
  constructor(
    destination: Subscriber<T>,
    private readonly state: MergeState,
  ) {
    super(destination);
  }

  protected override _complete(): void {
    const destination = this.destination!;
    this.unsubscribe();
    if (--this.state.active === 0) destination.complete();
  }
}

/** Every value of every source as it comes; completes when all have, errors with the first error. */
export function merge<A extends readonly unknown[]>(...sources: [...ObservableInputTuple<A>]): Observable<A[number]> {
  return new Observable<A[number]>((destination) => {
    const n = sources.length;
    if (n === 0) {
      destination.complete();
      return;
    }
    const state: MergeState = { active: n };
    for (let i = 0; i < n && !destination.closed; i++) {
      const inner = new MergeInner<A[number]>(destination, state);
      destination.add(inner);
      from(sources[i] as ObservableInput<A[number]>).subscribe(inner);
    }
  });
}

/** The state one `combineLatest` subscription shares between its inner subscribers. */
class CombineState {
  readonly values: unknown[];
  /** Sources that have not told a value yet. */
  waiting: number;
  active: number;
  readonly seen: boolean[];
  constructor(n: number) {
    this.values = new Array(n);
    this.seen = new Array(n).fill(false);
    this.waiting = n;
    this.active = n;
  }
}

class CombineInner<T> extends Subscriber<T> {
  constructor(
    destination: Subscriber<unknown[]>,
    private readonly state: CombineState,
    private readonly index: number,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    const state = this.state;
    state.values[this.index] = value;
    if (!state.seen[this.index]) {
      state.seen[this.index] = true;
      state.waiting--;
    }
    // A new array each time, as rxjs: a subscriber may keep the one it got.
    if (state.waiting === 0) this.destination!.next(state.values.slice());
  }

  protected override _complete(): void {
    const destination = this.destination!;
    this.unsubscribe();
    if (--this.state.active === 0) destination.complete();
  }
}

/**
 * An array of each source's latest value, from when every source has told one; completes when all have (at once when
 * `sources` is empty), errors with the first error.
 */
export function combineLatest<A extends readonly unknown[]>(
  sources: readonly [...ObservableInputTuple<A>],
): Observable<A> {
  return new Observable<A>((destination) => {
    const n = sources.length;
    if (n === 0) {
      destination.complete();
      return;
    }
    const state = new CombineState(n);
    for (let i = 0; i < n && !destination.closed; i++) {
      const inner = new CombineInner(destination as unknown as Subscriber<unknown[]>, state, i);
      destination.add(inner);
      from(sources[i] as ObservableInput<unknown>).subscribe(inner);
    }
  });
}
