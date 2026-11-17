/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from './observable.js';

// The shapes rxjs 7 names, so code and rxjs interoperate structurally (reactive §3.1).

export interface Unsubscribable {
  unsubscribe(): void;
}

export interface SubscriptionLike extends Unsubscribable {
  readonly closed: boolean;
}

export interface Observer<T> {
  next: (value: T) => void;
  error: (err: any) => void;
  complete: () => void;
}

export type PartialObserver<T> = Partial<Observer<T>>;

/** Anything with `subscribe(observer)`: an rxjs `Observable` is one. */
export interface Subscribable<T> {
  subscribe(observer: Partial<Observer<T>>): Unsubscribable;
}

export type TeardownLogic = Unsubscribable | (() => void) | void;

export type UnaryFunction<T, R> = (source: T) => R;
export type OperatorFunction<T, R> = UnaryFunction<Observable<T>, Observable<R>>;
export type MonoTypeOperatorFunction<T> = OperatorFunction<T, T>;

/**
 * What `from` and the operators that take inner sources accept: ours, any `Subscribable` (rxjs), a promise, an array
 * or an iterable. An object with `Symbol.observable` (or `'@@observable'`) is accepted at run time too.
 */
export type ObservableInput<T> = Observable<T> | Subscribable<T> | PromiseLike<T> | ArrayLike<T> | Iterable<T>;

export type ObservedValueOf<O> =
  O extends Observable<infer T>
    ? T
    : O extends PromiseLike<infer T>
      ? T
      : O extends ArrayLike<infer T>
        ? T
        : O extends Iterable<infer T>
          ? T
          : // An rxjs `Observable`: its `subscribe` is overloaded (inference reads the
            // last overload, which takes no observer); the last `forEach` still takes `next`.
            O extends { forEach(next: (value: infer T) => void, ...rest: any[]): any }
            ? T
            : O extends Subscribable<infer T>
              ? T
              : never;

export type ObservableInputTuple<T> = { [K in keyof T]: ObservableInput<T[K]> };

/** Something that tells of its changes with a plain callback. */
export interface Listenable<T = unknown> {
  listen(fn: (value: T) => void): () => void;
}

/**
 * A value that changes: its current `value`, `listen` (a plain callback, no `Subscription`), and a stream that replays
 * the current value then each new one. `State`, `derive`'s values and `BehaviorSubject` are each one.
 */
export type ReadonlyState<T> = Observable<T> & Listenable<T> & { readonly value: T };
