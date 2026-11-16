/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { reportUnhandledError } from './config.js';
import { Subscription } from './subscription.js';
import type {
  ObservableInput,
  ObservedValueOf,
  Observer,
  OperatorFunction,
  PartialObserver,
  Subscribable,
  TeardownLogic,
} from './types.js';

/** What holds a subscriber in its listener list (`Listened`): told when the subscriber ends. */
export interface SubscriberOwner {
  /** @internal */
  _remove(subscriber: Subscriber<any>): void;
}

/**
 * The observer side of one subscription, and the subscription itself. Stops at the first `error` or `complete`;
 * `unsubscribe` stops it too. Operators subclass it (`_next`, `_error`, `_complete`) with their state in fields: one
 * object per stage per subscription (reactive §5). When held in a listener list (`_owner`), it is its own slot.
 */
export class Subscriber<T> extends Subscription implements Observer<T> {
  /** @internal */
  isStopped = false;
  protected destination: Subscriber<any> | undefined;
  /** @internal The listener list holding this one, and its place there. */
  _owner: SubscriberOwner | undefined = undefined;
  /** @internal */
  _index = -1;

  constructor(destination?: Subscriber<any>) {
    super();
    this.destination = destination;
  }

  next(value: T): void {
    if (!this.isStopped) this._next(value);
  }

  error(err: any): void {
    if (this.isStopped) return;
    this.isStopped = true;
    this._error(err);
  }

  complete(): void {
    if (this.isStopped) return;
    this.isStopped = true;
    this._complete();
  }

  protected _next(value: T): void {
    this.destination!.next(value);
  }

  protected _error(err: any): void {
    try {
      this.destination!.error(err);
    } finally {
      this.unsubscribe();
    }
  }

  protected _complete(): void {
    try {
      this.destination!.complete();
    } finally {
      this.unsubscribe();
    }
  }

  protected override _dispose(): void {
    this.isStopped = true;
    const owner = this._owner;
    if (owner !== undefined) owner._remove(this);
    this.destination = undefined;
  }
}

/**
 * A subscriber at the end of a chain: the user's function or observer. What it throws is reported
 * (`reportUnhandledError`), not thrown to the producer, so the producer and the other listeners go on (rxjs).
 */
export class ConsumerSubscriber<T> extends Subscriber<T> {
  private _fn: ((value: T) => void) | undefined;
  private _observer: PartialObserver<T> | undefined;

  constructor(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null) {
    super();
    if (typeof observerOrNext === 'function') {
      this._fn = observerOrNext;
      this._observer = undefined;
    } else {
      this._fn = undefined;
      this._observer = observerOrNext ?? undefined;
    }
  }

  protected override _next(value: T): void {
    try {
      const fn = this._fn;
      if (fn !== undefined) fn(value);
      else this._observer?.next?.(value);
    } catch (e) {
      reportUnhandledError(e);
    }
  }

  protected override _error(err: any): void {
    const observer = this._observer;
    try {
      if (observer?.error) observer.error(err);
      else reportUnhandledError(err);
    } catch (e) {
      reportUnhandledError(e);
    }
    this.unsubscribe();
  }

  protected override _complete(): void {
    try {
      this._observer?.complete?.();
    } catch (e) {
      reportUnhandledError(e);
    }
    this.unsubscribe();
  }

  protected override _dispose(): void {
    super._dispose();
    this._fn = undefined;
    this._observer = undefined;
  }
}

/** A subscriber for `observerOrNext`: itself when it is one (an operator's), else a `ConsumerSubscriber`. */
export function toSubscriber<T>(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): Subscriber<T> {
  return observerOrNext instanceof Subscriber ? observerOrNext : new ConsumerSubscriber(observerOrNext);
}

/**
 * The interop key rxjs reads: `Symbol.observable` where defined, else `'@@observable'` (rxjs no longer polyfills it).
 * `Symbol.observable` is typed (below) only for the interop interface; it may be undefined at run time.
 */
const symbolObservable: string | symbol =
  (typeof Symbol === 'function' && (Symbol as { observable?: symbol }).observable) || '@@observable';

/**
 * A stream, cold: each `subscribe` runs the producer (`new Observable(subscriber => teardown)`). Structurally what
 * rxjs 7's `Observable` is to its users (`subscribe`, `pipe`), and an interop observable for rxjs's `from`.
 */
export class Observable<T> {
  constructor(subscribe?: (this: Observable<T>, subscriber: Subscriber<T>) => TeardownLogic) {
    if (subscribe) this._subscribe = subscribe;
  }

  /** The producer: a subclass overrides it, or the constructor's function is it. */
  protected _subscribe(_subscriber: Subscriber<T>): TeardownLogic {
    return;
  }

  subscribe(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): Subscription {
    const subscriber = toSubscriber(observerOrNext);
    try {
      subscriber.add(this._subscribe(subscriber));
    } catch (e) {
      if (subscriber.isStopped) reportUnhandledError(e);
      else subscriber.error(e);
    }
    return subscriber;
  }

  pipe(): Observable<T>;
  pipe<A>(op1: OperatorFunction<T, A>): Observable<A>;
  pipe<A, B>(op1: OperatorFunction<T, A>, op2: OperatorFunction<A, B>): Observable<B>;
  pipe<A, B, C>(op1: OperatorFunction<T, A>, op2: OperatorFunction<A, B>, op3: OperatorFunction<B, C>): Observable<C>;
  pipe<A, B, C, D>(
    op1: OperatorFunction<T, A>,
    op2: OperatorFunction<A, B>,
    op3: OperatorFunction<B, C>,
    op4: OperatorFunction<C, D>,
  ): Observable<D>;
  pipe<A, B, C, D, E>(
    op1: OperatorFunction<T, A>,
    op2: OperatorFunction<A, B>,
    op3: OperatorFunction<B, C>,
    op4: OperatorFunction<C, D>,
    op5: OperatorFunction<D, E>,
  ): Observable<E>;
  pipe<A, B, C, D, E, F>(
    op1: OperatorFunction<T, A>,
    op2: OperatorFunction<A, B>,
    op3: OperatorFunction<B, C>,
    op4: OperatorFunction<C, D>,
    op5: OperatorFunction<D, E>,
    op6: OperatorFunction<E, F>,
  ): Observable<F>;
  pipe<A, B, C, D, E, F, G>(
    op1: OperatorFunction<T, A>,
    op2: OperatorFunction<A, B>,
    op3: OperatorFunction<B, C>,
    op4: OperatorFunction<C, D>,
    op5: OperatorFunction<D, E>,
    op6: OperatorFunction<E, F>,
    op7: OperatorFunction<F, G>,
  ): Observable<G>;
  pipe<A, B, C, D, E, F, G, H>(
    op1: OperatorFunction<T, A>,
    op2: OperatorFunction<A, B>,
    op3: OperatorFunction<B, C>,
    op4: OperatorFunction<C, D>,
    op5: OperatorFunction<D, E>,
    op6: OperatorFunction<E, F>,
    op7: OperatorFunction<F, G>,
    op8: OperatorFunction<G, H>,
  ): Observable<H>;
  pipe<A, B, C, D, E, F, G, H, I>(
    op1: OperatorFunction<T, A>,
    op2: OperatorFunction<A, B>,
    op3: OperatorFunction<B, C>,
    op4: OperatorFunction<C, D>,
    op5: OperatorFunction<D, E>,
    op6: OperatorFunction<E, F>,
    op7: OperatorFunction<F, G>,
    op8: OperatorFunction<G, H>,
    op9: OperatorFunction<H, I>,
  ): Observable<I>;
  pipe(...operations: OperatorFunction<any, any>[]): Observable<unknown>;
  // A rest array per `pipe` call: pipes are built at setup, not per value.
  pipe(...operations: OperatorFunction<any, any>[]): Observable<any> {
    let result: Observable<any> = this;
    for (let i = 0; i < operations.length; i++) result = operations[i](result);
    return result;
  }
}

/** An operator applied at subscribe time: `make(destination)` is the stage's subscriber, subscribed to `source`. */
class LiftedObservable<T, R> extends Observable<R> {
  constructor(
    private readonly source: Observable<T>,
    private readonly make: (destination: Subscriber<R>) => Subscriber<T>,
  ) {
    super();
  }

  protected override _subscribe(destination: Subscriber<R>): TeardownLogic {
    return this.source.subscribe(this.make(destination));
  }
}

/** The stage `source.pipe(op)` makes: one `Observable` at pipe time, one subscriber (`make`) per subscription. */
export function lift<T, R>(source: Observable<T>, make: (destination: Subscriber<R>) => Subscriber<T>): Observable<R> {
  return new LiftedObservable(source, make);
}

/**
 * Whether `value` is a stream: ours, or anything with `subscribe` and `pipe` (an rxjs `Observable`), so code that takes
 * "any stream" takes both while the workspace moves off rxjs.
 */
export function isObservable(value: unknown): value is Observable<unknown> {
  if (value instanceof Observable) return true;
  return (
    !!value &&
    typeof (value as Observable<unknown>).subscribe === 'function' &&
    typeof (value as Observable<unknown>).pipe === 'function'
  );
}

/** Ours as is; an interop observable, a `Subscribable` (rxjs), a promise, an array or an iterable as ours. */
export function from<O extends ObservableInput<any>>(input: O): Observable<ObservedValueOf<O>> {
  if (input instanceof Observable) return input as Observable<ObservedValueOf<O>>;
  const any = input as any;
  if (any == null) throw new TypeError(`from: ${input} is not an observable input`);
  const interop = any[symbolObservable] ?? any['@@observable'];
  if (typeof interop === 'function') {
    return new Observable((subscriber) => (interop.call(any) as Subscribable<any>).subscribe(observerOf(subscriber)));
  }
  if (typeof any.subscribe === 'function') {
    return new Observable((subscriber) => (any as Subscribable<any>).subscribe(observerOf(subscriber)));
  }
  if (typeof any.then === 'function') {
    return new Observable((subscriber) => {
      (any as PromiseLike<any>).then(
        (value) => {
          if (subscriber.closed) return;
          subscriber.next(value);
          subscriber.complete();
        },
        (err) => subscriber.error(err),
      );
    });
  }
  if (typeof any.length === 'number') {
    return new Observable((subscriber) => {
      for (let i = 0; i < any.length && !subscriber.closed; i++) subscriber.next(any[i]);
      subscriber.complete();
    });
  }
  if (typeof any[Symbol.iterator] === 'function') {
    return new Observable((subscriber) => {
      for (const value of any as Iterable<any>) {
        if (subscriber.closed) return;
        subscriber.next(value);
      }
      subscriber.complete();
    });
  }
  throw new TypeError(`from: ${input} is not an observable input`);
}

/**
 * Our subscriber as a plain observer for a foreign `subscribe`: rxjs would treat a `Subscriber` of ours as a partial
 * observer anyway; the bound methods keep `this` whatever the foreign side does.
 */
function observerOf<T>(subscriber: Subscriber<T>): Observer<T> {
  return {
    next: (v) => subscriber.next(v),
    error: (e) => subscriber.error(e),
    complete: () => subscriber.complete(),
  };
}

declare global {
  interface SymbolConstructor {
    // As rxjs declares it (the same declaration merges): the key of its `InteropObservable`.
    readonly observable: symbol;
  }
}

// The interop by type: rxjs's `InteropObservable` is `{ [Symbol.observable]: () => Subscribable<T> }`, so rxjs's `from`,
// `combineLatest`, `switchMap` take ours. At run time the method is set below, under the key rxjs reads.
export interface Observable<T> {
  [Symbol.observable](): Observable<T>;
}

function interopSelf(this: Observable<unknown>) {
  return this;
}
// Both keys: whichever rxjs read when it loaded (a `Symbol.observable` a polyfill defines later).
(Observable.prototype as any)[symbolObservable] = interopSelf;
(Observable.prototype as any)['@@observable'] = interopSelf;
