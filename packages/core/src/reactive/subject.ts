/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Listened } from './listened.js';
import { Observable, type Subscriber, toSubscriber } from './observable.js';
import type { Subscription } from './subscription.js';
import type { Observer, PartialObserver, SubscriptionLike } from './types.js';

/**
 * A multicast stream with rxjs 7's semantics (reactive §3.2): `next` tells the listeners at once, nested; `error` and
 * `complete` end it (later `next` is ignored, a later subscriber gets the end at once); `unsubscribe` drops the
 * listeners and closes it.
 */
export class Subject<T> extends Listened<T> implements Observer<T>, SubscriptionLike {
  /** Whether `unsubscribe` was called (rxjs: not whether it ended). */
  closed = false;
  isStopped = false;
  hasError = false;
  thrownError: any = null;

  next(value: T): void {
    if (this.isStopped) return;
    this.emit(value);
  }

  error(err: any): void {
    if (this.isStopped) return;
    this.isStopped = true;
    this.hasError = true;
    this.thrownError = err;
    this.emitError(err);
  }

  complete(): void {
    if (this.isStopped) return;
    this.isStopped = true;
    this.emitComplete();
  }

  unsubscribe(): void {
    this.isStopped = true;
    this.closed = true;
    this.removeAll();
  }

  override subscribe(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): Subscription {
    const subscriber = toSubscriber(observerOrNext);
    if (this.hasError) subscriber.error(this.thrownError);
    else if (this.isStopped) subscriber.complete();
    else this._add(subscriber);
    return subscriber;
  }

  private _view: Observable<T> | undefined = undefined;

  /** A view without `next`, `error`, `complete`; made on first use and kept (rxjs makes one each call). */
  asObservable(): Observable<T> {
    return (this._view ??= new Observable<T>((subscriber: Subscriber<T>) => this.subscribe(subscriber)));
  }
}

/**
 * A subject with a current value: a subscriber gets it at once (nested, as rxjs), then each new one. `getValue` throws
 * the error once it has errored. It is a `ReadonlyState` (`value`, `listen`, a replaying stream), so `derive` follows it.
 */
export class BehaviorSubject<T> extends Subject<T> {
  constructor(private _value: T) {
    super();
  }

  get value(): T {
    return this.getValue();
  }

  getValue(): T {
    if (this.hasError) throw this.thrownError;
    return this._value;
  }

  override next(value: T): void {
    this._value = value;
    super.next(value);
  }

  override subscribe(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): Subscription {
    const subscriber = toSubscriber(observerOrNext);
    super.subscribe(subscriber);
    if (!subscriber.closed) subscriber.next(this._value);
    return subscriber;
  }
}
