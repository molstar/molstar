/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Listened } from '../listened.js';
import { type Observable, Subscriber } from '../observable.js';
import type { Subscription } from '../subscription.js';
import type { MonoTypeOperatorFunction, PartialObserver } from '../types.js';

class ShareSource<T> extends Subscriber<T> {
  constructor(private readonly shared: Shared<T>) {
    super();
  }

  protected override _next(value: T): void {
    this.shared.forward(value);
  }

  protected override _error(err: any): void {
    this.unsubscribe();
    this.shared.end(err, true);
  }

  protected override _complete(): void {
    this.unsubscribe();
    this.shared.end(undefined, false);
  }
}

/**
 * One upstream subscription for every subscriber: made with the first, ended with the last (ref-counted, as rxjs's
 * `share()`); an upstream error or completion goes to every subscriber and resets it, so the next subscriber
 * subscribes upstream anew. No replay.
 */
class Shared<T> extends Listened<T> {
  private upstream: ShareSource<T> | undefined = undefined;

  constructor(private readonly source: Observable<T>) {
    super();
  }

  override subscribe(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): Subscription {
    const subscriber = super.subscribe(observerOrNext);
    // Connect after the first subscriber is in the list, so a synchronous source reaches it.
    if (this.upstream === undefined && this.observed) {
      const upstream = new ShareSource(this);
      this.upstream = upstream;
      this.source.subscribe(upstream);
    }
    return subscriber;
  }

  /** @internal */
  forward(value: T): void {
    this.emit(value);
  }

  /** @internal */
  end(err: any, isError: boolean): void {
    this.upstream = undefined;
    if (isError) this.emitError(err);
    else this.emitComplete();
  }

  protected override onLast(): void {
    const upstream = this.upstream;
    this.upstream = undefined;
    upstream?.unsubscribe();
  }
}

/** Shares one subscription to the source among all subscribers (ref-counted, reset on end, no replay). */
export function share<T>(): MonoTypeOperatorFunction<T> {
  return (source) => new Shared(source);
}
