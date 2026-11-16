/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lift, Subscriber } from '../observable.js';
import type { MonoTypeOperatorFunction, PartialObserver } from '../types.js';

class TapSubscriber<T> extends Subscriber<T> {
  constructor(
    destination: Subscriber<T>,
    private readonly observer: PartialObserver<T>,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    try {
      this.observer.next?.(value);
    } catch (e) {
      this.error(e);
      return;
    }
    this.destination!.next(value);
  }

  protected override _error(err: any): void {
    try {
      this.observer.error?.(err);
    } catch (e) {
      err = e;
    }
    super._error(err);
  }

  protected override _complete(): void {
    try {
      this.observer.complete?.();
    } catch (e) {
      super._error(e);
      return;
    }
    super._complete();
  }
}

/** Calls `observerOrNext` with each notification, then passes it on; what it throws is the stream's error. */
export function tap<T>(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): MonoTypeOperatorFunction<T> {
  if (!observerOrNext) return (source) => source;
  const observer: PartialObserver<T> = typeof observerOrNext === 'function' ? { next: observerOrNext } : observerOrNext;
  return (source) => lift(source, (destination: Subscriber<T>) => new TapSubscriber(destination, observer));
}
