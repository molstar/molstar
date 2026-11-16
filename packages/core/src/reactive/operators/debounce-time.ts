/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lift, Subscriber } from '../observable.js';
import type { MonoTypeOperatorFunction } from '../types.js';

class DebounceSubscriber<T> extends Subscriber<T> {
  private handle: ReturnType<typeof setTimeout> | undefined = undefined;
  private hasValue = false;
  private last: T | undefined = undefined;
  private readonly fire = () => {
    this.handle = undefined;
    this.flush();
  };

  constructor(
    destination: Subscriber<T>,
    private readonly dueTime: number,
  ) {
    super(destination);
  }

  private flush() {
    if (!this.hasValue) return;
    const value = this.last as T;
    this.hasValue = false;
    this.last = undefined;
    this.destination!.next(value);
  }

  protected override _next(value: T): void {
    this.last = value;
    this.hasValue = true;
    if (this.handle !== undefined) clearTimeout(this.handle);
    this.handle = setTimeout(this.fire, this.dueTime);
  }

  protected override _complete(): void {
    // The pending value is told before the end (rxjs).
    if (this.handle !== undefined) clearTimeout(this.handle);
    this.handle = undefined;
    this.flush();
    super._complete();
  }

  protected override _dispose(): void {
    if (this.handle !== undefined) clearTimeout(this.handle);
    this.handle = undefined;
    this.last = undefined;
    super._dispose();
  }
}

/** The latest value once `dueTime` ms pass with no newer one; a pending value is told on complete. */
export function debounceTime<T>(dueTime: number): MonoTypeOperatorFunction<T> {
  return (source) => lift(source, (destination: Subscriber<T>) => new DebounceSubscriber(destination, dueTime));
}
