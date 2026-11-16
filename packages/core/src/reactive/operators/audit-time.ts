/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lift, Subscriber } from '../observable.js';
import type { MonoTypeOperatorFunction } from '../types.js';

class AuditSubscriber<T> extends Subscriber<T> {
  private handle: ReturnType<typeof setTimeout> | undefined = undefined;
  private hasValue = false;
  private last: T | undefined = undefined;
  private sourceDone = false;
  private readonly fire = () => {
    this.handle = undefined;
    if (this.hasValue) {
      const value = this.last as T;
      this.hasValue = false;
      this.last = undefined;
      this.destination!.next(value);
    }
    if (this.sourceDone) super._complete();
  };

  constructor(
    destination: Subscriber<T>,
    private readonly duration: number,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    this.last = value;
    this.hasValue = true;
    if (this.handle === undefined) this.handle = setTimeout(this.fire, this.duration);
  }

  protected override _complete(): void {
    // A value waiting for its window is told at the window's end, then the stream completes (rxjs's `audit`).
    if (this.hasValue && this.handle !== undefined) this.sourceDone = true;
    else super._complete();
  }

  protected override _dispose(): void {
    if (this.handle !== undefined) clearTimeout(this.handle);
    this.handle = undefined;
    this.last = undefined;
    super._dispose();
  }
}

/** A window of `duration` ms opens at a value; at its end, the latest value in it is told. */
export function auditTime<T>(duration: number): MonoTypeOperatorFunction<T> {
  return (source) => lift(source, (destination: Subscriber<T>) => new AuditSubscriber(destination, duration));
}
