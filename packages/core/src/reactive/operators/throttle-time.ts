/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { lift, Subscriber } from '../observable.js';
import type { MonoTypeOperatorFunction } from '../types.js';

export interface ThrottleConfig {
  leading?: boolean;
  trailing?: boolean;
}

/** rxjs 7's `throttle` with a timer for the duration, step for step (its trailing value opens a new window). */
class ThrottleSubscriber<T> extends Subscriber<T> {
  private handle: ReturnType<typeof setTimeout> | undefined = undefined;
  private hasValue = false;
  private sendValue: T | undefined = undefined;
  private sourceDone = false;
  private readonly endThrottling = () => {
    this.handle = undefined;
    if (this.trailing) {
      this.send();
      if (this.sourceDone) super._complete();
    }
  };

  constructor(
    destination: Subscriber<T>,
    private readonly duration: number,
    private readonly leading: boolean,
    private readonly trailing: boolean,
  ) {
    super(destination);
  }

  private startThrottle() {
    this.handle = setTimeout(this.endThrottling, this.duration);
  }

  private send() {
    if (!this.hasValue) return;
    this.hasValue = false;
    const value = this.sendValue as T;
    this.sendValue = undefined;
    this.destination!.next(value);
    if (!this.sourceDone && !this.closed) this.startThrottle();
  }

  protected override _next(value: T): void {
    this.hasValue = true;
    this.sendValue = value;
    if (this.handle !== undefined) return;
    if (this.leading) this.send();
    else this.startThrottle();
  }

  protected override _complete(): void {
    this.sourceDone = true;
    if (this.trailing && this.hasValue && this.handle !== undefined) return;
    super._complete();
  }

  protected override _dispose(): void {
    if (this.handle !== undefined) clearTimeout(this.handle);
    this.handle = undefined;
    this.sendValue = undefined;
    super._dispose();
  }
}

/**
 * Tells a value, then ignores the source for `duration` ms (`leading`, default true); with `trailing` (default
 * false), the latest value ignored is told at the window's end and opens the next window. The scheduler argument is
 * rxjs's position, unused: timers are `setTimeout`.
 */
export function throttleTime<T>(
  duration: number,
  _scheduler?: undefined,
  config?: ThrottleConfig,
): MonoTypeOperatorFunction<T> {
  const leading = config?.leading ?? true;
  const trailing = config?.trailing ?? false;
  return (source) =>
    lift(source, (destination: Subscriber<T>) => new ThrottleSubscriber(destination, duration, leading, trailing));
}
