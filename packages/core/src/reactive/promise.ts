/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { from, Subscriber } from './observable.js';
import type { ObservableInput } from './types.js';

/** A stream that completed with no value (rxjs's `EmptyError`). */
export class EmptyError extends Error {
  constructor() {
    super('no elements in sequence');
    this.name = 'EmptyError';
  }
}

export interface ValueFromConfig<D> {
  defaultValue: D;
}

class FirstValueSubscriber<T> extends Subscriber<T> {
  constructor(
    private readonly resolve: (value: T) => void,
    private readonly reject: (err: any) => void,
    private readonly fallback: { value: T } | undefined,
  ) {
    super();
  }

  protected override _next(value: T): void {
    this.resolve(value);
    this.unsubscribe();
  }

  protected override _error(err: any): void {
    this.reject(err);
    this.unsubscribe();
  }

  protected override _complete(): void {
    if (this.fallback) this.resolve(this.fallback.value);
    else this.reject(new EmptyError());
    this.unsubscribe();
  }
}

/** The first value, then unsubscribes; rejects with the error, or `EmptyError` on a completion with none. */
export function firstValueFrom<T, D>(source: ObservableInput<T>, config: ValueFromConfig<D>): Promise<T | D>;
export function firstValueFrom<T>(source: ObservableInput<T>): Promise<T>;
export function firstValueFrom<T, D>(source: ObservableInput<T>, config?: ValueFromConfig<D>): Promise<T | D> {
  return new Promise<T | D>((resolve, reject) => {
    const fallback = config ? { value: config.defaultValue as T | D } : undefined;
    from(source).subscribe(new FirstValueSubscriber<T | D>(resolve, reject, fallback));
  });
}

class LastValueSubscriber<T> extends Subscriber<T> {
  private hasValue = false;
  private last: T | undefined = undefined;

  constructor(
    private readonly resolve: (value: T) => void,
    private readonly reject: (err: any) => void,
    private readonly fallback: { value: T } | undefined,
  ) {
    super();
  }

  protected override _next(value: T): void {
    this.hasValue = true;
    this.last = value;
  }

  protected override _error(err: any): void {
    this.reject(err);
    this.unsubscribe();
  }

  protected override _complete(): void {
    if (this.hasValue) this.resolve(this.last as T);
    else if (this.fallback) this.resolve(this.fallback.value);
    else this.reject(new EmptyError());
    this.unsubscribe();
  }
}

/** The last value once the source completes; rejects with the error, or `EmptyError` on a completion with none. */
export function lastValueFrom<T, D>(source: ObservableInput<T>, config: ValueFromConfig<D>): Promise<T | D>;
export function lastValueFrom<T>(source: ObservableInput<T>): Promise<T>;
export function lastValueFrom<T, D>(source: ObservableInput<T>, config?: ValueFromConfig<D>): Promise<T | D> {
  return new Promise<T | D>((resolve, reject) => {
    const fallback = config ? { value: config.defaultValue as T | D } : undefined;
    from(source).subscribe(new LastValueSubscriber<T | D>(resolve, reject, fallback));
  });
}
