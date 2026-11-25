/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { from, lift, Subscriber } from '../observable.js';
import type { ObservableInput, ObservedValueOf, OperatorFunction } from '../types.js';

class SwitchInner<R> extends Subscriber<R> {
  constructor(
    destination: Subscriber<R>,
    private readonly outer: SwitchMapSubscriber<any, any>,
  ) {
    super(destination);
  }

  protected override _complete(): void {
    this.unsubscribe();
    this.outer.innerDone(this);
  }
}

class SwitchMapSubscriber<T, R> extends Subscriber<T> {
  private index = 0;
  private inner: SwitchInner<R> | undefined = undefined;
  private outerDone = false;

  constructor(
    destination: Subscriber<R>,
    private readonly project: (value: T, index: number) => ObservableInput<R>,
  ) {
    super(destination);
  }

  protected override _next(value: T): void {
    // The previous inner ends before the next starts.
    this.inner?.unsubscribe();
    this.inner = undefined;
    let input: ObservableInput<R>;
    try {
      input = this.project(value, this.index++);
    } catch (e) {
      this.error(e);
      return;
    }
    const destination = this.destination!;
    const inner = new SwitchInner<R>(destination, this);
    this.inner = inner;
    this.add(inner);
    from(input).subscribe(inner);
  }

  /** @internal */
  innerDone(inner: SwitchInner<R>): void {
    if (this.inner === inner) this.inner = undefined;
    if (this.outerDone && this.inner === undefined) super._complete();
  }

  protected override _complete(): void {
    // The source is done; the stream completes when the current inner does.
    this.outerDone = true;
    this.isStopped = true;
    if (this.inner === undefined) super._complete();
  }

  protected override _dispose(): void {
    super._dispose();
    this.inner = undefined;
  }
}

/** Each value to an inner source (`project`), following only the latest; completes when the source and it have. */
export function switchMap<T, O extends ObservableInput<any>>(
  project: (value: T, index: number) => O,
): OperatorFunction<T, ObservedValueOf<O>> {
  return (source) =>
    lift(
      source,
      (destination: Subscriber<ObservedValueOf<O>>) =>
        new SwitchMapSubscriber<T, ObservedValueOf<O>>(destination, project as any),
    );
}
