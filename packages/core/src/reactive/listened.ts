/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ConsumerSubscriber, Observable, type Subscriber, type SubscriberOwner, toSubscriber } from './observable.js';
import type { Subscription } from './subscription.js';
import type { PartialObserver } from './types.js';

function noop() {}

/**
 * The base of everything with listeners (reactive §2): `Subject`, `BehaviorSubject`, `State`, derived values,
 * `share`, and any source a package writes. Each listener is a `Subscriber` in a slot of a plain array: the subscriber
 * is the slot and the subscription at once, so `subscribe(fn)` allocates one object; removing it empties its slot (no
 * search), and the emptied slots are dropped once they are half the list and no delivery runs (O(1) amortized: rows
 * torn down together). The list is made on first use.
 */
export abstract class Listened<T> extends Observable<T> implements SubscriberOwner {
  private _subs: (Subscriber<T> | undefined)[] | undefined = undefined;
  private _count = 0;
  private _holes = 0;
  private _depth = 0;

  /** Whether anything listens. */
  get observed(): boolean {
    return this._count > 0;
  }

  /** The number of listeners. */
  protected get listenerCount(): number {
    return this._count;
  }

  /** The first listener came: start what feeds this (an upstream subscription, a DOM listener). */
  protected onFirst(): void {}

  /** The last listener went: stop it. */
  protected onLast(): void {}

  override subscribe(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): Subscription {
    const subscriber = toSubscriber(observerOrNext);
    this._add(subscriber);
    return subscriber;
  }

  /**
   * Calls `fn` with each value until the returned function is called: a plain callback, no replay of a current value.
   * What `fn` throws is reported (`reportUnhandledError`).
   */
  listen(fn: (value: T) => void): () => void {
    const subscriber = new ConsumerSubscriber(fn);
    if (!this._add(subscriber)) return noop;
    return () => subscriber.unsubscribe();
  }

  /** @internal Adds a listener; false when it was already closed. */
  _add(subscriber: Subscriber<T>): boolean {
    if (subscriber.closed) return false;
    const subs = (this._subs ??= []);
    subscriber._owner = this;
    subscriber._index = subs.length;
    subs.push(subscriber);
    if (++this._count === 1) this.onFirst();
    return true;
  }

  /** @internal */
  _remove(subscriber: Subscriber<any>): void {
    const subs = this._subs;
    if (subs === undefined || subscriber._owner !== this) return;
    subs[subscriber._index] = undefined;
    subscriber._owner = undefined;
    subscriber._index = -1;
    this._holes++;
    if (--this._count === 0) {
      if (this._depth === 0) this._clear();
      this.onLast();
    } else if (this._depth === 0 && this._holes * 2 > subs.length) {
      this._compact();
    }
  }

  /**
   * Tells every listener `value` at once, nested: a value emitted by a listener reaches the listeners before the
   * remaining ones get this one (rxjs's `Subject`). A listener added while it runs is not told this value; one removed
   * and not yet told is not told it. No allocation.
   */
  protected emit(value: T): void {
    const subs = this._subs;
    if (subs === undefined || this._count === 0) return;
    this._depth++;
    try {
      const n = subs.length;
      for (let i = 0; i < n; i++) {
        const s = subs[i];
        if (s !== undefined) s.next(value);
      }
    } finally {
      if (--this._depth === 0) this._settle();
    }
  }

  /** Tells every listener `err`; each ends (and leaves the list). */
  protected emitError(err: any): void {
    const subs = this._subs;
    if (subs === undefined) return;
    this._depth++;
    try {
      const n = subs.length;
      for (let i = 0; i < n; i++) {
        const s = subs[i];
        if (s !== undefined) s.error(err);
      }
    } finally {
      if (--this._depth === 0) this._settle();
    }
  }

  /** Tells every listener it is complete; each ends (and leaves the list). */
  protected emitComplete(): void {
    const subs = this._subs;
    if (subs === undefined) return;
    this._depth++;
    try {
      const n = subs.length;
      for (let i = 0; i < n; i++) {
        const s = subs[i];
        if (s !== undefined) s.complete();
      }
    } finally {
      if (--this._depth === 0) this._settle();
    }
  }

  /** Drops every listener without telling it (each subscription ends). */
  protected removeAll(): void {
    const subs = this._subs;
    if (subs === undefined) return;
    // As a delivery: no compaction while the loop walks the list, so none is skipped.
    this._depth++;
    try {
      const n = subs.length;
      for (let i = 0; i < n; i++) subs[i]?.unsubscribe();
    } finally {
      if (--this._depth === 0) this._settle();
    }
  }

  private _settle() {
    if (this._count === 0) this._clear();
    else if (this._holes * 2 > this._subs!.length) this._compact();
  }

  private _clear() {
    this._subs!.length = 0;
    this._holes = 0;
  }

  private _compact() {
    const subs = this._subs!;
    let j = 0;
    for (let i = 0; i < subs.length; i++) {
      const s = subs[i];
      if (s === undefined) continue;
      s._index = j;
      subs[j++] = s;
    }
    subs.length = j;
    this._holes = 0;
  }
}
