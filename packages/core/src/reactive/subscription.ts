/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { SubscriptionLike, TeardownLogic, Unsubscribable } from './types.js';

type Finalizer = (() => void) | Unsubscribable;

/**
 * What ends a subscription: `unsubscribe` runs its own teardown (`_dispose`, for subclasses), then each finalizer
 * `add`ed, once. Finalizers and parents are held as `undefined`, one, or an array from two, so a subscription with
 * one teardown allocates no array (reactive §5).
 */
export class Subscription implements SubscriptionLike {
  static EMPTY: Subscription = (() => {
    const empty = new Subscription();
    empty.closed = true;
    return empty;
  })();

  closed = false;
  private _finalizers: Finalizer | Finalizer[] | undefined = undefined;
  /** The subscriptions this one was `add`ed to: it leaves them when it ends first (no growth of a long-lived parent). */
  private _parents: Subscription | Subscription[] | undefined = undefined;

  constructor(initialTeardown?: () => void) {
    if (initialTeardown) this._finalizers = initialTeardown;
  }

  /** A subclass's own teardown, run first. */
  protected _dispose(): void {}

  /** Whether the finalizers run last added first (`Scope`); in the order added otherwise (rxjs). */
  protected get _lastFirst(): boolean {
    return false;
  }

  /** What `unsubscribe` throws when finalizers threw (several: an `AggregateError`). */
  protected _failure(errors: any[]): any {
    return errors.length === 1 ? errors[0] : new AggregateError(errors, 'Errors while unsubscribing');
  }

  unsubscribe(): void {
    if (this.closed) return;
    this.closed = true;

    const parents = this._parents;
    if (parents !== undefined) {
      this._parents = undefined;
      if (Array.isArray(parents)) for (const p of parents) p.remove(this);
      else parents.remove(this);
    }

    let errors: any[] | undefined;
    try {
      this._dispose();
    } catch (e) {
      errors = [e];
    }

    const finalizers = this._finalizers;
    if (finalizers !== undefined) {
      this._finalizers = undefined;
      if (Array.isArray(finalizers)) {
        const lastFirst = this._lastFirst;
        for (let i = 0, n = finalizers.length; i < n; i++) {
          try {
            runFinalizer(finalizers[lastFirst ? n - 1 - i : i]);
          } catch (e) {
            (errors ??= []).push(e);
          }
        }
      } else {
        try {
          runFinalizer(finalizers);
        } catch (e) {
          (errors ??= []).push(e);
        }
      }
    }

    if (errors) throw this._failure(errors);
  }

  /** Runs `teardown` when this ends (at once if it has ended). A `Subscription` added leaves this when it ends first. */
  add(teardown: TeardownLogic): void {
    if (!teardown || teardown === this) return;
    if (this.closed) {
      runFinalizer(teardown);
      return;
    }
    if (teardown instanceof Subscription) {
      if (teardown.closed) return;
      teardown._addParent(this);
    }
    const f = this._finalizers;
    if (f === undefined) this._finalizers = teardown;
    else if (Array.isArray(f)) f.push(teardown);
    else this._finalizers = [f, teardown];
  }

  remove(teardown: Exclude<TeardownLogic, void>): void {
    const f = this._finalizers;
    if (f === teardown) this._finalizers = undefined;
    else if (Array.isArray(f)) {
      const i = f.indexOf(teardown);
      if (i >= 0) f.splice(i, 1);
    }
    if (teardown instanceof Subscription) teardown._removeParent(this);
  }

  private _addParent(parent: Subscription) {
    const p = this._parents;
    if (p === undefined) this._parents = parent;
    else if (Array.isArray(p)) p.push(parent);
    else this._parents = [p, parent];
  }

  private _removeParent(parent: Subscription) {
    const p = this._parents;
    if (p === parent) this._parents = undefined;
    else if (Array.isArray(p)) {
      const i = p.indexOf(parent);
      if (i >= 0) p.splice(i, 1);
    }
  }
}

function runFinalizer(f: Finalizer) {
  if (typeof f === 'function') f();
  else f.unsubscribe();
}

export function isSubscription(value: unknown): value is Subscription {
  return value instanceof Subscription;
}
