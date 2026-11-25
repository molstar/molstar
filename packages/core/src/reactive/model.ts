/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Listened } from './listened.js';
import { ConsumerSubscriber, Observable } from './observable.js';
import { Subscription } from './subscription.js';
import type { Listenable, Subscribable, Unsubscribable } from './types.js';

/** Something that ends with `dispose`. */
export interface Disposable {
  dispose(): void;
}

/** What a scope owns: a subscription (ours or any `Unsubscribable`), a `Disposable`, or a teardown function. */
export type Owned = Unsubscribable | Disposable | (() => void);

/**
 * Owns teardowns and ends them together, last added first (later resources may depend on earlier ones); every one is
 * attempted and the first failure rethrown (reactive §7). A `Subscription` added leaves the scope when it ends first, so
 * a scope does not grow with subscriptions that end on their own. A `Subscription`: its storage holds one teardown
 * without an array.
 */
export class Scope extends Subscription implements Disposable {
  protected override get _lastFirst(): boolean {
    return true;
  }

  protected override _failure(errors: any[]): any {
    return errors[0];
  }

  /** Owns `resource`: ended with this scope (at once if it has ended). */
  own(resource: Owned): void {
    if (typeof resource === 'function' || typeof (resource as Unsubscribable).unsubscribe === 'function') {
      this.add(resource as Unsubscribable | (() => void));
    } else {
      this.add(() => (resource as Disposable).dispose());
    }
  }

  dispose(): void {
    this.unsubscribe();
  }
}

const Mounting = 1;
const Mounted = 2;
const Retired = 4;
const UnmountAfterMount = 8;
const WarnedTwoOwners = 16;

/** Subclasses whose lifecycle methods were checked (once per class, not per instance). */
const checkedClasses = new WeakSet<object>();

/**
 * A model with a lifecycle (reactive §7): `mount` (its owner shows it), `unmount` (hidden; it may mount again),
 * `retire` (its final owner is done; for good). What `onMount` subscribes to, listens to or owns ends with the mount;
 * what is owned permanently ends with `retire`. A subclass implements the hooks (`onMount`, `onUnmount`, `onRetire`),
 * never `mount`, `unmount`, `retire`.
 */
export abstract class ReactiveModel {
  private _flags = 0;
  private _mountScope: Scope | undefined = undefined;
  private _permanentScope: Scope | undefined = undefined;
  private _mountController: AbortController | undefined = undefined;

  constructor() {
    const cls = this.constructor;
    if (!checkedClasses.has(cls)) {
      const base = ReactiveModel.prototype;
      if (this.mount !== base.mount || this.unmount !== base.unmount || this.retire !== base.retire) {
        throw new Error('ReactiveModel lifecycle methods cannot be overridden; use the protected lifecycle hooks');
      }
      checkedClasses.add(cls);
    }
  }

  mount(): void {
    const flags = this._flags;
    if (flags & Retired) return;
    if (flags & Mounted) {
      if (!(flags & WarnedTwoOwners)) {
        this._flags |= WarnedTwoOwners;
        console.warn(
          `${this.constructor.name}: mounted again while mounted. A model has one owner: the first unmount will tear it down under the other.`,
        );
      }
      return;
    }
    const scope = new Scope();
    this._mountScope = scope;
    this._flags |= Mounted | Mounting;
    try {
      this.onMount();
    } catch (error) {
      this._flags &= ~(Mounted | UnmountAfterMount | Mounting);
      if (this._mountScope === scope) this._mountScope = undefined;
      try {
        scope.dispose();
      } catch {
        // The mount failed first: its cleanup is best effort and must not hide that failure.
      }
      this._abortMount();
      throw error;
    }
    this._flags &= ~Mounting;
    if (this._flags & UnmountAfterMount) {
      this._flags &= ~UnmountAfterMount;
      this.unmount();
    }
  }

  unmount(): void {
    const flags = this._flags;
    if (!(flags & Mounted)) return;
    if (flags & Mounting) {
      // Carried out when `onMount` returns.
      this._flags |= UnmountAfterMount;
      return;
    }
    this._flags &= ~Mounted;
    const scope = this._mountScope;
    this._mountScope = undefined;
    let failed = false;
    let failure: unknown;
    try {
      this.onUnmount();
    } catch (error) {
      failed = true;
      failure = error;
    }
    try {
      scope?.dispose();
    } catch (error) {
      if (!failed) {
        failed = true;
        failure = error;
      }
    }
    this._abortMount();
    if (failed) throw failure;
  }

  /** Ends this model for good. Only its final owner calls this; a view that only hides it calls `unmount`. */
  retire(): void {
    if (this._flags & Retired) return;
    this._flags |= Retired;
    let failed = false;
    let failure: unknown;
    try {
      this.unmount();
    } catch (error) {
      failed = true;
      failure = error;
    }
    try {
      this.onRetire();
    } catch (error) {
      if (!failed) {
        failed = true;
        failure = error;
      }
    }
    try {
      this._permanentScope?.dispose();
    } catch (error) {
      if (!failed) {
        failed = true;
        failure = error;
      }
    }
    if (failed) throw failure;
  }

  /** Whether it is mounted now. */
  protected get isMounted(): boolean {
    return (this._flags & Mounted) !== 0;
  }

  /** Whether it was retired. */
  protected get isRetired(): boolean {
    return (this._flags & Retired) !== 0;
  }

  /**
   * Aborted when this mount ends (`unmount`, `retire`, a throwing `onMount`); a remount gets a fresh one. For work that
   * outlives a call but not the mount: fetches, timers, async loops. Read outside a mount, it is already aborted.
   */
  protected get mountSignal(): AbortSignal {
    if (!(this._flags & Mounted)) return AbortSignal.abort();
    return (this._mountController ??= new AbortController()).signal;
  }

  /** Runs `dispose` when this mount ends; the returned function runs it now instead (once). Only in `onMount`. */
  protected registerDispose(dispose: () => void): () => void {
    return removable(this._duringMount('register dispose'), dispose);
  }

  /** Runs `dispose` on `retire`; the returned function runs it now instead (once). */
  protected registerPermanentDispose(dispose: () => void): () => void {
    return removable(this._permanent('register permanent dispose'), dispose);
  }

  /** Owns `resource` until this mount ends. Only in `onMount`. */
  protected own(resource: Owned): void {
    this._duringMount('own a resource').own(resource);
  }

  /** Owns `resource` until `retire`. */
  protected ownPermanently(resource: Owned): void {
    this._permanent('own permanently').own(resource);
  }

  /**
   * Calls `action` with each value of `source` until this mount ends (a state or a `BehaviorSubject` gives its current
   * value first). One object per subscription for a stream of ours. Only in `onMount`.
   */
  protected subscribe<T>(source: Observable<T> | Subscribable<T>, action: (value: T) => void): void {
    const scope = this._duringMount('subscribe');
    scope.add(subscribeTo(source, action));
  }

  /** Calls `action` with each value of `source` until `retire`. */
  protected subscribePermanently<T>(source: Observable<T> | Subscribable<T>, action: (value: T) => void): void {
    const scope = this._permanent('subscribe permanently');
    // Added after subscribing: a synchronous value that retires this model ends the subscription at once (the scope
    // has ended, `add` runs it).
    scope.add(subscribeTo(source, action));
  }

  /**
   * Calls `fn` with each change of `source` (a state, a subject, a keyed value) until this mount ends: no current value
   * first, where `subscribe` would give one. Only in `onMount`.
   */
  protected listen<T>(source: Listenable<T>, fn: (value: T) => void): void {
    this._duringMount('listen').add(listenTo(source, fn));
  }

  /** Calls `fn` with each change of `source` until `retire`. */
  protected listenPermanently<T>(source: Listenable<T>, fn: (value: T) => void): void {
    this._permanent('listen permanently').add(listenTo(source, fn));
  }

  protected onMount(): void {}

  protected onUnmount(): void {}

  protected onRetire(): void {}

  private _duringMount(what: string): Scope {
    const scope = this._mountScope;
    if (!(this._flags & Mounting) || scope === undefined) throw new Error(`Cannot ${what} outside onMount`);
    return scope;
  }

  private _permanent(what: string): Scope {
    if (this._flags & Retired) throw new Error(`Cannot ${what} after retirement`);
    return (this._permanentScope ??= new Scope());
  }

  private _abortMount() {
    const controller = this._mountController;
    this._mountController = undefined;
    controller?.abort();
  }
}

// The lifecycle methods cannot be replaced on an instance either (an assignment throws in strict code).
for (const name of ['mount', 'unmount', 'retire'] as const) {
  Object.defineProperty(ReactiveModel.prototype, name, { writable: false, configurable: false });
}

function removable(scope: Scope, dispose: () => void): () => void {
  scope.add(dispose);
  let active = true;
  return () => {
    if (!active) return;
    active = false;
    // Already run when the scope ended.
    if (scope.closed) return;
    scope.remove(dispose);
    dispose();
  };
}

function subscribeTo<T>(source: Observable<T> | Subscribable<T>, action: (value: T) => void): Unsubscribable {
  // Ours takes the function (one subscriber object); a foreign `Subscribable` takes an observer.
  return source instanceof Observable ? source.subscribe(action) : source.subscribe({ next: action });
}

function listenTo<T>(source: Listenable<T>, fn: (value: T) => void): Unsubscribable | (() => void) {
  if (source instanceof Listened) {
    // Its listener list directly: the subscriber is the listener and the subscription (one object).
    const subscriber = new ConsumerSubscriber(fn);
    source._add(subscriber);
    return subscriber;
  }
  return source.listen(fn);
}
