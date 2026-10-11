/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isObservable, type Observable, type Unsubscribable } from 'rxjs';
import type { ObservableState } from '../reactive/observable-state.js';

/**
 * Owns what a piece of rendered DOM needs to clean up: subscriptions of its bindings, child scopes of its dynamic
 * regions, and anything registered with `onCleanup`. Disposing a scope disposes its children first, then runs its own
 * cleanups in reverse order.
 */
export class Scope {
  readonly #parent: Scope | undefined;
  // Allocated on first use: most scopes (a list's rows) own a few subscriptions and no children.
  #children: Set<Scope> | undefined;
  #cleanups: ((() => void) | Unsubscribable)[] | undefined;
  #disposed = false;
  /** Set while disposing: children then skip leaving the parent's set (it is being cleared anyway). */
  #disposing = false;

  constructor(parent?: Scope) {
    this.#parent = parent;
  }

  get disposed(): boolean {
    return this.#disposed;
  }

  /** A scope disposed with this one (or on its own, earlier). */
  child(): Scope {
    const scope = new Scope(this);
    if (this.#disposed) scope.dispose();
    else (this.#children ??= new Set()).add(scope);
    return scope;
  }

  /** Runs `fn` when the scope is disposed; at once when it already is. */
  onCleanup(fn: () => void): void {
    if (this.#disposed) {
      runCleanup(fn);
      return;
    }
    (this.#cleanups ??= []).push(fn);
  }

  /** Unsubscribes `subscription` when the scope is disposed (no closure per binding). */
  own(subscription: Unsubscribable): void {
    if (this.#disposed) {
      subscription.unsubscribe();
      return;
    }
    (this.#cleanups ??= []).push(subscription);
  }

  dispose(): void {
    if (this.#disposed) return;
    this.#disposed = true;
    this.#disposing = true;
    const children = this.#children;
    if (children) {
      // Newest first; a Set iterates in insertion order, so walk a copy backwards.
      const list = [...children];
      for (let i = list.length - 1; i >= 0; i--) list[i].dispose();
      this.#children = undefined;
    }
    const cleanups = this.#cleanups;
    this.#cleanups = undefined;
    if (cleanups) {
      for (let i = cleanups.length - 1; i >= 0; i--) {
        const c = cleanups[i];
        if (typeof c === 'function') runCleanup(c);
        else c.unsubscribe();
      }
    }
    const parent = this.#parent;
    if (parent && !parent.#disposing) parent.#children?.delete(this);
  }
}

function runCleanup(fn: () => void) {
  try {
    fn();
  } catch (e) {
    // A failing cleanup must not keep the others from running.
    console.error(e);
  }
}

let current: Scope | undefined;

/** The scope that owns what is rendered now, if any. */
export function getScope(): Scope | undefined {
  return current;
}

/** The scope that owns what is rendered now; rendering outside one is a mistake (its bindings would leak). */
export function requireScope(what = 'A binding'): Scope {
  if (!current) throw new Error(`${what} was created outside a render scope; render with \`render\` or in an element.`);
  return current;
}

/** Runs `fn` with `scope` owning what it renders. */
export function runInScope<T>(scope: Scope, fn: () => T): T {
  const previous = current;
  current = scope;
  try {
    return fn();
  } finally {
    current = previous;
  }
}

/** Runs `fn` when the current scope is disposed. */
export function onCleanup(fn: () => void): void {
  requireScope('onCleanup').onCleanup(fn);
}

/** Runs `fn` in a new root scope (owned by no other); it gets the scope's dispose. */
export function createRoot<T>(fn: (dispose: () => void) => T): T {
  const scope = new Scope();
  return runInScope(scope, () => fn(() => scope.dispose()));
}

/** A value, or a stream of values the DOM follows: an observable, a `ReadonlyState`, or a writable `observableState`. */
export type MaybeReactive<T> = T | Observable<T> | ObservableState<T>;

/** The stream behind a reactive value (a writable state's readonly view), or undefined for a plain value. */
export function toStream<T>(value: MaybeReactive<T>): Observable<T> | undefined {
  if (isObservable(value)) return value;
  if (
    value &&
    typeof value === 'object' &&
    'asReadonly' in value &&
    isObservable((value as ObservableState<T>).asReadonly)
  )
    return (value as ObservableState<T>).asReadonly;
  return undefined;
}

/**
 * Calls `fn` with the value, and again with each new value while the current scope lives. A plain value is applied
 * once. A state applies its current value at once.
 */
export function bind<T>(source: MaybeReactive<T>, fn: (value: T) => void): void {
  const stream = toStream(source);
  if (!stream) {
    fn(source as T);
    return;
  }
  const scope = requireScope();
  scope.own(stream.subscribe(fn));
}
