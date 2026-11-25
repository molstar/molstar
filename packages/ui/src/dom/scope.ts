/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { reportUnhandledError } from '@molstar/core/reactive/config';
import { Scope } from '@molstar/core/reactive/model';
import { isObservable, type Observable } from '@molstar/core/reactive/observable';
import type { ReadonlyState, Unsubscribable } from '@molstar/core/reactive/types';

/** A value, or something the DOM follows: any `Observable` (a state is one), or a writable state's `asReadonly`. */
export type MaybeReactive<T> = T | Observable<T> | ReadonlyState<T> | { readonly asReadonly: ReadonlyState<T> };

/** The stream behind a reactive value, or undefined for a plain value. */
export function toStream<T>(value: MaybeReactive<T>): Observable<T> | undefined {
  if (isObservable(value)) return value as Observable<T>;
  if (value && typeof value === 'object' && 'asReadonly' in value) {
    const view = (value as { asReadonly: unknown }).asReadonly;
    if (isObservable(view)) return view as Observable<T>;
  }
  return undefined;
}

/** @internal What a scope stops and starts with its activity: a binding, a watch. */
export interface Live {
  _start(): void;
  _stop(): void;
}

/** What `dispose` throws when cleanups failed: they were reported one by one. */
const Reported = {};

/**
 * The registry of root scopes a test setup puts on `globalThis` (`vitest.setup.ts`: it disposes the roots a test file
 * leaves, so that a view left listening on the document cannot wake up in a later file of the worker); absent in apps,
 * where this is one lookup per root scope
 */
const TestRootScopes = Symbol.for('molstar.test.rootScopes');
function testRootScopes(): { add(scope: { dispose(): void }): void; delete(scope: object): void } | undefined {
  return (globalThis as { [TestRootScopes]?: ReturnType<typeof testRootScopes> })[TestRootScopes];
}

/**
 * What a piece of rendered DOM owns (dom §1): its bindings and watches, the child scopes of its regions, its cleanups.
 * Core's `Scope` (cleanups last added first, every one attempted) plus a parent, a depth (a watch's place in a flush:
 * parents first), children (an intrusive list: no allocation per child, O(1) removal) and **activity**: a paused scope
 * (`setActive(false)`, or a paused parent) holds no listener for its bindings and watches; resumed, they catch up once
 * (each binding writes its source's current value, each watch renders). Nothing is rebuilt. Disposing disposes the
 * children first (newest first), then the cleanups; a failing cleanup is reported and the rest run.
 */
export class ViewScope extends Scope {
  readonly depth: number;
  private readonly _parent: ViewScope | undefined;
  private _firstChild: ViewScope | undefined = undefined;
  private _lastChild: ViewScope | undefined = undefined;
  private _prev: ViewScope | undefined = undefined;
  private _next: ViewScope | undefined = undefined;
  /** What `setActive` asked, and what holds with the parents' (an inactive one pauses them all). */
  private _asked = true;
  private _active: boolean;
  /** The parts that pause (bindings, watches), apart from what only ends, so a pause visits only them. */
  private _live: Live[] | undefined = undefined;
  private _hooks: ((active: boolean) => void)[] | undefined = undefined;

  constructor(parent?: ViewScope) {
    super();
    this._parent = parent;
    this.depth = parent ? parent.depth + 1 : 0;
    this._active = parent ? parent._active : true;
    if (!parent) testRootScopes()?.add(this);
  }

  /** Whether the scope and all its parents are active; bindings update the DOM only while it is. */
  get active(): boolean {
    return this._active;
  }

  get disposed(): boolean {
    return this.closed;
  }

  /**
   * Pauses the scope (false) or resumes it (true), with everything under it. Idempotent; a child made while it is
   * paused starts paused. What decides the hiding is bound outside the scope it pauses, or it would stop with it.
   */
  setActive(active: boolean): void {
    if (this._asked === active) return;
    this._asked = active;
    this._sync();
  }

  private _sync(): void {
    if (this.closed) return;
    const active = this._asked && (this._parent ? this._parent._active : true);
    if (active === this._active) return;
    this._active = active;
    // Resuming: this scope's bindings, then its hooks, then the children (a list reconciles before its rows' bindings
    // start, so a row it removes costs nothing). Pausing: hooks first, so imperative code stops before what it reads.
    if (active) this._lives(true);
    this._runHooks(active);
    if (!active) this._lives(false);
    // The next child is read after each one syncs, never saved before (dom §1): a child's resume can flush a change that
    // makes this scope's region dispose a sibling, or the child itself (then the walk starts over: those synced already
    // return at once).
    for (let c = this._firstChild; c !== undefined; ) {
      c._sync();
      // A nested `setActive` flipped this scope meanwhile: that call synced the children.
      if (this._active !== active) return;
      c = c.closed ? this._firstChild : c._next;
    }
  }

  private _lives(start: boolean): void {
    const live = this._live;
    if (live === undefined) return;
    // By index: a binding that starts may add one.
    for (let i = 0; i < live.length; i++) {
      if (start) live[i]._start();
      else live[i]._stop();
    }
  }

  private _runHooks(active: boolean): void {
    const hooks = this._hooks;
    if (hooks === undefined) return;
    for (let i = 0; i < hooks.length; i++) {
      try {
        hooks[i](active);
      } catch (e) {
        reportUnhandledError(e);
      }
    }
  }

  /**
   * Calls `fn` with whether the scope is active now, and on every change: for what bindings cannot pause (an observer, a
   * frame loop, a timer). Cleaning up for good is still `onCleanup`.
   */
  onActive(fn: (active: boolean) => void): void {
    if (this.closed) return;
    (this._hooks ??= []).push(fn);
    try {
      fn(this._active);
    } catch (e) {
      reportUnhandledError(e);
    }
  }

  /** The scope is active while `active` is true (a stream: its current value, then each change; never paused). */
  follow(active: MaybeReactive<boolean>): void {
    const stream = toStream(active);
    if (!stream) this.setActive(!!active);
    else this.add(stream.subscribe((v) => this.setActive(!!v)));
  }

  /** @internal A part that pauses with the scope (and ends with it); started now if the scope is active. */
  _addLive(live: Live & Unsubscribable): void {
    if (this.closed) {
      live.unsubscribe();
      return;
    }
    (this._live ??= []).push(live);
    this.add(live);
    if (this._active) live._start();
  }

  /** Runs `fn` when the scope is disposed; at once when it is already. */
  onCleanup(fn: () => void): void {
    if (this.closed) {
      runCleanup(fn);
      return;
    }
    this.add(fn);
  }

  /** A scope disposed with this one (or on its own, earlier); active while this one is. */
  child(): ViewScope {
    const scope = new ViewScope(this);
    if (this.closed) {
      scope.dispose();
      return scope;
    }
    scope._prev = this._lastChild;
    if (this._lastChild) this._lastChild._next = scope;
    else this._firstChild = scope;
    this._lastChild = scope;
    return scope;
  }

  override dispose(): void {
    try {
      this.unsubscribe();
    } catch (e) {
      if (e !== Reported) throw e;
    }
  }

  protected override _dispose(): void {
    if (!this._parent) testRootScopes()?.delete(this);
    this._hooks = undefined;
    // Newest first, each child unlinked before it is disposed and the list read again after it (dom §1): this scope's
    // bindings are still live while its children go, so a child's cleanup can make one of its regions dispose another
    // child; a saved link would end the walk there.
    for (let c = this._lastChild; c !== undefined; c = this._lastChild) {
      this._unlink(c);
      c.dispose();
    }
    this._live = undefined;
    this._parent?._unlink(this);
  }

  /** Takes `c` out of the children (nothing if it is not among them). */
  private _unlink(c: ViewScope): void {
    if (c._prev === undefined && this._firstChild !== c) return;
    if (c._prev) c._prev._next = c._next;
    else this._firstChild = c._next;
    if (c._next) c._next._prev = c._prev;
    else this._lastChild = c._prev;
    c._prev = c._next = undefined;
  }

  protected override _failure(errors: any[]): any {
    for (const e of errors) reportUnhandledError(e);
    return Reported;
  }
}

function runCleanup(fn: () => void) {
  try {
    fn();
  } catch (e) {
    reportUnhandledError(e);
  }
}

let current: ViewScope | undefined;

/** The scope that owns what is rendered now, if any. */
export function getScope(): ViewScope | undefined {
  return current;
}

/** The scope that owns what is rendered now; rendering outside one is a mistake (its bindings would leak). */
export function requireScope(what = 'A binding'): ViewScope {
  if (!current) throw new Error(`${what} was created outside a render scope; render with \`render\` or in an element.`);
  return current;
}

/** Runs `fn` with `scope` owning what it renders (`undefined`: nothing owns it, as a `watch` render). */
export function runInScope<T>(scope: ViewScope | undefined, fn: () => T): T {
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
  const scope = new ViewScope();
  return runInScope(scope, () => fn(() => scope.dispose()));
}

/**
 * Runs `fn` with whether the current scope is active now, and again on every change: the hook for what bindings cannot
 * pause. See `ViewScope.onActive`.
 */
export function onActive(fn: (active: boolean) => void): void {
  requireScope('onActive').onActive(fn);
}
