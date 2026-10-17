/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isObservable, type Observable, type Subscription, type Unsubscribable } from 'rxjs';
import { type Listenable, listen } from '../reactive/listen.js';
import type { ObservableState } from '../reactive/observable-state.js';
import type { ReadonlyState } from '../reactive/readonly-state.js';

/** What a scope stops and starts with its activity: a binding, a watch. */
abstract class Live implements Unsubscribable {
  abstract start(): void;
  abstract stop(): void;
  abstract unsubscribe(): void;
}

/**
 * A binding that can be stopped and started again with its scope's activity (`Scope.setActive`): stopped it holds
 * nothing (nothing runs for it), started it writes the current value and follows the changes. A state (or a derived
 * value) is followed by its `listen`, a plain callback; any other stream (the exception: rendering.md §2.2) by an rxjs
 * subscription, which replays its current value.
 */
class Binding<T> extends Live {
  #stop: (() => void) | Subscription | undefined;
  #live = false;
  #dead = false;
  constructor(
    readonly stream: Observable<T>,
    readonly fn: (value: T) => void,
  ) {
    super();
  }

  start() {
    if (this.#live || this.#dead) return;
    this.#live = true;
    const source = this.stream as Observable<T> & Partial<ReadonlyState<T>>;
    if (typeof source.listen === 'function') {
      // `fn` may stop the binding meanwhile.
      this.fn(source.value as T);
      if (this.#live) this.#stop = source.listen(this.fn);
      return;
    }
    // The stream replays its current value on subscribe (a state does), and `fn` may stop the binding meanwhile.
    const sub = this.stream.subscribe(this.fn);
    if (this.#live) this.#stop = sub;
    else sub.unsubscribe();
  }

  stop() {
    this.#live = false;
    const stop = this.#stop;
    this.#stop = undefined;
    if (typeof stop === 'function') stop();
    else stop?.unsubscribe();
  }

  unsubscribe() {
    this.#dead = true;
    this.stop();
  }
}

// --- watch ---------------------------------------------------------------------------------------

/** The watches marked by a change since the last flush (rendering.md §2.3): rendered together, in a microtask. */
let dirty = new Set<Watcher>();
let scheduled = false;
let flushing = false;
/** Flushes in a row before a watch that keeps marking itself (a render that sets what it watches) is dropped. */
const MaxPasses = 100;

function schedule(w: Watcher) {
  dirty.add(w);
  if (scheduled) return;
  scheduled = true;
  queueMicrotask(flush);
}

const byDepth = (a: Watcher, b: Watcher) => a.depth - b.depth;

/**
 * Renders now what `watch`es are waiting for (they render in a microtask otherwise): for code that sets state and then
 * reads the DOM in the same task (scrolling to a row just opened, focusing one, a test). Parents render before their
 * children, a disposed or paused watch not at all; renders that mark others are rendered in the same flush. Nested in a
 * flush (a render that calls it), it does nothing: the running flush goes on.
 */
export function flush(): void {
  scheduled = false;
  if (flushing) return;
  flushing = true;
  try {
    for (let pass = 0; dirty.size > 0; pass++) {
      const list = [...dirty];
      dirty = new Set();
      if (pass === MaxPasses) {
        console.error('watch: renders kept marking watches (a render sets what it watches?); dropped', list.length);
        break;
      }
      // A Set keeps the order of marking; the sort is stable, so siblings render in the order they changed.
      if (list.length > 1) list.sort(byDepth);
      for (let i = 0; i < list.length; i++) list[i].run();
    }
  } finally {
    flushing = false;
  }
}

/**
 * What `watch` makes: one listener per source (`listen`, no operators), `render` once when it starts and once per flush
 * after a change; stopped (its scope paused), it holds no listener.
 */
class Watcher extends Live {
  readonly #sources: readonly (Observable<unknown> | ReadonlyState<unknown> | Listenable)[];
  readonly #render: () => void;
  #stops: (() => void)[] | undefined;
  #live = false;
  #dead = false;
  constructor(
    readonly depth: number,
    sources: readonly (Observable<unknown> | ReadonlyState<unknown> | Listenable)[],
    render: () => void,
  ) {
    super();
    this.#sources = sources;
    this.#render = render;
  }

  readonly #mark = () => {
    if (this.#live) schedule(this);
  };

  start() {
    if (this.#live || this.#dead) return;
    this.#live = true;
    const sources = this.#sources;
    const stops: (() => void)[] = new Array(sources.length);
    for (let i = 0; i < sources.length; i++) stops[i] = listen(sources[i], this.#mark);
    this.#stops = stops;
    this.run();
  }

  /** Renders (outside any scope: a render writes, it makes nothing to own), unless stopped meanwhile. */
  run() {
    if (!this.#live) return;
    const previous = current;
    current = undefined;
    try {
      this.#render();
    } catch (e) {
      // One failing render must not keep the others in the flush from running; the watch goes on.
      console.error(e);
    } finally {
      current = previous;
    }
  }

  stop() {
    this.#live = false;
    dirty.delete(this);
    const stops = this.#stops;
    this.#stops = undefined;
    if (stops) for (let i = 0; i < stops.length; i++) stops[i]();
  }

  unsubscribe() {
    this.#dead = true;
    this.stop();
  }
}

/**
 * Owns what a piece of rendered DOM needs to clean up: subscriptions of its bindings, child scopes of its dynamic
 * regions, and anything registered with `onCleanup`. Disposing a scope disposes its children first, then runs its own
 * cleanups in reverse order.
 *
 * A scope is **active** or not (hidden UI does no UI updates, styling.md §5): while it is not, the bindings made in it
 * and in its children (`bind`, so every reactive prop and child) hold no subscription, and `onActive` hooks tell
 * imperative code to pause; when it is active again they subscribe once more, and a state replays its current value, so
 * the DOM catches up with one write per binding. Nothing is rebuilt: the DOM and the child scopes stay as they are.
 * The shell sets the scope of a hidden panel or a folded section inactive (`activity`).
 */
export class Scope {
  readonly #parent: Scope | undefined;
  /** How many scopes are above it: a watch's place in a flush (parents first). */
  readonly depth: number;
  // Allocated on first use: most scopes (a list's rows) own a few subscriptions and no children.
  #children: Set<Scope> | undefined;
  #cleanups: ((() => void) | Unsubscribable)[] | undefined;
  #disposed = false;
  /** Set while disposing: children then skip leaving the parent's set (it is being cleared anyway). */
  #disposing = false;
  /** What was asked of this scope (`setActive`), and what holds with its parents' (an inactive one pauses them all). */
  #own = true;
  #active = true;
  #hooks: ((active: boolean) => void)[] | undefined;

  constructor(parent?: Scope) {
    this.#parent = parent;
    this.depth = parent ? parent.depth + 1 : 0;
    this.#active = parent ? parent.#active : true;
  }

  /** Whether the scope and all its parents are active; bindings update the DOM only while it is. */
  get active(): boolean {
    return this.#active;
  }

  /**
   * Pauses the scope (false) or resumes it (true), with everything under it: bindings unsubscribe and subscribe again,
   * `onActive` hooks run. Idempotent; a child scope made while it is paused starts paused. Whatever decides the
   * hiding (a flag, the element's `hidden`) is bound outside the scope it pauses, or it would stop with it.
   */
  setActive(active: boolean): void {
    if (this.#own === active) return;
    this.#own = active;
    this.#sync();
  }

  #sync(): void {
    if (this.#disposed) return;
    const active = this.#own && (this.#parent ? this.#parent.#active : true);
    if (active === this.#active) return;
    this.#active = active;
    // Resuming: this scope's bindings, then its hooks, then the children (a list reconciles before its rows' own
    // bindings subscribe, so a row it removes costs nothing). Pausing: hooks first, so imperative code stops before the
    // bindings it may read.
    if (active) this.#bindings(true);
    this.#runHooks(active);
    if (!active) this.#bindings(false);
    const children = this.#children;
    if (children) for (const child of children) child.#sync();
  }

  #bindings(start: boolean): void {
    const cleanups = this.#cleanups;
    if (!cleanups) return;
    // By index: a binding that starts may add one (the list stays the same length otherwise).
    for (let i = 0; i < cleanups.length; i++) {
      const c = cleanups[i];
      if (c instanceof Live) {
        if (start) c.start();
        else c.stop();
      }
    }
  }

  #runHooks(active: boolean): void {
    const hooks = this.#hooks;
    if (!hooks) return;
    for (let i = 0; i < hooks.length; i++) {
      try {
        hooks[i](active);
      } catch (e) {
        console.error(e);
      }
    }
  }

  /**
   * Calls `fn` with whether the scope is active now, and on every change: for what bindings cannot pause, imperative
   * views (a ResizeObserver, a rAF loop, a canvas that redraws on a timer, a store subscription of a React view): start
   * it when active, stop it when not. Cleaning up for good is still `onCleanup`.
   */
  onActive(fn: (active: boolean) => void): void {
    if (this.#disposed) return;
    (this.#hooks ??= []).push(fn);
    try {
      fn(this.#active);
    } catch (e) {
      console.error(e);
    }
  }

  /**
   * The scope is active while `active` is true (a plain value is applied once; a stream's current value now, then each
   * change). Its own subscription never pauses.
   */
  follow(active: MaybeReactive<boolean>): void {
    const stream = toStream(active);
    if (!stream) this.setActive(!!active);
    else this.own(stream.subscribe((v) => this.setActive(!!v)));
  }

  /**
   * Subscribes `fn` to `stream` for as long as the scope is active (`bind` does this): the subscription goes while it
   * is paused and comes back, with the stream's current value, when it resumes.
   */
  bind<T>(stream: Observable<T>, fn: (value: T) => void): void {
    const binding = new Binding(stream, fn);
    this.own(binding);
    if (this.#active) binding.start();
  }

  /** `watch` in this scope. */
  watch(sources: readonly WatchSource[], render: () => void): void {
    const streams: (Observable<unknown> | ReadonlyState<unknown> | Listenable)[] = [];
    for (const s of sources) {
      const stream =
        toStream(s as MaybeReactive<unknown>) ??
        (typeof (s as Listenable).listen === 'function' ? (s as Listenable) : undefined);
      if (stream) streams.push(stream);
    }
    const watcher = new Watcher(this.depth, streams, render);
    this.own(watcher);
    if (this.#active) watcher.start();
  }

  get disposed(): boolean {
    return this.#disposed;
  }

  /** A scope disposed with this one (or on its own, earlier); active while this one is (and `setActive` on it allows). */
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
    this.#hooks = undefined;
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
 * Runs `fn` with whether the current scope is active now (it starts active), and again on every change: the hook for
 * what bindings cannot pause (a ResizeObserver, a rAF loop, a timer, a canvas that redraws on its own). Hidden UI does
 * no updates: start the work when active, stop it when not. See `Scope.onActive`.
 */
export function onActive(fn: (active: boolean) => void): void {
  requireScope('onActive').onActive(fn);
}

/**
 * Calls `fn` with the value, and again with each new value while the current scope lives and is active. A plain value
 * is applied once. A state applies its current value at once; the source should replay its current value on subscribe
 * (a state does): while the scope is paused (`Scope.setActive`) the subscription is dropped, and subscribing again on
 * resume is what brings the DOM up to date. A scope that is paused when the binding is made subscribes on its resume.
 */
export function bind<T>(source: MaybeReactive<T>, fn: (value: T) => void): void {
  const stream = toStream(source);
  if (!stream) {
    fn(source as T);
    return;
  }
  requireScope().bind(stream, fn);
}

/**
 * What a `watch` follows: a state (its `listen`), a writable state, anything else with `listen` (a keyed selector's
 * key), or any stream (subscribed).
 */
export type WatchSource = ReadonlyState<unknown> | ObservableState<unknown> | Observable<unknown> | Listenable;

/**
 * A view as a function of its state (rendering.md §2): `render` reads the current values (`state.value`) and writes
 * what differs from the DOM (`dom/writes.ts`). It runs now, and again after any of `sources` changes: the change marks
 * the watch, and every marked watch renders once in a microtask (`flush`), parents before children, so many changes in
 * a task are one render. While the scope is paused (`Scope.setActive`) it listens to nothing; resumed, it renders once.
 * A render owns nothing (it runs outside the scope: structure goes in regions, `Show`, `For`); one that throws is
 * logged and the watch goes on.
 */
export function watch(sources: readonly WatchSource[], render: () => void): void {
  requireScope('watch').watch(sources, render);
}
