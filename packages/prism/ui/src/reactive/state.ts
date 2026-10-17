/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Observable, type Subscriber } from 'rxjs';
import type { ReadonlyState } from './readonly-state.js';

// The reactive base (rendering.md §2.2, view-models.md §3): a state is its latest value plus a plain array of
// listeners. Nothing on its path is rxjs: `listen` is a callback, `derive` follows other states by `listen`. rxjs is the
// boundary: a state is also an `Observable` (`asObservable`, made by subscribing: one listener per subscriber), and
// `toState` turns a stream at a model's edge into a state, once.

/**
 * Counts for the benches (Internals › Performance › Trees, `__mspRows.updates`): listener calls (every state's and
 * derived value's), derived recomputes, `watch` renders. A plain increment each.
 */
export const ReactiveCounts = { listener: 0, derive: 0, render: 0 };

/**
 * For tests and benches (the waste audit, view-models.md §3): set `on`, and it is told of every state and derived value
 * made (`created`), every derived recompute and whether it changed (`derive`); `source` is the state whose change is
 * being told now (the outermost), for what follows from it to be put down to it. Unset, nothing is told.
 */
export const ReactiveProbe: {
  on: boolean;
  source: object | undefined;
  created?: (made: object, kind: 'state' | 'derive') => void;
  derive?: (d: object, changed: boolean) => void;
  /** A change told (the outermost), with the listener calls it took, its own and what followed. */
  told?: (state: object, listenerCalls: number) => void;
  /**
   * A listener added to `state` (the subscription inventory: what listens to what, view-models.md §3); what it returns
   * is given back to `unlistened` when that listener goes.
   */
  listened?: (state: object) => unknown;
  unlistened?: (token: unknown) => void;
} = { on: false, source: undefined };

/** A listener: told each value set, in order (one set inside a delivery is told after it). */
export type Listener<T> = (value: T) => void;

/**
 * The listeners of a state or a derived value: a plain array, removals leave a hole (cleared when the next delivery
 * ends), so a listener removed or added while others are told is safe.
 */
/** A listener's place in a list (`Listened`): emptied when it goes. */
interface Slot<T> {
  fn: Listener<T> | undefined;
}

export abstract class Listened<T> extends Observable<T> {
  /**
   * Each listener in a slot of its own: removed, its slot is emptied (no search of the list: a list of rows torn down
   * together removed each by a scan, N rows × L listeners), and the emptied slots are dropped once they are half the
   * list (and none is being told), so a removal costs O(1) amortized.
   */
  #listeners: Slot<T>[] | undefined;
  #holes = 0;
  #telling = 0;
  /** Values told while the listeners are told (a listener's set): told after, in order (made on first need). */
  #queued: T[] | undefined;

  abstract get value(): T;

  get asObservable(): Observable<T> {
    return this;
  }

  /** Whether anything listens (a derived value follows its sources only then). */
  protected get listened(): boolean {
    return !!this.#listeners && this.#listeners.length > this.#holes;
  }

  /** Calls `fn` with each new value until the returned function is called. */
  listen(fn: Listener<T>): () => void {
    const first = !this.listened;
    const slot: Slot<T> = { fn };
    (this.#listeners ??= []).push(slot);
    if (first) this.onFirst();
    const token = ReactiveProbe.on ? ReactiveProbe.listened?.(this) : undefined;
    return () => {
      if (!slot.fn) return;
      slot.fn = undefined;
      if (token !== undefined) ReactiveProbe.unlistened?.(token);
      this.#holes++;
      if (this.#telling === 0) this.#compact();
      if (!this.listened) this.onLast();
    };
  }

  /** Drops the emptied slots once they are half the list (all of them when none is left). */
  #compact() {
    const list = this.#listeners;
    if (!list || this.#holes === 0) return;
    if (this.#holes === list.length) {
      list.length = 0;
      this.#holes = 0;
      return;
    }
    if (this.#holes * 2 < list.length) return;
    let k = 0;
    for (let i = 0; i < list.length; i++) if (list[i].fn) list[k++] = list[i];
    list.length = k;
    this.#holes = 0;
  }

  /** Before the first listener (a derived value starts following its sources). */
  protected onFirst(): void {}
  /** After the last listener went. */
  protected onLast(): void {}

  /**
   * Tells every listener `value`. One told while the listeners are being told (a listener sets the state, a key set
   * again, a derived value recomputed) is queued and told after, so every listener gets the values in the order told and
   * ends on the latest (rendering.md §2.2).
   */
  protected tell(value: T): void {
    const list = this.#listeners;
    if (!list || list.length === this.#holes) return;
    if (this.#telling > 0) {
      (this.#queued ??= []).push(value);
      return;
    }
    // The probe: the outermost state told is what the rest follows from.
    if (ReactiveProbe.on && ReactiveProbe.source === undefined) {
      ReactiveProbe.source = this;
      const before = ReactiveCounts.listener;
      try {
        this.#deliver(list, value);
      } finally {
        ReactiveProbe.source = undefined;
        ReactiveProbe.told?.(this, ReactiveCounts.listener - before);
      }
      return;
    }
    this.#deliver(list, value);
  }

  #deliver(list: Slot<T>[], value: T): void {
    this.#tellAll(list, value);
    const q = this.#queued;
    if (!q) return;
    while (q.length > 0) this.#tellAll(list, q.shift() as T);
  }

  #tellAll(list: Slot<T>[], value: T): void {
    this.#telling++;
    try {
      // Those added meanwhile are told the next value, not this one.
      for (let i = 0, n = list.length; i < n; i++) {
        const fn = list[i].fn;
        if (!fn) continue;
        try {
          ReactiveCounts.listener++;
          fn(value);
        } catch (e) {
          console.error(e);
        }
      }
    } finally {
      if (--this.#telling === 0) this.#compact();
    }
  }

  /** The rxjs boundary: a subscriber gets the current value, then each new one (one listener per subscriber). */
  protected _subscribe(subscriber: Subscriber<T>) {
    return replayThenListen(this, subscriber);
  }
}

/**
 * A subscriber gets the current value, then each new one; one set while the current value is given (the subscriber
 * sets the state) is given after it, not missed.
 */
function replayThenListen<T>(state: ReadonlyState<T> & Pick<Listened<T>, 'listen'>, subscriber: Subscriber<T>) {
  let replaying = true;
  let missed = false;
  const stop = state.listen((v) => {
    if (replaying) missed = true;
    else subscriber.next(v);
  });
  const v = state.value;
  subscriber.next(v);
  replaying = false;
  if (missed && !Object.is(state.value, v)) subscriber.next(state.value);
  return stop;
}

// --- Batching ------------------------------------------------------------------------------------

let batchDepth = 0;
/** States set in the running batch, each with its delivery (in the order first set). */
const batched = new Map<object, () => void>();

/**
 * Runs `fn` with state changes held: each state set inside is told once, with its final value, when the outermost
 * batch ends (in the order first set); one set back to where it started tells nothing. Reads inside see the new
 * values. Synchronous: nothing is deferred past `fn`. For code that sets several states (or one several times) for one
 * change, so what follows them sees it once.
 */
export function batch<T>(fn: () => T): T {
  batchDepth++;
  try {
    return fn();
  } finally {
    if (--batchDepth === 0) {
      // A delivery may set states (and batch) again: drain until nothing is held.
      while (batched.size > 0) {
        const deliveries = [...batched.values()];
        batched.clear();
        for (const deliver of deliveries) deliver();
      }
    }
  }
}

// --- State ---------------------------------------------------------------------------------------

/**
 * A writable state. `set` compares with the latest value (`compare`), and `value` is always the latest set; a set while
 * the state's listeners are told is told after (`Listened.tell`), so every listener gets the values in the order they
 * were set and ends on the latest; in a `batch`, it is held and told once when the batch ends.
 */
export class State<T> extends Listened<T> implements ReadonlyState<T> {
  #value: T;
  readonly #compare: (a: T, b: T) => boolean;
  #held: { value: T } | undefined;

  constructor(initial: T, compare: (a: T, b: T) => boolean = Object.is) {
    super();
    this.#value = initial;
    this.#compare = compare;
    if (ReactiveProbe.on) ReactiveProbe.created?.(this, 'state');
  }

  /** The latest value set (in a batch, the one held), also while a set before it is still being told. */
  get value(): T {
    return this.#held ? this.#held.value : this.#value;
  }

  #view: ReadonlyState<T> | undefined;

  /** The state as its readers see it (made on first use): its value, `listen`, a stream; no `set`. */
  get asReadonly(): ReadonlyState<T> {
    return (this.#view ??= new View(this));
  }

  readonly #deliverHeld = () => {
    const h = this.#held;
    this.#held = undefined;
    if (h) this.set(h.value);
  };

  set(value: T): void {
    if (batchDepth > 0) {
      if (this.#held) this.#held.value = value;
      else if (!this.#compare(this.#value, value)) {
        this.#held = { value };
        batched.set(this, this.#deliverHeld);
      }
      return;
    }
    if (this.#compare(this.#value, value)) return;
    this.#value = value;
    this.tell(value);
  }
}

/** A record's state: `update` sets the fields given, when any differs. */
export class RecordState<T extends Record<string, unknown>> extends State<T> {
  update(partial: Partial<T>): void {
    const current = this.value as Record<string, unknown>;
    for (const key of Object.keys(partial)) {
      if (!Object.is(current[key], (partial as Record<string, unknown>)[key])) {
        this.set({ ...current, ...partial } as T);
        return;
      }
    }
  }
}

/** A state's readonly view: everything but `set`, delegated (it holds nothing of its own). */
class View<T> extends Observable<T> implements ReadonlyState<T> {
  readonly #state: State<T>;
  constructor(state: State<T>) {
    super();
    this.#state = state;
  }
  get value(): T {
    return this.#state.value;
  }
  get asObservable(): Observable<T> {
    return this;
  }
  listen(fn: Listener<T>): () => void {
    return this.#state.listen(fn);
  }
  protected _subscribe(subscriber: Subscriber<T>) {
    return replayThenListen(this, subscriber);
  }
}

// --- derive --------------------------------------------------------------------------------------

type Values<S extends readonly ReadonlyState<unknown>[]> = {
  [K in keyof S]: S[K] extends ReadonlyState<infer V> ? V : never;
};

/**
 * A value derived from states: lazy and hot. While nothing listens, `value` computes it from the sources' current
 * values and nothing is held; while something listens, it follows its sources (`listen`), keeps the value, and tells
 * its listeners only when it differs (`same`, `Object.is` by default). One derived value is shared by everything that
 * uses it. `derive(a, f)` replaces `a.pipe(map(f), distinctUntilChanged())`; `derive([a, b], f)` replaces
 * `combineLatest([a, b]).pipe(map(f))`.
 */
export function derive<S, T>(
  source: ReadonlyState<S>,
  fn: (value: S) => T,
  same?: (a: T, b: T) => boolean,
): ReadonlyState<T>;
export function derive<const S extends readonly ReadonlyState<unknown>[], T>(
  sources: S,
  fn: (...values: Values<S>) => T,
  same?: (a: T, b: T) => boolean,
): ReadonlyState<T>;
export function derive(
  sources: ReadonlyState<unknown> | readonly ReadonlyState<unknown>[],
  fn: (...values: any[]) => unknown,
  same: (a: unknown, b: unknown) => boolean = Object.is,
): ReadonlyState<unknown> {
  return Array.isArray(sources)
    ? new Derived(sources as readonly ReadonlyState<unknown>[], fn, same)
    : new Derived1(sources as ReadonlyState<unknown>, fn, same);
}

/** What `derive` follows its sources with: their `listen`, or (a stream without one) a subscription. */
function follow(source: ReadonlyState<unknown>, fn: () => void): () => void {
  if (typeof source.listen === 'function') return source.listen(fn);
  let subscribing = true;
  const sub = source.subscribe(() => {
    if (!subscribing) fn();
  });
  subscribing = false;
  return () => sub.unsubscribe();
}

/** One source (the common case: no array, no spread). */
class Derived1<S, T> extends Listened<T> implements ReadonlyState<T> {
  #cached!: T;
  /** The source's value the cache was computed from: a source whose value runs ahead of its telling is still read right. */
  #from: S | undefined;
  #stop: (() => void) | undefined;
  constructor(
    readonly source: ReadonlyState<S>,
    readonly fn: (value: S) => T,
    readonly same: (a: T, b: T) => boolean,
  ) {
    super();
    if (ReactiveProbe.on) ReactiveProbe.created?.(this, 'derive');
  }

  get value(): T {
    const v = this.source.value;
    if (this.#stop && Object.is(v, this.#from)) return this.#cached;
    return this.fn(v);
  }

  readonly #changed = () => {
    ReactiveCounts.derive++;
    const v = this.source.value;
    const next = this.fn(v);
    this.#from = v;
    const same = this.same(this.#cached, next);
    if (ReactiveProbe.on) ReactiveProbe.derive?.(this, !same);
    if (same) return;
    this.#cached = next;
    this.tell(next);
  };

  protected override onFirst() {
    const v = this.source.value;
    this.#from = v;
    this.#cached = this.fn(v);
    this.#stop = follow(this.source, this.#changed);
  }

  protected override onLast() {
    const stop = this.#stop;
    this.#stop = undefined;
    this.#from = undefined;
    stop?.();
  }
}

class Derived<T> extends Listened<T> implements ReadonlyState<T> {
  #cached!: T;
  #from: unknown[] = [];
  #stops: (() => void)[] | undefined;
  constructor(
    readonly sources: readonly ReadonlyState<unknown>[],
    readonly fn: (...values: unknown[]) => T,
    readonly same: (a: T, b: T) => boolean,
  ) {
    super();
    if (ReactiveProbe.on) ReactiveProbe.created?.(this, 'derive');
  }

  #compute(): T {
    const s = this.sources;
    switch (s.length) {
      case 1:
        return this.fn(s[0].value);
      case 2:
        return this.fn(s[0].value, s[1].value);
      case 3:
        return this.fn(s[0].value, s[1].value, s[2].value);
      default:
        return this.fn(...s.map((x) => x.value));
    }
  }

  #current(): boolean {
    const s = this.sources;
    for (let i = 0; i < s.length; i++) if (!Object.is(s[i].value, this.#from[i])) return false;
    return true;
  }

  #remember() {
    const s = this.sources;
    for (let i = 0; i < s.length; i++) this.#from[i] = s[i].value;
  }

  get value(): T {
    return this.#stops && this.#current() ? this.#cached : this.#compute();
  }

  readonly #changed = () => {
    ReactiveCounts.derive++;
    this.#remember();
    const next = this.#compute();
    const same = this.same(this.#cached, next);
    if (ReactiveProbe.on) ReactiveProbe.derive?.(this, !same);
    if (same) return;
    this.#cached = next;
    this.tell(next);
  };

  protected override onFirst() {
    this.#remember();
    this.#cached = this.#compute();
    this.#stops = this.sources.map((s) => follow(s, this.#changed));
  }

  protected override onLast() {
    const stops = this.#stops;
    this.#stops = undefined;
    if (stops) for (const stop of stops) stop();
  }
}

// --- toState -------------------------------------------------------------------------------------

/** A state fed by a stream, with what stops it. */
export type FedState<T> = ReadonlyState<T> & { dispose(): void };

/**
 * The boundary from rxjs: a state that holds `initial`, then each value `source` emits (one subscription, made now, at
 * a model's or a source's edge: Mol*'s events, a `switchMap` over models). Views get the state, never the stream.
 * `dispose` ends the subscription.
 */
export function toState<T>(source: Observable<T>, initial: T, compare?: (a: T, b: T) => boolean): FedState<T> {
  const state = new State(initial, compare);
  const sub = source.subscribe((v) => state.set(v));
  return Object.assign(state as ReadonlyState<T>, { dispose: () => sub.unsubscribe() });
}

/** Whether a value is a state (or derived): something with a current `value` and `listen`. */
export function isState(
  value: unknown,
): value is ReadonlyState<unknown> & { listen: (fn: Listener<unknown>) => () => void } {
  return value instanceof Listened;
}
