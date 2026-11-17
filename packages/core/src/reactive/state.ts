/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { reportUnhandledError } from './config.js';
import { Listened } from './listened.js';
import { from, Observable, Subscriber, toSubscriber } from './observable.js';
import { BehaviorSubject } from './subject.js';
import { isPlainRecord, sameShallow } from '../util/object.js';
import type { Subscription } from './subscription.js';
import type {
  Listenable,
  ObservableInput,
  PartialObserver,
  ReadonlyState,
  Subscribable,
  Unsubscribable,
} from './types.js';

// States (reactive §4): a state is its latest value plus listeners; `derive` is a value computed from states,
// shared, lazy, glitch-free. Push dirty, pull values: `set` bumps the state's version and marks what derives from it;
// the flush then tells each listened node whose value differs from what its listeners were last told, pulling derived
// values in dependency order on read, so a diamond computes and tells once. Sources are explicit (no tracking).

/** What a node's listeners were told before anything was. */
const NOTHING: any = Object.freeze({ nothing: true });

/** Bumped by every change anywhere: a cold derived value that saw this version returns its cache unchecked. */
let globalVersion = 0;
let batchDepth = 0;
let flushing = false;
/**
 * Moves at the end of each flush: the sets between two flushes are one change (a batch, or one set). A state remembers
 * its value and version at the first set of an epoch, so a set back to that value restores the version (reactive §4).
 */
let epoch = 0;
/** A state's versions come from one counter: a version names one value of one state, never reused for another. */
let nextVersion = 0;
/**
 * Listened nodes that may have changed and marked sinks, in the order marked; a reused array (no allocation per
 * change).
 */
const pending: (StateNode<any> | Sink | undefined)[] = [];
/** The queue's length: emptied by count (an array's `length = 0` frees its store, and the next push allocates one). */
let pendingCount = 0;
/** A flush that tells more than this many nodes is a cycle (a listener that sets its own source forever). */
const FlushLimit = 1_000_000;

function endBatch() {
  if (--batchDepth === 0) flush();
}

function flush() {
  if (flushing) return;
  flushing = true;
  let i = 0;
  try {
    for (; i < pendingCount; i++) {
      if (i === FlushLimit) {
        reportUnhandledError(new Error(`reactive: changes did not settle after ${FlushLimit} deliveries (a cycle?)`));
        break;
      }
      const node = pending[i]!;
      pending[i] = undefined;
      node._pending = false;
      if (node instanceof Sink) {
        try {
          node._flush();
        } catch (e) {
          reportUnhandledError(e);
        }
        continue;
      }
      if (!node._hasListeners) continue;
      let value: unknown;
      try {
        value = node.value;
      } catch (e) {
        reportUnhandledError(e);
        continue;
      }
      if (node._told !== NOTHING && node._equals(node._told, value)) continue;
      node._told = value;
      node._deliver(value);
    }
  } finally {
    // What a broken-off flush left (the cycle guard): unmarked and dropped.
    for (let j = i; j < pendingCount; j++) {
      pending[j]!._pending = false;
      pending[j] = undefined;
    }
    pendingCount = 0;
    flushing = false;
    epoch++;
  }
}

/**
 * Runs `fn` with deliveries held: each state changed inside is told once, its final value, when the outermost batch
 * ends; one set back to where it started tells nothing. Reads inside see the new values. Synchronous.
 */
export function batch<T>(fn: () => T): T {
  batchDepth++;
  try {
    return fn();
  } finally {
    endBatch();
  }
}

/**
 * @internal For node kinds outside this file (`keyed.ts`): something changed that `nodes` read; they are marked, and
 * outside a batch the change is told now.
 */
export function signalChange(nodes: readonly Derived<any>[] | Derived<any> | undefined): void {
  globalVersion++;
  if (nodes !== undefined) {
    if (Array.isArray(nodes)) for (let i = 0; i < nodes.length; i++) nodes[i]._mark();
    else (nodes as Derived<any>)._mark();
  }
  if (batchDepth === 0) flush();
}

/** What a node tells when it changes: a derived value, or a sink. */
export interface Dependent {
  _mark(): void;
}

/** A node of the state graph: a `State` or a derived value. */
/** @internal */
export abstract class StateNode<T> extends Listened<T> implements ReadonlyState<T> {
  /** @internal Bumped when the value changes (by `_equals`). */
  _version = 0;
  /** @internal Listened (1 while it has listeners) plus each derived value following it: hot while above 0. */
  _refs = 0;
  /** @internal The derived values following this one while they are hot. */
  _dependents: Dependent[] | undefined = undefined;
  /** @internal In the flush's queue. */
  _pending = false;
  /** @internal The value the listeners were last told. */
  _told: any = NOTHING;

  constructor(
    /** @internal */
    readonly _equals: (a: T, b: T) => boolean,
  ) {
    super();
  }

  abstract get value(): T;

  /**
   * Whether anything follows the node: a listener, a subscriber, or a dependent in the graph (a derived value, a sink:
   * `@molstar/ui/dom`'s bindings and watches), reactive §4.1.
   */
  override get observed(): boolean {
    return this._refs > 0;
  }

  /** @internal Whether it has listeners to tell (the flush's delivery). */
  get _hasListeners(): boolean {
    return this.listenerCount > 0;
  }

  /** @internal Brings the value up to date and returns its version. */
  abstract _pull(): number;

  /** @internal Started being followed (listened or derived from). */
  protected _hot(): void {}
  /** @internal Stopped being followed. */
  protected _cold(): void {}

  protected override onFirst(): void {
    if (this._refs++ === 0) this._hot();
    try {
      this._told = this.value;
    } catch {
      this._told = NOTHING;
    }
  }

  protected override onLast(): void {
    this._told = NOTHING;
    if (--this._refs === 0) this._cold();
  }

  /** @internal */
  _addDependent(d: Dependent): void {
    (this._dependents ??= []).push(d);
    if (this._refs++ === 0) this._hot();
  }

  /** @internal */
  _removeDependent(d: Dependent): void {
    const ds = this._dependents!;
    const i = ds.lastIndexOf(d);
    if (i >= 0) {
      ds[i] = ds[ds.length - 1];
      ds.pop();
    }
    if (--this._refs === 0) this._cold();
  }

  /** @internal Queues this for the flush when something listens. */
  protected _enqueue(): void {
    if (!this._pending && this.listenerCount > 0) {
      this._pending = true;
      pending[pendingCount++] = this;
    }
  }

  protected _markDependents(): void {
    const ds = this._dependents;
    if (ds !== undefined) for (let i = 0; i < ds.length; i++) ds[i]._mark();
  }

  /** @internal */
  _deliver(value: T): void {
    this.emit(value);
  }

  /**
   * The current value, then each change. A set made while the current value is given (the subscriber sets the state)
   * is told after it, not nested.
   */
  override subscribe(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): Subscription {
    const subscriber = toSubscriber(observerOrNext);
    if (!this._add(subscriber)) return subscriber;
    batchDepth++;
    try {
      let value: T;
      try {
        value = this.value;
      } catch (e) {
        subscriber.error(e);
        return subscriber;
      }
      subscriber.next(value);
    } finally {
      endBatch();
    }
    return subscriber;
  }
}

/**
 * A writable state: `set` compares with the latest value (`equals`, `Object.is` by default) and, outside a batch,
 * tells the listeners at once; a set while listeners are told is queued and told after, in order, so every listener
 * ends on the latest. Never ends.
 */
export class State<T> extends StateNode<T> {
  private _value: T;
  private _view: StateView<T> | undefined = undefined;
  private _baseEpoch = -1;
  private _baseValue: T | undefined = undefined;
  private _baseVersion = 0;

  constructor(initial: T, equals: (a: T, b: T) => boolean = Object.is) {
    super(equals);
    this._value = initial;
  }

  /** The latest value set (also inside a batch, and while an earlier set is still being told). */
  get value(): T {
    return this._value;
  }

  /** The state as its readers see it: everything but `set` (made on first use). */
  get asReadonly(): ReadonlyState<T> {
    return (this._view ??= new StateView(this));
  }

  set(value: T): void {
    if (this._equals(this._value, value)) return;
    // Only a held change (a batch, a flush running) can be set back before it is told: remember where it started.
    if (batchDepth > 0 || flushing) {
      if (this._baseEpoch !== epoch) {
        this._baseEpoch = epoch;
        this._baseValue = this._value;
        this._baseVersion = this._version;
      }
    } else if (this._baseEpoch !== -1) {
      this._baseEpoch = -1;
      this._baseValue = undefined;
    }
    this._value = value;
    // Set back to where it was at the start of this change: its version too, so what compares versions (a derived
    // value's stamps, a binding, a watch) sees no change.
    this._version =
      this._baseEpoch === epoch && this._equals(this._baseValue as T, value) ? this._baseVersion : ++nextVersion;
    globalVersion++;
    this._enqueue();
    this._markDependents();
    if (batchDepth === 0) flush();
  }

  /** @internal */
  _pull(): number {
    return this._version;
  }
}

/**
 * A record's state: `update(partial)` sets the fields given when any differs (`Object.is`), building on the latest value.
 * Compared by `sameShallow` by default, so a set of an equal copy tells nothing.
 */
export class RecordState<T extends Record<string, unknown>> extends State<T> {
  constructor(initial: T, equals: (a: T, b: T) => boolean = sameShallow) {
    super(initial, equals);
  }

  update(partial: Partial<T>): void {
    const current = this.value as Record<string, unknown>;
    for (const key in partial) {
      if (!Object.is(current[key], (partial as Record<string, unknown>)[key])) {
        this.set({ ...current, ...partial } as T);
        return;
      }
    }
  }
}

/** A state's readonly view: delegates everything (holds nothing of its own). */
class StateView<T> extends Observable<T> implements ReadonlyState<T> {
  constructor(
    /** @internal */
    readonly _state: State<T>,
  ) {
    super();
  }

  get value(): T {
    return this._state.value;
  }

  listen(fn: (value: T) => void): () => void {
    return this._state.listen(fn);
  }

  override subscribe(observerOrNext?: PartialObserver<T> | ((value: T) => void) | null): Subscription {
    return this._state.subscribe(observerOrNext);
  }
}

/** A graph node for `source`, or the source itself when it is foreign (a `BehaviorSubject`, another state). */
function unwrap(source: ReadonlyState<unknown>): ReadonlyState<unknown> {
  return source instanceof StateView ? source._state : source;
}

/** What a derived value compares to know a source changed: a node's version, a foreign source's value. */
function stampOf(source: ReadonlyState<unknown>): unknown {
  return source instanceof StateNode ? source._pull() : source.value;
}

/** Tells a derived value that a foreign source changed. */
class ForeignSubscriber extends Subscriber<unknown> {
  constructor(private readonly derived: Derived<any>) {
    super();
  }

  protected override _next(): void {
    globalVersion++;
    this.derived._mark();
    if (batchDepth === 0) flush();
  }

  // A foreign source that ends leaves the derived value at its last value.
  protected override _error(): void {
    this.unsubscribe();
  }

  protected override _complete(): void {
    this.unsubscribe();
  }
}

function followForeign(source: ReadonlyState<unknown>, derived: Derived<any>): Unsubscribable {
  if (source instanceof Listened) {
    // Its listener list directly: no replay of a `BehaviorSubject`'s value.
    const subscriber = new ForeignSubscriber(derived);
    source._add(subscriber);
    return subscriber;
  }
  if (typeof source.listen === 'function') {
    const stop = source.listen(() => {
      globalVersion++;
      derived._mark();
      if (batchDepth === 0) flush();
    });
    return { unsubscribe: stop };
  }
  // A stream with a value but no `listen` (an rxjs-based state): what it replays while subscribing is no change.
  const subscriber = new ForeignSubscriber(derived);
  subscriber.isStopped = true;
  const sub = source.subscribe(subscriber);
  subscriber.isStopped = false;
  return sub;
}

/**
 * A value derived from explicit sources. Lazy: unlistened, `value` computes from the sources (cached until anything
 * changes) and nothing follows them. Hot while listened or derived from: it follows its sources, recomputes when read
 * after one changed, and its version moves only when the result differs (`equals`), so what derives from it stops
 * there.
 */
/** @internal The base of computed nodes (`derive`, and `keyed.ts`'s keys). */
export abstract class Derived<T> extends StateNode<T> {
  private _dirty = true;
  private _valid = false;
  private _checkedAt = -1;
  private _markedAt = -1;
  private _cached: T = undefined as T;
  private _foreign: Unsubscribable[] | undefined = undefined;
  protected _hasForeign: boolean;

  constructor(
    protected readonly _sources: readonly ReadonlyState<unknown>[],
    equals: (a: T, b: T) => boolean,
  ) {
    super(equals);
    let foreign = false;
    for (let i = 0; i < _sources.length; i++) if (!(_sources[i] instanceof StateNode)) foreign = true;
    this._hasForeign = foreign;
  }

  get value(): T {
    this._ensure();
    return this._cached;
  }

  /** @internal */
  _pull(): number {
    this._ensure();
    return this._version;
  }

  /** Whether a source changed since the last compute (records the sources' stamps). */
  protected abstract _updateStamps(valid: boolean): boolean;
  protected abstract _compute(): T;

  private _ensure() {
    if (this._refs > 0) {
      if (!this._dirty) return;
    } else if (!this._hasForeign && this._valid && this._checkedAt === globalVersion) {
      return;
    }
    const changed = this._updateStamps(this._valid);
    if (changed || !this._valid) {
      let next: T;
      try {
        next = this._compute();
      } catch (e) {
        this._valid = false;
        throw e;
      }
      if (!this._valid || !this._equals(this._cached, next)) {
        this._cached = next;
        this._version++;
      }
      this._valid = true;
    }
    this._dirty = false;
    this._checkedAt = globalVersion;
  }

  /** @internal A source changed: possibly stale, queued if listened, and so is what derives from this. */
  _mark(): void {
    if (this._markedAt === globalVersion) return;
    this._markedAt = globalVersion;
    this._dirty = true;
    this._enqueue();
    this._markDependents();
  }

  protected override _hot(): void {
    const sources = this._sources;
    for (let i = 0; i < sources.length; i++) {
      const s = sources[i];
      if (s instanceof StateNode) s._addDependent(this);
      else (this._foreign ??= []).push(followForeign(s, this));
    }
    this._dirty = true;
  }

  protected override _cold(): void {
    const sources = this._sources;
    for (let i = 0; i < sources.length; i++) {
      const s = sources[i];
      if (s instanceof StateNode) s._removeDependent(this);
    }
    const foreign = this._foreign;
    if (foreign !== undefined) {
      for (const f of foreign) f.unsubscribe();
      foreign.length = 0;
    }
    this._checkedAt = -1;
  }
}

class Derived1<S, T> extends Derived<T> {
  private _stamp: unknown = NOTHING;
  private readonly _source: ReadonlyState<S>;

  constructor(
    source: ReadonlyState<S>,
    private readonly _fn: (value: S) => T,
    equals: (a: T, b: T) => boolean,
  ) {
    super([source], equals);
    this._source = source;
  }

  protected _updateStamps(valid: boolean): boolean {
    const stamp = stampOf(this._source);
    if (valid && Object.is(stamp, this._stamp)) return false;
    this._stamp = stamp;
    return true;
  }

  protected _compute(): T {
    return this._fn(this._source.value);
  }
}

class DerivedN<T> extends Derived<T> {
  private readonly _stamps: unknown[];

  constructor(
    sources: readonly ReadonlyState<unknown>[],
    private readonly _fn: (...values: any[]) => T,
    equals: (a: T, b: T) => boolean,
  ) {
    super(sources, equals);
    this._stamps = new Array(sources.length).fill(NOTHING);
  }

  protected _updateStamps(valid: boolean): boolean {
    const sources = this._sources;
    const stamps = this._stamps;
    let changed = !valid;
    for (let i = 0; i < sources.length; i++) {
      const stamp = stampOf(sources[i]);
      if (!Object.is(stamp, stamps[i])) {
        stamps[i] = stamp;
        changed = true;
      }
    }
    return changed;
  }

  protected _compute(): T {
    const s = this._sources;
    // No spread for the common counts.
    switch (s.length) {
      case 0:
        return this._fn();
      case 1:
        return this._fn(s[0].value);
      case 2:
        return this._fn(s[0].value, s[1].value);
      case 3:
        return this._fn(s[0].value, s[1].value, s[2].value);
      case 4:
        return this._fn(s[0].value, s[1].value, s[2].value, s[3].value);
      default: {
        const values = new Array(s.length);
        for (let i = 0; i < s.length; i++) values[i] = s[i].value;
        return this._fn(...values);
      }
    }
  }
}

/** Named sources: `fn` gets a new record of their values each compute (it may keep it). */
class DerivedRecord<T> extends DerivedN<T> {
  constructor(
    private readonly _keys: readonly string[],
    sources: readonly ReadonlyState<unknown>[],
    private readonly _project: (values: Record<string, unknown>) => T,
    equals: (a: T, b: T) => boolean,
  ) {
    super(sources, () => undefined as T, equals);
  }

  protected override _compute(): T {
    const values: Record<string, unknown> = {};
    const keys = this._keys;
    const sources = this._sources;
    for (let i = 0; i < keys.length; i++) values[keys[i]] = sources[i].value;
    return this._project(values);
  }
}

type RecordValues<S extends Record<string, ReadonlyState<unknown>>> = {
  [K in keyof S]: S[K] extends ReadonlyState<infer V> ? V : never;
};

type Values<S extends readonly ReadonlyState<unknown>[]> = {
  [K in keyof S]: S[K] extends ReadonlyState<infer V> ? V : never;
};

/**
 * A value derived from states (or any `ReadonlyState`: a `BehaviorSubject`): lazy and hot, shared by everything that
 * reads it, told only when it differs (`equals`, `Object.is` by default). `derive(a, f)` replaces
 * `a.pipe(map(f), distinctUntilChanged(), shareReplay(1))`; `derive([a, b], f)` replaces
 * `combineLatest([a, b]).pipe(map(f))`; `derive({ a, b }, ({ a, b }) => …)` takes named sources.
 */
export function derive<S, T>(
  source: ReadonlyState<S>,
  fn: (value: S) => T,
  equals?: (a: T, b: T) => boolean,
): ReadonlyState<T>;
export function derive<const S extends readonly ReadonlyState<unknown>[], T>(
  sources: S,
  fn: (...values: Values<S>) => T,
  equals?: (a: T, b: T) => boolean,
): ReadonlyState<T>;
export function derive<const S extends Record<string, ReadonlyState<unknown>>, T>(
  sources: S,
  fn: (values: RecordValues<S>) => T,
  equals?: (a: T, b: T) => boolean,
): ReadonlyState<T>;
export function derive(
  sources: ReadonlyState<unknown> | readonly ReadonlyState<unknown>[] | Record<string, ReadonlyState<unknown>>,
  fn: (...values: any[]) => unknown,
  equals: (a: unknown, b: unknown) => boolean = Object.is,
): ReadonlyState<unknown> {
  if (Array.isArray(sources)) return new DerivedN(sources.map(unwrap), fn, equals);
  if (isPlainRecord(sources)) {
    const keys = Object.keys(sources);
    const list = keys.map((k) => unwrap((sources as Record<string, ReadonlyState<unknown>>)[k]));
    return new DerivedRecord(keys, list, fn, equals);
  }
  return new Derived1(unwrap(sources as ReadonlyState<unknown>), fn, equals);
}

/**
 * @internal The third kind of graph node (reactive §4.1): sources and no value, for what renders or writes from states
 * (`@molstar/ui/dom`'s `watch` and bindings). It is a dependent of each graph source it follows (a slot in the source's array,
 * no object of its own), marked with them, and run in the flush (`_flush`), after the change has spread, so what it
 * reads is consistent; once per flush however many of its sources changed.
 */
export abstract class Sink implements Dependent {
  /** @internal */
  _pending = false;
  private _markedAt = -1;

  /** @internal A source changed. */
  _mark(): void {
    if (this._markedAt === globalVersion) return;
    this._markedAt = globalVersion;
    if (this._pending) return;
    this._pending = true;
    pending[pendingCount++] = this;
  }

  /** @internal Run in the flush. */
  abstract _flush(): void;
}

/**
 * @internal The graph node behind `source` (a `State`, a derived value, a keyed value, a state's readonly view), or
 * undefined for anything else (a `BehaviorSubject`, a foreign stream).
 */
export function graphNode(source: unknown): StateNode<any> | undefined {
  const s = source instanceof StateView ? source._state : source;
  return s instanceof StateNode ? s : undefined;
}

/** @internal Makes `dependent` a follower of `node` (it is marked with it); the node is hot while followed. */
export function follow(node: StateNode<any>, dependent: Dependent): void {
  node._addDependent(dependent);
}

/** @internal Undoes `follow`. */
export function unfollow(node: StateNode<any>, dependent: Dependent): void {
  node._removeDependent(dependent);
}

/** @internal A node's version (brought up to date: a derived value computes if it must). */
export function versionOf(node: StateNode<any>): number {
  return node._pull();
}

/** A state fed by a stream, with what stops it. */
export type FedState<T> = ReadonlyState<T> & { dispose(): void };

class Fed<T> extends State<T> {
  /** @internal */
  _subscription: Unsubscribable | undefined = undefined;

  dispose(): void {
    const s = this._subscription;
    this._subscription = undefined;
    s?.unsubscribe();
  }
}

/**
 * A state that holds `initial`, then each value `source` tells (one subscription, made now): a stream with no current
 * value at a model's edge. A `BehaviorSubject` needs none: `derive` follows it directly.
 */
export function toState<T>(
  source: Observable<T> | Subscribable<T>,
  initial: T,
  equals?: (a: T, b: T) => boolean,
): FedState<T> {
  const state = new Fed(initial, equals);
  state._subscription = (from(source as ObservableInput<T>) as Observable<T>).subscribe((v) => state.set(v));
  return state;
}

/** Whether `value` is a state: a `State`, a derived value, a readonly view, or a `BehaviorSubject`. */
export function isState(value: unknown): value is ReadonlyState<unknown> {
  return value instanceof StateNode || value instanceof StateView || value instanceof BehaviorSubject;
}

/**
 * Calls `fn` after `source` changes, until the returned function is called: its `listen` when it has one, else a
 * subscription, where what it tells while subscribing (a replayed current value) is no change.
 */
export function listen(
  source: Listenable<any> | Observable<unknown> | Subscribable<unknown>,
  fn: () => void,
): () => void {
  const own = (source as Listenable<unknown>).listen;
  if (typeof own === 'function') return own.call(source, () => fn());
  let subscribing = true;
  const sub = (source as Subscribable<unknown>).subscribe({
    next: () => {
      if (!subscribing) fn();
    },
  });
  subscribing = false;
  return () => sub.unsubscribe();
}
