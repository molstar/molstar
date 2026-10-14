/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// State primitives (seqspec:viewer#state): states with a synchronous value, live records, events. Emission is a plain
// loop over subscribers and allocates nothing; rxjs 7's `Subject.next` is not used.

import { Observable, type Subscriber, type TeardownLogic } from 'rxjs';

/** Who made a change, carried so a host skips its own echoes. */
export type Origin = 'pointer' | 'keyboard' | 'api' | 'behavior' | 'graph' | 'snapshot' | `host:${string}`;

/** An rxjs `Observable` with a synchronous `value` (structurally any-ui's ReadonlyState); subscribing replays it. */
export interface SeqState<T> extends Observable<T> {
  readonly value: T;
  /** The origin of the value being delivered (read it in a subscriber), else of the last set. */
  readonly origin: Origin | undefined;
}
export interface WritableSeqState<T> extends SeqState<T> {
  set(v: T, origin?: Origin): void;
  update(partial: Partial<T>): void;
  readonly asReadonly: SeqState<T>;
  /** Completes every subscriber; later sets do nothing. */
  complete(): void;
}

/** A reused mutable record owned by its manager, and a version that bumps on change; for pointer- and frame-rate values. */
export interface Live<T> {
  readonly current: T;
  readonly version: SeqState<number>;
}
export interface WritableLive<T> extends Live<T> {
  bump(): void;
  complete(): void;
}

export interface Unsubscribable {
  unsubscribe(): void;
}
/** A happening, not state: nothing replays; subscribers by priority (higher first, then subscription order). */
export interface SeqEvent<T> {
  subscribe(fn: (e: T) => void, opts?: { priority?: number }): Unsubscribable;
}
export interface SeqEventEmitter<T> extends SeqEvent<T> {
  emit(e: T): void;
  complete(): void;
  readonly subscriberCount: number;
}

function isPlain(x: unknown): x is Record<string, unknown> {
  if (x === null || typeof x !== 'object') return false;
  const p = Object.getPrototypeOf(x);
  return p === Object.prototype || p === null;
}

/** Identity, and shallow for plain records (edge-rate: UI-rate states hold immutable values). */
function defaultEquals(a: unknown, b: unknown): boolean {
  if (a === b) return true;
  if (!isPlain(a) || !isPlain(b)) return false;
  let n = 0;
  for (const k in a) {
    if (!Object.hasOwn(b, k) || a[k] !== b[k]) return false;
    n++;
  }
  for (const k in b) if (Object.hasOwn(b, k)) n--;
  return n === 0;
}

const EMPTY: readonly never[] = [];

class StateImpl<T> extends Observable<T> implements WritableSeqState<T> {
  private latest: T;
  private lastOrigin: Origin | undefined = undefined;
  private readonly equals: (a: T, b: T) => boolean;
  /** Copy on write: delivery iterates the array it started with. */
  private subs: readonly Subscriber<T>[] = EMPTY;
  private delivering = false;
  /** The value being delivered: a subscriber added meanwhile replays it, then gets the queued ones. */
  private delivered: T;
  private closed = false;
  /** Sets made during delivery (U017): a ring of values and origins, grown only when full. */
  private qValues: T[] = [];
  private qOrigins: (Origin | undefined)[] = [];
  private qHead = 0;
  private qSize = 0;
  /** A value made on demand (`pulled`): called before the value is read or replayed; it sets the value. */
  private readonly pull: (() => void) | null;

  constructor(initial: T, equals?: (a: T, b: T) => boolean, pull?: () => void) {
    // rxjs calls this with the observable as `this` on every subscribe.
    super(function (this: Observable<T>, subscriber: Subscriber<T>) {
      return (this as StateImpl<T>).attach(subscriber);
    });
    this.latest = initial;
    this.delivered = initial;
    this.equals = equals ?? defaultEquals;
    this.pull = pull ?? null;
  }

  get value() {
    this.pull?.();
    return this.latest;
  }
  get origin() {
    this.pull?.();
    return this.lastOrigin;
  }
  get asReadonly(): SeqState<T> {
    return this;
  }
  get subscriberCount() {
    return this.subs.length;
  }

  /** Replays the current value, then delivers until unsubscribed. */
  private attach(subscriber: Subscriber<T>): TeardownLogic {
    if (this.closed) {
      subscriber.complete();
      return;
    }
    if (!this.delivering) this.pull?.();
    this.subs = [...this.subs, subscriber];
    subscriber.next(this.delivering ? this.delivered : this.latest);
    return () => {
      this.subs = this.subs.filter((s) => s !== subscriber);
    };
  }

  set(v: T, origin?: Origin) {
    if (this.closed || this.equals(this.latest, v)) return;
    this.latest = v;
    if (this.delivering) {
      this.enqueue(v, origin);
      return;
    }
    this.delivering = true;
    try {
      this.deliver(v, origin);
      while (this.qSize > 0) {
        const i = this.qHead;
        const value = this.qValues[i];
        const o = this.qOrigins[i];
        this.qValues[i] = undefined as T;
        this.qHead = (i + 1) % this.qValues.length;
        this.qSize--;
        this.deliver(value, o);
      }
    } finally {
      this.delivering = false;
      this.delivered = undefined as T;
      this.qHead = 0;
      this.qSize = 0;
    }
  }

  update(partial: Partial<T>) {
    this.set({ ...this.latest, ...partial });
  }

  complete() {
    if (this.closed) return;
    this.closed = true;
    const subs = this.subs;
    this.subs = EMPTY;
    for (let i = 0; i < subs.length; i++) subs[i].complete();
  }

  private deliver(v: T, origin: Origin | undefined) {
    this.lastOrigin = origin;
    this.delivered = v;
    const subs = this.subs;
    // A closed subscriber is skipped: rxjs allocates a notification for a `next` after unsubscribe.
    for (let i = 0; i < subs.length; i++) if (!subs[i].closed) subs[i].next(v);
  }

  private enqueue(v: T, origin: Origin | undefined) {
    const cap = this.qValues.length;
    if (this.qSize === cap) {
      const next = Math.max(4, cap * 2);
      const values = new Array<T>(next);
      const origins = new Array<Origin | undefined>(next);
      for (let k = 0; k < this.qSize; k++) {
        values[k] = this.qValues[(this.qHead + k) % cap];
        origins[k] = this.qOrigins[(this.qHead + k) % cap];
      }
      this.qValues = values;
      this.qOrigins = origins;
      this.qHead = 0;
    }
    const at = (this.qHead + this.qSize) % this.qValues.length;
    this.qValues[at] = v;
    this.qOrigins[at] = origin;
    this.qSize++;
  }
}

/** A state: `set` compares (identity, shallow for plain records) and emits only on change. */
export function state<T>(initial: T, compare?: (a: T, b: T) => boolean): WritableSeqState<T> {
  return new StateImpl(initial, compare);
}

/**
 * A state whose owner sets values on demand: `pull` is called before the value or origin is read and before a subscriber
 * is replayed, and sets the pending value, if any; the owner sets at once while `subscriberCount` is above 0. A change
 * that nobody reads then costs no value (a viewport's `committed` per programmatic move).
 */
export function pulled<T>(
  initial: T,
  pull: () => void,
  compare?: (a: T, b: T) => boolean,
): WritableSeqState<T> & { readonly subscriberCount: number } {
  return new StateImpl(initial, compare, pull);
}

/** A live record: `current` is mutated in place by its owner, then `bump()`. */
export function live<T>(current: T): WritableLive<T> {
  const version = new StateImpl(0);
  return {
    current,
    version,
    bump: () => version.set(version.value + 1),
    complete: () => version.complete(),
  };
}

interface EventSub<T> {
  readonly fn: (e: T) => void;
  readonly priority: number;
  active: boolean;
}

class EventImpl<T> implements SeqEventEmitter<T> {
  private subs: readonly EventSub<T>[] = EMPTY;
  private closed = false;
  private readonly onError: ((error: unknown) => void) | undefined;

  constructor(onError?: (error: unknown) => void) {
    this.onError = onError;
  }

  get subscriberCount() {
    return this.subs.length;
  }

  subscribe(fn: (e: T) => void, opts?: { priority?: number }): Unsubscribable {
    const sub: EventSub<T> = { fn, priority: opts?.priority ?? 0, active: !this.closed };
    if (this.closed) return { unsubscribe() {} };
    // Stable: after every subscriber of the same or a higher priority.
    let at = this.subs.length;
    while (at > 0 && this.subs[at - 1].priority < sub.priority) at--;
    this.subs = [...this.subs.slice(0, at), sub, ...this.subs.slice(at)];
    return {
      unsubscribe: () => {
        if (!sub.active) return;
        sub.active = false;
        this.subs = this.subs.filter((s) => s !== sub);
      },
    };
  }

  emit(e: T) {
    const subs = this.subs;
    let failed = false;
    let error: unknown;
    for (let i = 0; i < subs.length; i++) {
      const s = subs[i];
      if (!s.active) continue;
      try {
        s.fn(e);
      } catch (err) {
        if (this.onError) this.onError(err);
        else if (!failed) {
          failed = true;
          error = err;
        }
      }
    }
    if (failed) throw error;
  }

  complete() {
    this.closed = true;
    for (const s of this.subs) s.active = false;
    this.subs = EMPTY;
  }
}

/**
 * An event. The emitter reuses its record, so an emission allocates nothing. A subscriber that throws does not stop the
 * others: `onError` gets it, else the first error is rethrown after every subscriber ran.
 */
export function event<T>(onError?: (error: unknown) => void): SeqEventEmitter<T> {
  return new EventImpl<T>(onError);
}
