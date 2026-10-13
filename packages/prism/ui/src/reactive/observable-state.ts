/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject } from 'rxjs';

import { isPlainRecord, sameShallow } from './same-shallow.js';
import type { ReadonlyState } from './readonly-state.js';

export type ObservableStateOptions<T> = {
  readonly compare?: (left: T, right: T) => boolean;
};

export type ObservableState<T> = {
  readonly value: T;
  readonly asReadonly: ReadonlyState<T>;
  set(value: T): void;
};

export type ObservableRecordState<T extends Record<string, unknown>> = ObservableState<T> & {
  update(partial: Partial<T>): void;
};

function hasPartialChange<T extends Record<string, unknown>>(current: T, partial: Partial<T>): boolean {
  for (const key of Object.keys(partial) as (keyof T)[]) {
    if (!Object.is(current[key], partial[key])) return true;
  }
  return false;
}

// --- Batching ------------------------------------------------------------------------------------

let batchDepth = 0;
/** States set in the running batch, each with its delivery (in the order first set). */
const batched = new Map<object, () => void>();

/**
 * Runs `fn` with state changes held: each state set inside emits once, with its final value, when the outermost batch
 * ends (in the order first set); one set back to where it started emits nothing. Reads inside see the new values.
 * Synchronous: nothing is deferred past `fn`. For code that sets several states (or one several times) for one change,
 * so views render it once.
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

export function observableState<T extends Record<string, unknown>>(
  initial: T,
  options?: ObservableStateOptions<T>,
): ObservableRecordState<T>;
export function observableState<T>(initial: T, options?: ObservableStateOptions<T>): ObservableState<T>;
export function observableState<T>(initial: T, options?: ObservableStateOptions<T>): ObservableState<T> {
  const subject = new BehaviorSubject(initial);
  const asReadonly = subject.asObservable() as ReadonlyState<T>;
  Object.defineProperty(asReadonly, 'value', { get: () => state.value });
  const compare = options?.compare ?? sameShallow;
  // A set during a delivery (from a subscriber) is queued and delivered after
  // it, so every subscriber gets the values in the order they were set and ends on the latest; a
  // set compares with the latest value, queued or delivered.
  const queued: T[] = [];
  let delivering = false;
  // Held by a batch: the value set last, delivered when the batch ends.
  let held: { value: T } | undefined;
  const latest = () => (held ? held.value : queued.length > 0 ? (queued[queued.length - 1] as T) : subject.getValue());
  const deliverHeld = () => {
    const h = held;
    held = undefined;
    if (h) state.set(h.value);
  };
  const state: ObservableState<T> = {
    get value() {
      return held ? held.value : subject.getValue();
    },
    asReadonly,
    set(value) {
      if (batchDepth > 0) {
        if (held) held.value = value;
        else if (!compare(latest(), value)) {
          held = { value };
          batched.set(state, deliverHeld);
        }
        return;
      }
      if (compare(latest(), value)) return;
      if (delivering) {
        queued.push(value);
        return;
      }
      delivering = true;
      try {
        subject.next(value);
        while (queued.length > 0) subject.next(queued.shift() as T);
      } finally {
        delivering = false;
        queued.length = 0;
      }
    },
  };
  if (!isPlainRecord(initial)) return state;
  return Object.assign(state, {
    update(partial: Partial<T>) {
      const current = latest();
      if (!isPlainRecord(current) || !hasPartialChange(current, partial as Partial<Record<string, unknown>>)) return;
      state.set({ ...current, ...partial } as T);
    },
  });
}
