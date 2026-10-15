/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The one debounce of the viewer's forms (a color picker or a slider fires a change per pointer move; one edit is one
// history step, so it waits for the next one). DOM-free: plain timers. A view releases it with `dispose` (or `cancel`)
// when it goes, so no edit lands after its form is gone.

/** The delay before a form edit is applied (ms): short enough to follow a drag, long enough to make it one step. */
export const EditDelay = 150;

/** Delayed runs by key: scheduling a key again restarts its wait, and each key runs once its wait is over. */
export class Debouncer<K = void> {
  readonly #timers = new Map<K, ReturnType<typeof setTimeout>>();

  constructor(
    private readonly run: (key: K) => void,
    readonly delay: number = EditDelay,
  ) {}

  /** Runs `key` after the delay, from now (an earlier wait of the key is dropped). */
  schedule(key: K): void {
    clearTimeout(this.#timers.get(key));
    this.#timers.set(
      key,
      setTimeout(() => {
        this.#timers.delete(key);
        this.run(key);
      }, this.delay),
    );
  }

  /** Whether `key` waits to run. */
  has(key: K): boolean {
    return this.#timers.has(key);
  }

  /** Whether anything waits. */
  get waiting(): boolean {
    return this.#timers.size > 0;
  }

  /** Runs `key` now if it waits. */
  flush(key: K): void {
    if (!this.#timers.has(key)) return;
    clearTimeout(this.#timers.get(key));
    this.#timers.delete(key);
    this.run(key);
  }

  /** Runs everything that waits now. */
  flushAll(): void {
    for (const key of [...this.#timers.keys()]) this.flush(key);
  }

  /** Drops `key`'s wait. */
  cancel(key: K): void {
    clearTimeout(this.#timers.get(key));
    this.#timers.delete(key);
  }

  /** The view went: nothing runs after this. */
  dispose(): void {
    for (const t of this.#timers.values()) clearTimeout(t);
    this.#timers.clear();
  }
}

/** A function that runs once the calls stop for `delay` ms, with the last call's arguments. */
export interface Debounced<A extends unknown[]> {
  (...args: A): void;
  /** Whether a call waits. */
  readonly pending: boolean;
  /** Runs the waiting call now. */
  flush(): void;
  /** Drops the waiting call. */
  cancel(): void;
  /** The view went: as `cancel`. */
  dispose(): void;
}

export function debounce<A extends unknown[]>(fn: (...args: A) => void, delay: number = EditDelay): Debounced<A> {
  let last: A | undefined;
  const timer = new Debouncer<void>(() => {
    const args = last as A;
    last = undefined;
    fn(...args);
  }, delay);
  const call = ((...args: A) => {
    last = args;
    timer.schedule();
  }) as Debounced<A>;
  Object.defineProperty(call, 'pending', { get: () => timer.waiting });
  call.flush = () => timer.flushAll();
  call.cancel = () => {
    last = undefined;
    timer.dispose();
  };
  call.dispose = call.cancel;
  return call;
}
