/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// `use(promise)` of React 19, for the catalog's React 18: a component reads a promise's value and suspends until it
// has one. The settled state is kept per promise, so a component may call it again with the same promise on every
// render (the stories' document builders cache theirs). No React import: Suspense catches a thrown promise.

type Settled<T> = { status: 'pending' } | { status: 'done'; value: T } | { status: 'failed'; reason: unknown };

const settled = new WeakMap<Promise<unknown>, Settled<unknown>>();

/** The value of `promise`: suspends the calling component while it is pending, throws its rejection. */
export function readPromise<T>(promise: Promise<T>): T {
  let s = settled.get(promise) as Settled<T> | undefined;
  if (!s) {
    s = { status: 'pending' };
    settled.set(promise, s);
    promise.then(
      (value) => settled.set(promise, { status: 'done', value }),
      (reason) => settled.set(promise, { status: 'failed', reason }),
    );
  }
  if (s.status === 'done') return s.value;
  if (s.status === 'failed') throw s.reason;
  throw promise;
}
