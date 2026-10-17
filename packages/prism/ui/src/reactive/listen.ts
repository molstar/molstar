/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { isObservable, type Observable } from 'rxjs';
import type { ReadonlyState } from './readonly-state.js';

/** Anything that tells of its changes with a plain callback (a state; a keyed selector's key). */
export interface Listenable {
  listen(fn: () => void): () => void;
}

/**
 * Calls `fn` after `source` changes, until the returned function is called. A state's own `listen` is used (a plain
 * callback); any other stream is subscribed, and what it emits while subscribing (a replayed current value) is not a
 * change.
 */
export function listen(source: Observable<unknown> | ReadonlyState<unknown> | Listenable, fn: () => void): () => void {
  const own = (source as ReadonlyState<unknown>).listen;
  if (typeof own === 'function') return own.call(source, fn);
  if (!isObservable(source)) throw new Error('listen: not a state or a stream');
  let subscribing = true;
  const sub = source.subscribe(() => {
    if (!subscribing) fn();
  });
  subscribing = false;
  return () => sub.unsubscribe();
}
