/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from 'rxjs';
import type { ReadonlyState } from './readonly-state.js';

/**
 * Calls `fn` after `source` changes, until the returned function is called. A state's own `listen` is used (a plain
 * callback); any other stream is subscribed, and what it emits while subscribing (a replayed current value) is not a
 * change.
 */
export function listen(source: Observable<unknown> | ReadonlyState<unknown>, fn: () => void): () => void {
  const own = (source as ReadonlyState<unknown>).listen;
  if (typeof own === 'function') return own.call(source, fn);
  let subscribing = true;
  const sub = source.subscribe(() => {
    if (!subscribing) fn();
  });
  subscribing = false;
  return () => sub.unsubscribe();
}
