/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { distinctUntilChanged, map } from 'rxjs';
import { listen } from './listen.js';
import type { ReadonlyState } from './readonly-state.js';
import { sameShallow } from './same-shallow.js';

/**
 * A state derived from another (prism-viewer README decision 22): its `value` is computed from the source's on read,
 * and it emits when the derived value changes (`compare`, shallow by default). Nothing is held: no subscription until
 * someone subscribes. Its `listen` is the source's.
 */
export function mapState<S, T>(
  source: ReadonlyState<S>,
  fn: (value: S) => T,
  compare: (a: T, b: T) => boolean = sameShallow,
): ReadonlyState<T> {
  const out = source.pipe(map(fn), distinctUntilChanged(compare)) as ReadonlyState<T>;
  Object.defineProperty(out, 'value', { get: () => fn(source.value), enumerable: true });
  // A listener is told when the source changes (the derived value may be the same: `listen` says "may have changed").
  Object.defineProperty(out, 'listen', { value: (f: () => void) => listen(source, f) });
  return out;
}
