/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '@molstar/prism-ui';
import { derive } from '@molstar/prism-ui/reactive/state';

/**
 * `source`'s value while `keep` holds for it, else the last one that did (what a `filter` did on a stream): a row's
 * record by key, gone only as the row goes; the card's turn while it is the latest. `source` holds a kept value when
 * this is made.
 */
export function held<S, T extends S>(source: ReadonlyState<S>, keep: (v: S) => v is T): ReadonlyState<T> {
  let last = source.value as T;
  return derive(source, (v) => (keep(v) ? (last = v) : last));
}

export const present = <T>(x: T | undefined): x is T => x !== undefined;
