/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from './readonly-state.js';
import { sameShallow } from './same-shallow.js';
import { derive } from './state.js';

/**
 * A state derived from another (prism-viewer README decision 22): `derive` with a shallow comparison by default.
 * Nothing is held while nothing listens.
 */
export function mapState<S, T>(
  source: ReadonlyState<S>,
  fn: (value: S) => T,
  compare: (a: T, b: T) => boolean = sameShallow,
): ReadonlyState<T> {
  return derive(source, fn, compare);
}
