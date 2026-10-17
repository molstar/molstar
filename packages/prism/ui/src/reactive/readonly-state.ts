/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from 'rxjs';

/** State with a current value replayed to new subscribers. */
export type ReadonlyState<T> = Observable<T> & {
  readonly value: T;
  /**
   * Calls `fn` after the state may have changed (no value: read `value`), until the returned function is called. A
   * plain callback, not an rxjs subscription: what `watch` listens with (rendering.md §2). States made by
   * `observableState` and `mapState` have it; `listen` (reactive/listen.ts) takes any stream.
   */
  readonly listen?: (fn: () => void) => () => void;
};
