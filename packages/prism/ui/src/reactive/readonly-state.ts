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
   * Calls `fn` with each new value until the returned function is called: a plain callback, not an rxjs subscription
   * (rendering.md §2.2). States (`observableState`, `State`) and derived values (`derive`, `mapState`) have it;
   * `listen` (reactive/listen.ts) takes any stream.
   */
  readonly listen?: (fn: (value: T) => void) => () => void;
};
