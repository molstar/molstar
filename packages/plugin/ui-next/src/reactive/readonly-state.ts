/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from 'rxjs';

/** State with a current value replayed to new subscribers. */
export type ReadonlyState<T> = Observable<T> & {
  readonly value: T;
};
