/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from 'rxjs';

/** A field's value as a control sees it: a stream, the current value, and how to change it. */
export interface FieldBinding {
  value: Observable<any>;
  get(): any;
  set(value: any): void;
}
