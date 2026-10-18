/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from '../../reactive/readonly-state.js';

/** A field's value as a control sees it: a state of it, the current value, and how to change it. */
export interface FieldBinding {
  value: ReadonlyState<any>;
  get(): any;
  set(value: any): void;
}
