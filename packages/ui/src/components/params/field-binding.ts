/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive, type ReadonlyState } from '@molstar/core/reactive';

/**
 * A field's value as a control sees it: a state of it (narrow: told only when this field changes; `.value` is the
 * current value) and how to change it.
 */
export interface FieldBinding<T = any> {
  readonly value: ReadonlyState<T>;
  set(value: T): void;
}

/** One key of a record field (a group's param, a mapped value's `params`), as a binding of its own. */
export function keyBinding(parent: FieldBinding, key: string): FieldBinding {
  return {
    value: derive(parent.value, (v) => v?.[key]),
    set: (x) => parent.set({ ...parent.value.value, [key]: x }),
  };
}
