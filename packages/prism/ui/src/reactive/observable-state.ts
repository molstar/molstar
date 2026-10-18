/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ReadonlyState } from './readonly-state.js';
import { isPlainRecord, sameShallow } from './same-shallow.js';
import { RecordState, State } from './state.js';

export { batch } from './state.js';

export type ObservableStateOptions<T> = {
  readonly compare?: (left: T, right: T) => boolean;
};

/** A writable state (`State`, reactive/state.ts): its value, its readers' view, `set`. */
export type ObservableState<T> = {
  readonly value: T;
  readonly asReadonly: ReadonlyState<T>;
  set(value: T): void;
};

export type ObservableRecordState<T extends Record<string, unknown>> = ObservableState<T> & {
  update(partial: Partial<T>): void;
};

/**
 * A writable state (rendering.md §2.2): the latest value plus its listeners, no rxjs on its path. `compare` decides
 * what is no change (shallow by default). A record's state has `update`.
 */
export function observableState<T extends Record<string, unknown>>(
  initial: T,
  options?: ObservableStateOptions<T>,
): ObservableRecordState<T>;
export function observableState<T>(initial: T, options?: ObservableStateOptions<T>): ObservableState<T>;
export function observableState<T>(initial: T, options?: ObservableStateOptions<T>): ObservableState<T> {
  const compare = options?.compare ?? sameShallow;
  return isPlainRecord(initial)
    ? (new RecordState(initial as Record<string, unknown>, compare as never) as unknown as ObservableState<T>)
    : new State(initial, compare);
}
