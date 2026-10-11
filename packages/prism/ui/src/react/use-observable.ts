/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from 'rxjs';

import type { ReadonlyState } from '../reactive/readonly-state.js';
import { type ObservableOptions, useObservableSelection } from './use-observable-selection.js';

export type { ObservableOptions } from './use-observable-selection.js';

function identity<T>(value: T): T {
  return value;
}

export function useObservable<T>(source: ReadonlyState<T>, options?: ObservableOptions<T>): T;
export function useObservable<T>(
  source: Observable<T> | null | undefined,
  options: ObservableOptions<T> & { readonly defaultValue: T },
): T;
export function useObservable<T>(
  source: Observable<T> | null | undefined,
  options?: ObservableOptions<T>,
): T | undefined;
export function useObservable<T>(
  source: Observable<T> | null | undefined,
  options?: ObservableOptions<T>,
): T | undefined {
  return useObservableSelection(source, identity<T>, options);
}
