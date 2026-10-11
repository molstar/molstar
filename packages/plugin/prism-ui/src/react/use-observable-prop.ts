/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Observable } from 'rxjs';

import type { ReadonlyState } from '../reactive/readonly-state.js';
import { type ObservableOptions, useObservableSelection } from './use-observable-selection.js';

export function useObservableProp<T, U>(
  source: ReadonlyState<T>,
  selector: (value: T) => U,
  options?: ObservableOptions<T, U>,
): U;
export function useObservableProp<T, U>(
  source: Observable<T> | null | undefined,
  selector: (value: T) => U,
  options: ObservableOptions<T, U> & { readonly defaultValue: T },
): U;
export function useObservableProp<T, U>(
  source: Observable<T> | null | undefined,
  selector: (value: T) => U,
  options?: ObservableOptions<T, U>,
): U | undefined;
export function useObservableProp<T, U>(
  source: Observable<T> | null | undefined,
  selector: (value: T) => U,
  options?: ObservableOptions<T, U>,
): U | undefined {
  return useObservableSelection(source, selector, options);
}
