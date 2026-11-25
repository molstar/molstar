/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
import type { ReadonlyState, Subscribable } from '@molstar/core/reactive';
import { type ReactiveOptions, useReactiveSelection } from './use-reactive-selection.js';

/**
 * A part of a source's value (`selector`), rendering again only when that part changes (`compare`, shallow by
 * default): a component that shows one field of a record state is not rendered for the others (react §2).
 */
export function useReactiveProp<T, U>(
  source: ReadonlyState<T>,
  selector: (value: T) => U,
  options?: ReactiveOptions<T, U>,
): U;
export function useReactiveProp<T, U>(
  source: Subscribable<T> | null | undefined,
  selector: (value: T) => U,
  options: ReactiveOptions<T, U> & { readonly defaultValue: T },
): U;
export function useReactiveProp<T, U>(
  source: Subscribable<T> | null | undefined,
  selector: (value: T) => U,
  options?: ReactiveOptions<T, U>,
): U | undefined;
export function useReactiveProp<T, U>(
  source: Subscribable<T> | null | undefined,
  selector: (value: T) => U,
  options?: ReactiveOptions<T, U>,
): U | undefined {
  return useReactiveSelection(source, selector, options);
}
