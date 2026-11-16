/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
import type { ReadonlyState, Subscribable } from '@molstar/core/reactive';
import { type ReactiveOptions, useReactiveSelection } from './use-reactive-selection.js';

function identity<T>(value: T): T {
  return value;
}

/**
 * A source's value in a React component (react §1): a state's current value (rendering again when it changes), a
 * stream's latest (`defaultValue` until its first). Subscribing never replays a state's value into React.
 */
export function useReactive<T>(source: ReadonlyState<T>, options?: ReactiveOptions<T>): T;
export function useReactive<T>(
  source: Subscribable<T> | null | undefined,
  options: ReactiveOptions<T> & { readonly defaultValue: T },
): T;
export function useReactive<T>(source: Subscribable<T> | null | undefined, options?: ReactiveOptions<T>): T | undefined;
export function useReactive<T>(
  source: Subscribable<T> | null | undefined,
  options?: ReactiveOptions<T>,
): T | undefined {
  return useReactiveSelection(source, identity<T>, options);
}
