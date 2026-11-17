/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */
import type { ReactiveModel } from '@molstar/core/reactive';
import { useEffect } from 'react';

/**
 * Mounts a model while the component is mounted (`mount` / `unmount`; React's Strict Mode remount mounts it again),
 * never retires it: retiring is its creator's (react §3).
 */
export function useReactiveModel<T extends ReactiveModel>(model: T): T;
export function useReactiveModel<T extends ReactiveModel>(model: T | null | undefined): T | undefined;
export function useReactiveModel<T extends ReactiveModel>(model: T | null | undefined): T | undefined {
  useEffect(() => {
    if (!model) return;
    model.mount();
    return () => model.unmount();
  }, [model]);
  return model ?? undefined;
}
