/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { useEffect } from 'react';

import type { ReactiveModel } from '../reactive/reactive-model.js';

export function useReactiveModel<T extends ReactiveModel>(model: T): T;
export function useReactiveModel<T extends ReactiveModel>(model: T | null | undefined): T | undefined;
export function useReactiveModel<T extends ReactiveModel>(model: T | null | undefined): T | undefined {
  // React owns a repeatable mount scope, never the creator's final model retirement.
  useEffect(() => {
    if (!model) {
      return;
    }
    model.mount();
    return () => model.unmount();
  }, [model]);
  return model ?? undefined;
}
