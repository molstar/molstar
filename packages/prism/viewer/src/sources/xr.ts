/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { combineLatest, of, type Subscription, switchMap } from 'rxjs';

// WebXR as the page reads it (state.md §3; parity.md's "illumination, XR" row): whether the canvas can enter XR and
// whether it is in it, for the toolbar's XR button (shown only where it can). Follows the canvas once it exists.

export interface XrState {
  readonly supported: boolean;
  readonly presenting: boolean;
}

export interface XrSource {
  readonly state: ReadonlyState<XrState>;
  dispose(): void;
}

export function createXrSource(plugin: PluginContext): XrSource {
  const state = observableState<XrState>({ supported: false, presenting: false });
  const sub: Subscription = plugin.behaviors.canvas3d.initialized
    .pipe(
      switchMap((ready) => {
        const xr = ready ? plugin.canvas3d?.xr : undefined;
        return xr ? combineLatest([xr.isSupported, xr.isPresenting]) : of([false, false] as const);
      }),
    )
    .subscribe(([supported, presenting]) => {
      const current = state.value;
      if (current.supported !== supported || current.presenting !== presenting) state.set({ supported, presenting });
    });
  return { state: state.asReadonly, dispose: () => sub.unsubscribe() };
}
