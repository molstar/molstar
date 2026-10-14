/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ComponentType, createElement, type ReactNode } from 'react';
import type { ViewRenderer } from '../dom/view.js';

/** A React (or `preact/compat`) root, as `createRoot` from `react-dom/client` makes it. */
export interface ReactRoot {
  render(node: ReactNode): void;
  unmount(): void;
}

/**
 * A contributed view written with React or Preact. The caller passes `createRoot` (from `react-dom/client`, or
 * `preact/compat/client`), so prism-ui itself depends on neither.
 *
 * Hidden UI does no UI updates: the component gets `active` (false while the shell has the view hidden, a folded
 * section, a panel behind a card; true otherwise) and renders again when it changes. A component that subscribes to
 * models can skip that (`useObservable(active ? model.x : undefined)`) or stop its timers; React keeps it mounted.
 */
export function reactView<Ctx>(
  component: ComponentType<{ ctx: Ctx; active?: boolean }>,
  options: { createRoot: (container: Element) => ReactRoot },
): ViewRenderer<Ctx> {
  return {
    mount(host, ctx, scope) {
      const container = host.ownerDocument.createElement('div');
      container.style.display = 'contents';
      host.appendChild(container);
      const root = options.createRoot(container);
      let unmounted = false;
      // Called now with the state it starts in (the first render), then on each change.
      scope.onActive((active) => {
        if (!unmounted) root.render(createElement(component, { ctx, active }));
      });
      return () => {
        unmounted = true;
        root.unmount();
        container.remove();
      };
    },
  };
}
