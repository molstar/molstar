/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
 * A view written with React or Preact, mounted where the UI mounts views (react §4). The caller passes `createRoot`
 * (from `react-dom/client`, or `preact/compat/client`), so this package depends on neither.
 *
 * Hidden UI does no updates: the component gets `active` (false while the view is hidden but kept mounted) and renders
 * again when it changes; it can stop following its sources meanwhile (`useReactive(active ? model.x : undefined)`).
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
        container.remove();
        // After the current task: a view disposed from inside a React render or commit (a React app hosting a UI
        // element whose view is this) must not unmount a root synchronously there (react §4). Detached first, so it is
        // gone from the page at once.
        queueMicrotask(() => root.unmount());
      };
    },
  };
}
