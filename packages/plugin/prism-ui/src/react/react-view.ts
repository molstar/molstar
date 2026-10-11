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
 * A contributed view written with React or Preact. The caller passes `createRoot` (from `react-dom/client`, or
 * `preact/compat/client`), so prism-ui itself depends on neither.
 */
export function reactView<Ctx>(
  component: ComponentType<{ ctx: Ctx }>,
  options: { createRoot: (container: Element) => ReactRoot },
): ViewRenderer<Ctx> {
  return {
    mount(host, ctx) {
      const container = host.ownerDocument.createElement('div');
      container.style.display = 'contents';
      host.appendChild(container);
      const root = options.createRoot(container);
      root.render(createElement(component, { ctx }));
      return () => {
        root.unmount();
        container.remove();
      };
    },
  };
}
