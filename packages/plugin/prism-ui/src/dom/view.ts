/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Child } from './jsx-types.js';
import { render } from './render.js';

/**
 * How the shell shows a contributed view (a tab, a section, an overlay, a command's custom UI): render into `host`,
 * return a disposer. Any technology can implement it; `domView` and `elementView` are the built-in ones, `reactView`
 * (in `react/`) renders with React or Preact.
 */
export interface ViewRenderer<Ctx> {
  mount(host: HTMLElement, ctx: Ctx): () => void;
}

/** A view written with the built-in JSX runtime. */
export function domView<Ctx>(component: (props: { ctx: Ctx }) => Child): ViewRenderer<Ctx> {
  return { mount: (host, ctx) => render(() => component({ ctx }), host) };
}

/** A view that is a custom element: the shell creates it and sets the context on `property` (default `ctx`). */
export function elementView<Ctx>(tag: string, property = 'ctx'): ViewRenderer<Ctx> {
  return {
    mount(host, ctx) {
      const el = document.createElement(tag);
      (el as any)[property] = ctx;
      host.appendChild(el);
      return () => el.remove();
    },
  };
}
