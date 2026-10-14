/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Child } from './jsx-types.js';
import { render } from './render.js';
import type { Scope } from './scope.js';

/**
 * How the shell shows a contributed view (a tab, a section, an overlay, a command's custom UI): render into `host`,
 * return a disposer. Any technology can implement it; `domView` and `elementView` are the built-in ones, `reactView`
 * (in `react/`) renders with React or Preact.
 *
 * Hidden UI does no UI updates (styling.md §5): `scope` is the scope the view is mounted in, which the shell pauses
 * while the view is hidden (a folded section, a panel behind a card; the view stays mounted). `domView` follows it by
 * itself (its bindings unsubscribe and resume). A view that does its own subscribing, timers or observers pauses them
 * with `scope.onActive(active => …)`, which also tells it the state it starts in; it is not called after the disposer.
 */
export interface ViewRenderer<Ctx> {
  mount(host: HTMLElement, ctx: Ctx, scope: Scope): () => void;
}

/** A view written with the built-in JSX runtime. */
export function domView<Ctx>(component: (props: { ctx: Ctx }) => Child): ViewRenderer<Ctx> {
  return { mount: (host, ctx, scope) => render(() => component({ ctx }), host, scope) };
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
